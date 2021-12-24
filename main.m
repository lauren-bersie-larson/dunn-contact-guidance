%%% Main file for contact guidance code
% Written by Lauren Bersie-Larson and Victor Lai
clc;
close all;
clear all;

addpath('Delaunay Networks')
%% Network and cylinder geometry info
% RVE boundaries in computational units 
% format: [-x +x -y +y -z +z]
boundaries = [-0.5 +0.5 -0.5 +0.5 -0.5 +0.5];

% How many fiber section numbers (sect_no) I want to add to each fiber 
% intersecting the cylinder, aka the level of discretization
sect_no = 2; 

% Define material parameters for the fiber network
% Here, the product fiber_A_val * fiber_area = fiber A in our paper,
% based on the values from Lai et al. 2012 J. Biomech. Eng.
fiber_A_val = 1.0823e+07; 
fiber_B = 2.5; 
fiber_area = pi*((100e-9)^2); 

% Define parameters for the protrusion
direction = 2; % Defines the axis the cylinder will extend parallel to 
L = 0.5; % Cylinder length in computational units
R = 0.125; % Cylinder radius in computational units

if direction==1
    h = -0.25;
    k = 0;
    l = 0;
    prot_node1 = [h k l]; % coordinates of center of cylinder at one end
    prot_node2 = [h+L k l]; % coordinates of center of cylinder at other end
    rad1 = 2;  %rad2 and rad3 are used to define which coordinates to use to solve for the intersection between fiber and protrusion (s1 and s2 below)
    rad2 = 3;
end

if direction==2
    h = 0;
    k = -0.25;
    l = 0;
    prot_node1 = [h k l]; % coordinates of center of cylinder at one end
    prot_node2 = [h k+L l]; % coordinates of center of cylinder at other end
    rad1 = 1;
    rad2 = 3;
end

if direction==3
    h = 0;
    k = 0;
    l = -0.25;
    prot_node1 = [h k l]; % coordinates of center of cylinder at one end
    prot_node2 = [h k l+L]; % coordinates of center of cylinder at other end
    rad1 = 1;
    rad2 = 2;
end

addpath('Delaunay Networks'); % Directory where networks are saved

fibseg_nodes1 = [];
network_list = [1,22,43,64,85,106,127,148,169,190,211,232,253,274,295,316,337,358,379,400,421,442,463,484]; %[1:504]; 

%% Loop through networks to be solved
 for i = 1:length(network_list)
    close all
    tic
    
    % Load network
    n = network_list(i);
    network_file = ['Del_',num2str(n),'.mat'];
    load(network_file);

    % Calculate network orientation before protrusion
    omega_bef = calc_orient(nodes, fibers)   
    
    % Both fibers and nodes are in 2D - convert to linear
    fibers = conv_2D_2_lin(fibers);   
    nodes = conv_2D_2_lin(nodes);

    % Round nodal coordinates! When reading a .mat file, boundary nodes may 
    % not get recognized due to round-off error
    nodes = round(nodes, 8);

    bnd_nodenum = find_boundary_nodes(nodes,boundaries);

    % Running list to track boundary nodes within pseudopod 
    rem_bndnode = [];  

    for numnode = 1:length(bnd_nodenum)
        bndnode = [nodes(3*bnd_nodenum(numnode) - 2) nodes(3*bnd_nodenum(numnode) - 1) nodes(3*bnd_nodenum(numnode))];
        r = sqrt(bndnode(rad1)^2 + bndnode(rad2)^2);

        if r <= R
            rem_bndnode = [rem_bndnode bnd_nodenum(numnode)];
        end
    end
    
    % Find boundary nodes
    [network.bnd_node_nums] = find_boundary_nodes(nodes,boundaries);
    % Find interior nodes
    [network.int_node_nums] = find_int_nodes(nodes, boundaries); 
    % Calculate fiber lengths
    [network.init_lens] = calc_lens(nodes,fibers);

    network.nodes = nodes;
    network.fibers = fibers;

    plot_network(network.nodes, network.fibers, direction, h, k, l, R, L); 
    figure_file1 = ['Network_before_',num2str(n),'.fig'];
    savefig(figure_file1);

    num_fibers = length(fibers)/2;
        
    % Discretize fibers inside the cylinder. These discretized fibers 
    % are  listed in the variable intfib_index.  
    [nodes, fibers, intfib_index, num_fibers, cyl_bnd_nodes] = discretize_fibs(num_fibers, fibers, nodes, sect_no, direction, R, prot_node1, prot_node2, rad1, rad2);
    
    % Calculate initial lengths of fibers before they get pushed out
    initlens = calc_lens(nodes,fibers);
    intfib_initlens = zeros(size(intfib_index,1),1);
    intfib_nodes = [];
    fibseg_nodes = [];
    
    % Determine the node numbers on the discretized fibers that will need
    % to be pushed out. These node numbers will be called intfib_nodes.
    for j = 1 : size(intfib_index,1)
        fib_segment_nos = intfib_index(j,:);
        fib_segment_nos(isnan(fib_segment_nos)) = [];

        for p = 1:length(fib_segment_nos)
            fibseg_nodes = [fibseg_nodes fibers(2*fib_segment_nos(p)-1) fibers(2*fib_segment_nos(p))];
        end

        % Track all nodes of intersecting fibers
        intfib_nodes = [intfib_nodes fibseg_nodes];  
        seg_lens = calc_lens(nodes,fibseg_nodes);
        intfib_initlens(j) = sum(seg_lens);
    end

    intfib_nodes = unique(intfib_nodes); % Remove any repeated nodes

    % Update nodes and fibers lists for calculating forces and stress
    orig_node_coords = nodes; % Will use later for applying tethering force!
    network.nodes = nodes;
    network.fibers = fibers;
    [network.init_lens] = calc_lens(nodes,fibers);
    [network.int_node_nums] = find_int_nodes(nodes, boundaries);

    % Plot network after adding new nodes and save
    plot_network(network.nodes, network.fibers, direction, h, k, l, R, L);
    figure_file = ['Network_during_',num2str(n),'.fig'];
    savefig(figure_file);       

    % Solve
    [network, magF, residuals, diff, forces, conv_flag, geo_flag] = solve_fixed_surf(intfib_nodes, orig_node_coords, network, fiber_A_val, fiber_area, fiber_B, L, R, direction, [], 1);
    plot_network_final(network.nodes, network.fibers, intfib_nodes, direction, h, k, l, R, L);

    err = norm(residuals);
    
    % Plot final network and save
    plot_network(network.nodes, network.fibers, direction, h, k, l, R, L);
    figure_file = ['Network_after_',num2str(n),'.fig'];
    savefig(figure_file);
    
    % Calculate length of fibers in contact with cylinder
    [contact_len2, surface_nodes] = calc_contact_len(network.nodes, network.fibers, direction, intfib_nodes, h, k, l, R, L, initlens, prot_node1, prot_node2);
    fprintf('contact length = %f \n', contact_len2);
    
    %sum up all protrusion forces on all nodes
    prot_force = 0;

    % Sum up all protrusion forces on all nodes for total pushing
    % force
    for p = 1:length(network.nodes)/3
        prot_force = prot_force + magF(p); 
    end
    fprintf(1,'protrusion force = %e \n' , prot_force);

    % Calculate network orientation after protrusion
    omega_aft = calc_orient(network.nodes, network.fibers);  

    % Retract protrusion and calculate local stiffness sensed
    [loc_stiffness, tot_cyl_forces] = retract_cyl(network, initlens, fiber_A_val, fiber_area, fiber_B, direction, h, k, l, R, L, surface_nodes);

    tot_time = toc
    
    % Write results to a text file
    fnm = sprintf('%s%s', 'Results', '.txt');
    filename = fnm;
    fileid = fopen(filename , 'a'); % appends any existing file

    fprintf(fileid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n', ...  % file line filters
    'Network #', ...
    'Direction', ...
    'Protrusion Force', ...
    'Contact Length',...
    'Local Stiffness',...
    'Convergence',...
    'Geo Flag',...
    'Final Error',...
    'Omega xx before',...
    'Omega yy before',...
    'Omega zz before',...
    'Omega xx after',...
    'Omega yy after',... 
    'Omega zz after', ...
    'Total simulation time');

    fprintf(fileid, '%i %i %e %e %e %s %s %e %f %f %f %f %f %f %f \n', ... 
    n, ...
    direction, ...
    prot_force, ...
    contact_len2, ...
    loc_stiffness, ...
    conv_flag, ...
    geo_flag, ...
    err,...
    omega_bef(1,1),...
    omega_bef(2,2), ...
    omega_bef(3,3), ...
    omega_aft(1,1), ...
    omega_aft(2,2), ...
    omega_aft(3,3), ...
    tot_time);

    fclose(fileid);
 end

