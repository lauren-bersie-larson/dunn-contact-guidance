function [network, magF, residuals, diff, forces, conv_flag, geo_flag] = solve_fixed_surf(intfib_nodes, orig_node_coords, network, modulus, fiber_area, fiber_B, L, R, direction, surface_nodes, bc_flag) 
% DESCRIPTION                                                                             
% ===========
% Solve for the position of interior nodes after a deformation
%
% INPUTS
% =========
% intfib_nodes = the numbers of the nodes on fibers within the cylinder to
% push out
%
% orig_node_coords = starting xyz coordinates of the nodes 
%
% network = structure containing network information such as node list,
% fiber list, boundary node numbers, interior node numbers, and initial
% fiber lengths
%
% modulus = fiber E value in code. This value * fiber_area equals fiber A
% value provided in paper, based on Lai et al. 2012 J. Biomech. Eng.
%
% fiber_area = fiber cross-sectional area
%
% fiber_B = fiber B parameter (nonlinearity parameter). Based on 
% Lai et al. 2012 J. Biomech. Eng.
%
% L = length of pseudopod cylinder (computational units)
%
% R = radius of pseudopod cylinder (computational units)
%
% direction = direction pseudopod is extended into (1 = x, 2 = y, etc.)



% COLLECT NETWORK INFO
% ====================
nodes = network.nodes;
fibers = network.fibers;
bnd_node_nums = network.bnd_node_nums;
int_node_nums = network.int_node_nums;
init_lens = network.init_lens;
[modulii] = modulus * ones( 1, length(init_lens) );

% POPULATE INPUT VECTORS
% ======================
% Only want to solve for positions of interior nodes- the boundary nodes
% are fixed. These will be known as the "free nodes"
free_node_x_indices = int_node_nums * 3 - 2 ;
free_node_y_indices = int_node_nums * 3 - 1 ;
free_node_z_indices = int_node_nums * 3 - 0 ;

if bc_flag == 1 %  For initial pseudopod extension- Do nothing

elseif bc_flag == 2 % For pseudopod retraction
    % Fix nodes on cylinder surface while retracting the cylinder by 
    % removing the surface node DOFs from the free_node_indices arrays    
    surface_nodes_x = surface_nodes*3 - 2;
    surface_nodes_y = surface_nodes*3 - 1;
    surface_nodes_z = surface_nodes*3;
    
    free_node_x_indices = setdiff(free_node_x_indices, surface_nodes_x);
    free_node_y_indices = setdiff(free_node_y_indices, surface_nodes_y);
    free_node_z_indices = setdiff(free_node_z_indices, surface_nodes_z);
end

free_node_indices = sort( [free_node_x_indices free_node_y_indices free_node_z_indices] );

X0 = nodes( free_node_indices );

% SOLVE USING EXPLICIT METHOD
% Find initial forces and the error
[forces,~] = calc_forces(nodes, fibers, init_lens, modulii, fiber_area, fiber_B, L, R, direction, int_node_nums);

INITIAL_ERR_VAL = norm(forces(free_node_indices))

[dXdt] = dampedFiber(1, nodes(free_node_indices), free_node_indices,...
        nodes, fibers, init_lens, modulii, fiber_area, fiber_B, L, R, direction, intfib_nodes, orig_node_coords, int_node_nums);

 nodes0=nodes;

 opts = '';
 [t, X0] = ode15s(@(t, X2)dampedFiber(t, X2, free_node_indices,...
        nodes, fibers, init_lens, modulii, fiber_area, fiber_B, L, R, direction, intfib_nodes, orig_node_coords, int_node_nums), [0 0.5], nodes(free_node_indices), opts);

% OPTIONAL: plot change in position of nodes during solve    
% figure;
% plot(t, (X0-nodes0(free_node_indices)));


 % OPTIONAL: Plot network as it displaces over time
%   if flag == 100 
%       videoName = VideoWriter('cyl_push_del.mp4', 'MPEG-4'); % Specify the output file name
%       videoName.Quality = 100;
%       open(videoName) 
% 
%      for iter = 1:4:length(t)
%          % Plot network
%          temp_nodes = nodes;
%          temp_nodes(free_node_indices) = X0(iter, :);
%          plot_network(temp_nodes, fibers, direction, h, k, l, R, L) %, network.bnd_node_nums, 2);
% 
%          % Set view of network
%          switch direction
%              case 1
%                  view(90, 0)
%              case 2
%                  view(0, 0)
%              case 3
%                  view(0, 90)
%          end
% 
%          % Plot settings
%         axis ([-0.23 0.23 -0.5 0.5 -0.5 0.5]); 
%         axis off 
% 
%          % Save figure to gif
%          fig_handle = getframe(gcf);   
%          filename = 'cyl_push.gif'; % Specify the output file name
%          im{iter} = frame2im(fig_handle);
% 
%          writeVideo(videoName,im{iter})
%          close;
%      end
% 
%   end
     
% apply X0 to nodal coordinates
nodes(free_node_indices) = X0(end,:);
network.nodes = nodes;

% get updated nodal force components
[forces, fiber_forces, network.fiber_stretch, magF, diff] = calc_forces(nodes, fibers, init_lens, modulii, fiber_area, fiber_B, L, R, direction, int_node_nums);

err = norm(forces(free_node_indices))

residuals = forces(free_node_indices);

% CHECK CONVERGENCE CRITERIA TO SEE IF CONDITIONS ARE MET 
% ======================
% 1. Force equilibrium 
if err <= 1e-10   
    fprintf('Converged! \n');
    conv_flag = 'y';
else 
    fprintf('Convergence failed. \n');
    conv_flag = 'n';
end

% 2. Geometric constraint 
% Ensure that nodes are pushed out of cylinder within certain tolerance
% Only check interior nodes of network
int_nodes_x = nodes(intfib_nodes * 3 - 2);  
int_nodes_y = nodes(intfib_nodes * 3 - 1); 
int_nodes_z = nodes(intfib_nodes * 3); 

% Calculate radial position of all nodes
switch direction
    case 1
        r = sqrt(int_nodes_y.^2+int_nodes_z.^2);

    case 2
        r = sqrt(int_nodes_x.^2+int_nodes_z.^2);

    case 3
        r = sqrt(int_nodes_x.^2+int_nodes_y.^2);
end

tol = 1e-2; 
if isempty(find(r <= R-tol, 1))
    fprintf('Geometric constraints met! All nodes pushed out \n');
    geo_flag = 'y';
else
    fprintf('Failed to push out all nodes. \n');
    geo_flag = 'n';
end


end
