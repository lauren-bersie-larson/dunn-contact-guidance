function [loc_stiffness, tot_cyl_forces] = retract_cyl(network, init_lens, modulus, fiber_area, fiber_B, direction, h, k, l, R, L, surface_nodes) %, cyl_bnd_nodes)
% Written by Lauren Bersie-Larson, 4-13-19
% This function takes the solved network state, where nodes have been 
% pushed out of the cylinder, and then retracts the cylinder to measure the
% so-called local stiffness that the cylinder senses. 
% This is done in the following steps: 
% 1. Nodes that are within a certain tolerance from the cylinder radius are
% considered to be in contact with the cylinder surface, and are then 
% retracted with the cylinder by changing the nodal coordinates a certain
% amount (this assumes they are fixed to cylinder surface)
% 2. The fiber forces are then calculated based on the displacement
% 3. Nodal forces for nodes on the cylinder surface are calculated
% 4. Nodal forces are the opposite of the force the cylinder is applying
% 5. Displacement vs. cylinder force is tracked, and the slope is the
% "stiffness"

modulii = modulus * ones(1, length(init_lens));
num_disp_steps = 4; % Number of displacement steps to calculate forces
disp_step_size = 0.05; % Size of displacement step
tot_cyl_forces = [];
disp_steps = [disp_step_size:disp_step_size:(disp_step_size*num_disp_steps)];

fprintf('================================ \n');
fprintf('RETRACTING CYLINDER \n');

for step = 1:num_disp_steps
    fprintf('RETRACTION STEP: %i', step); 
    % Displace cylinder nodes
    switch direction
        case 1
            network.nodes(surface_nodes*3-2) = network.nodes(surface_nodes*3-2) - disp_step_size;
        case 2
            network.nodes(surface_nodes*3-1) = network.nodes(surface_nodes*3-1) - disp_step_size;
        case 3
            network.nodes(surface_nodes*3) = network.nodes(surface_nodes*3) - disp_step_size;
    end
    
    % Solve for equilibrated network 
    [network, ~, residuals, ~, forces, conv_flag, geo_flag] = solve_fixed_surf([], [], network, modulus, fiber_area, fiber_B, L, R, direction, surface_nodes, 2); %, h, k, l, 2, [], 2, surface_nodes, cyl_bnd_nodes);
    
    err = norm(residuals);
    % Print final error of retraction step
    fprintf('%s%i%s%e \n', 'Retraction step: ', step, ' Residual: ', err);

    % Only collect the nodal forces for the nodes on the cylinder surface
    xforces = forces(surface_nodes.*3-2);
    yforces= forces(surface_nodes.*3-1);
    zforces = forces(surface_nodes.*3);
    
    % Change the sign on the forces to make it the cylinder force
    cyl_forcesx = xforces .* -1;
    cyl_forcesy = yforces .* -1;
    cyl_forcesz = zforces .* -1;
        
    % Calculate and sum the forces on the cylinder nodes
    result_forces = [];
    switch direction
        case 1
                result_forces = cyl_forcesx;
        case 2
                result_forces = cyl_forcesy;
        case 3
                result_forces = cyl_forcesz;
    end
  
    tot_cyl_forces = [tot_cyl_forces; abs(sum(result_forces))];
    
    % Plot network with retracted cylinder and save
    fig_handle = plot_network_retract(network.nodes, network.fibers, network.fiber_stretch, direction, h, k, l, R, L, step, disp_step_size);
    figure_file = ['Network_retract_',num2str(step),'.fig'];
    im{step} = frame2im(fig_handle);
    fprintf('-------------------------------- \n');
end

% OPTIONAL: Write all retraction figures to a gif
% filename = 'cyl_retract.gif'; % Specify the output file name
% for idx = 1:num_disp_steps
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     end
% end

% OPTIONAL: Plot animated stiffness plot
% disp_steps = disp_steps * 10; % Dimensionalize from CU to um
% filename = 'stiffness.gif';
% p = figure;
% title('Local Stiffness')
% xlabel('Pseudopod Displacement (um)') 
% ylabel('Force on Pseudopod (N)') 
% xlim([0 disp_steps(end)])
% ylim([0 tot_cyl_forces(end)])
% set(gcf, 'color', 'white');
% set(gca,'FontSize',16)
% h = animatedline('LineWidth',2);
% for n = 1:length(disp_steps)
%     % Draw plot
%     x = disp_steps;
%     y = tot_cyl_forces;
%     addpoints(h, x(n),y(n));
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(p); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',1); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',1); 
%       end 
% end

% OPTIONAL: Plot static displacement plot
% figure;
% title('Local Stiffness')
% xlabel('Pseudopod Displacement (um)') 
% ylabel('Force on Pseudopod (N)') 
% xlim([0 disp_steps(end)])
% ylim([0 tot_cyl_forces(end)])
% set(gcf, 'color', 'white');
% set(gca,'FontSize',16)
% plot(disp_steps, tot_cyl_forces);
% savefig('Stiffness_plot.fig');    

% Calculate the local stiffness
loc_stiffness = (tot_cyl_forces(2) - tot_cyl_forces(1))/(disp_step_size);
fprintf('Local stiffness is: %e \n', loc_stiffness)

end