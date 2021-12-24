function [forces, force, fiber_stretch, magF, diff] = calc_forces(nodes, fibers, init_lens, modulii, fiber_area, fiber_B, L, R, direction, interior_nodes)
% IN
% ==
% NODES         1-by-3N array of xyzxyz... nodal coordinates for N nodes
% FIBERS        1-by-2N array of ababab... fiber nodes for N fibers
% INIT_LENS     1-by-N array of initial lengths for N fibers
% MODULII       1-by-N array of fiber modulii
%
% OUT
% ===
% FORCES        1-by-3N array of xyzxyz... force components acting on N nodes

num_fibers = length( fibers ) / 2;
nodes_x = nodes(1:3:end);
nodes_y = nodes(2:3:end);
nodes_z = nodes(3:3:end);

forces = zeros( 1, length(nodes) );
force = zeros(1, num_fibers);
fiber_stretch = zeros(1,num_fibers);
  
for j = 1 : num_fibers
    
    node_1_num = fibers(j*2-1);
    node_2_num = fibers(j*2-0);
    
    x_span = nodes_x( node_2_num ) - nodes_x( node_1_num );
    y_span = nodes_y( node_2_num ) - nodes_y( node_1_num );
    z_span = nodes_z( node_2_num ) - nodes_z( node_1_num );
    
    fiber_current_length = sqrt( x_span^2 + y_span^2 + z_span^2 );
    
    % EXP FORCE RELATION F = E A ( exp(B GS) - 1 ) / B
    
    lambda_limit = 1.3;
    
    lambda = fiber_current_length / init_lens(j);
    fiber_stretch(j) = lambda;
    
    if lambda > lambda_limit
        
        gs = 0.5*(lambda_limit^2-1);
        force_exp = modulii(j) * fiber_area * ( exp(fiber_B*gs) -1) / fiber_B;
        slope_at_lambda_limit = modulii(j) * fiber_area * lambda_limit * exp(fiber_B*gs);
        force(j) = force_exp + slope_at_lambda_limit * (lambda - lambda_limit);
        
    else
        gs = 0.5 * ( lambda^2 - 1 );
        force(j) = modulii(j) * fiber_area * (( exp( fiber_B * gs) - 1 ) / fiber_B);
        
    end
    
    cosine_A = x_span / fiber_current_length;
    cosine_B = y_span / fiber_current_length;
    cosine_C = z_span / fiber_current_length;
    
    node_1_force_x = force(j) * cosine_A * (+1);
    node_1_force_y = force(j) * cosine_B * (+1);
    node_1_force_z = force(j) * cosine_C * (+1);
    
    node_2_force_x = force(j) * cosine_A * (-1);
    node_2_force_y = force(j) * cosine_B * (-1);
    node_2_force_z = force(j) * cosine_C * (-1);
    
    forces( node_1_num*3-2 ) = forces( node_1_num*3-2 ) + node_1_force_x;
    forces( node_1_num*3-1 ) = forces( node_1_num*3-1 ) + node_1_force_y;
    forces( node_1_num*3-0 ) = forces( node_1_num*3-0 ) + node_1_force_z;
    
    forces( node_2_num*3-2 ) = forces( node_2_num*3-2 ) + node_2_force_x;
    forces( node_2_num*3-1 ) = forces( node_2_num*3-1 ) + node_2_force_y;
    forces( node_2_num*3-0 ) = forces( node_2_num*3-0 ) + node_2_force_z;
    
end

magF = zeros(length(nodes_x), 1);
cylR = R.* ones(length(nodes_x), 1);
P = L-0.5; %define x-coord of protrusion end - between -0.5 and 0.5

switch direction
    case 1
        r = sqrt(nodes_y.^2+nodes_z.^2);
        
    case 2
        r = sqrt(nodes_x.^2+nodes_z.^2);
        
    case 3
        r = sqrt(nodes_x.^2+nodes_y.^2);
end
   
Ap = 1e-3; 

diff = cylR-r';
diff(diff < 0) = 0; % If nodes are outside cylinder, set diff = 0
[~, pushed_nodes] = find(diff'); % Nodes that need to be pushed will have a nonzero diff
pushed_nodes = pushed_nodes';
magF = Ap.*diff;
magF = magF';

% Calculate components for force- ensuring that it does not change
% sign because of any small perturbations!
x_comp = nodes_x./r;
y_comp = nodes_y./r;
z_comp = nodes_z./r;

switch direction
    case 1
        forces( pushed_nodes.*3-1 ) = forces( pushed_nodes.*3-1 ) +((magF(pushed_nodes) .*  (y_comp(pushed_nodes))));
        forces( pushed_nodes.*3-0 ) = forces( pushed_nodes.*3-0 ) + ((magF(pushed_nodes) .* (z_comp(pushed_nodes))));
                        
    case 2
        forces(pushed_nodes*3-2 ) = forces( pushed_nodes*3-2 ) + (magF(pushed_nodes) .*  (x_comp(pushed_nodes)));
        forces( pushed_nodes*3-0 ) = forces( pushed_nodes*3-0 ) + (magF(pushed_nodes) .* (z_comp(pushed_nodes)));
        
    case 3
        forces(pushed_nodes*3-2 ) = forces(pushed_nodes*3-2 ) + (magF(pushed_nodes) .* (x_comp(pushed_nodes)));
        forces(pushed_nodes*3-1 ) = forces(pushed_nodes*3-1 ) + (magF(pushed_nodes) .* (y_comp(pushed_nodes)));
end

end