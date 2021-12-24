function [contact_len2, surface_nodes] = calc_contact_len(nodes, fibers, direction, fibseg_nodes, h, k, l, R, L, initlens, prot_node1, prot_node2)
    % Written by Lauren Bersie-Larson
    num_fibers = length(fibers)/2;
    nodes_x = nodes(1:3:end);
    nodes_y = nodes(2:3:end);
    nodes_z = nodes(3:3:end);
    
    % Calculate radius value for every node
    switch direction
        case 1
            radii = sqrt((nodes_y.^2) + (nodes_z.^2));
        case 2
            radii = sqrt((nodes_x.^2) + (nodes_z.^2));
        case 3
            radii = sqrt((nodes_x.^2) + (nodes_y.^2));
    end
    
    % Find nodes that are on cylinder surface
    % Because nodes may not be exactly on surface, there will be an
    % acceptable tolerance +/- from the cylinder radius that will be
    % considered on the surface
    tol = 1e-2; % Check to make sure this is good enough
    surface_nodes = find((radii <= R+tol) & (radii >= R-tol));
    % Ensure nodes are inside the cylinder length!
    nodes3d = conv_lin_2_2D(nodes, 3);
    surface_nodes = surface_nodes(find((nodes3d(surface_nodes, direction) <= prot_node2(direction)+tol) & (nodes3d(surface_nodes, direction) >= prot_node1(direction)-tol))); 
    
    % Find fibers with both nodes that are surface nodes. Fibers that only
    % have one surface node will not be considered to have a fiber length
    % that is in contact with the surface in any meaningful amount
    
    % Find fibers in contact with surface
    fibers = conv_lin_2_2D(fibers, 2);
    surface_fibs = [];
    
    for i = 1 : num_fibers
        if (ismember(fibers(i,1), surface_nodes) && ismember(fibers(i,2), surface_nodes))
            surface_fibs = [surface_fibs i]; % Fiber numbers
        end
    end
    
    fibers2 = fibers(surface_fibs, :);
    fibers2 = conv_2D_2_lin(fibers2);
    surface_lens = calc_lens(nodes, fibers2);
    
    % Only count fiber lengths that are below stretch of 1.4
    stretches = surface_lens ./ initlens(surface_fibs);
    contact_len2 = sum(surface_lens);
    
    % Plot to show which fibers are being counted  
    fibers = conv_2D_2_lin(fibers);
    plot_network_final_lens(nodes, fibers, fibseg_nodes, direction, h, k, l, R, L, surface_fibs)

end