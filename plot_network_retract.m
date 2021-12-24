function [fig_handle] = plot_network_retract(nodes, fibers, fiber_stretch, direction, h, k, l, R, L, disp_step, disp_step_size)

%--------------------------------------------------------------------------
%
% PLOT NETWORK
%
% nodes         N x 3 array of xyz coordinates for N nodes
% fibers        N x 2 array of start-end nodes for N fibers
%
%--------------------------------------------------------------------------
disp = disp_step * disp_step_size;
fiber_size = 1.5;

figure;

num_nodes = length(nodes)/3;

for i = 1:num_nodes      
    nodes_x = nodes(3*i-2);
    nodes_y = nodes(3*i-1);
    nodes_z = nodes(3*i);
    
    plot3( nodes_x , nodes_y, nodes_z,'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6);

hold on;
end

num_fibers = length(fibers)/2;

nodes = cat(2,nodes(1:3:end)', nodes(2:3:end)', nodes(3:3:end)');
fibers = cat(2,fibers(1:2:end)', fibers(2:2:end)'); 

for n = 1 : num_fibers
    
    x(1) = nodes( fibers(n,1),1 ); % node 1 x coord
    y(1) = nodes( fibers(n,1),2 ); % node 1 y coord
    z(1) = nodes( fibers(n,1),3 ); % node 1 z coord
        
    x(2) = nodes( fibers(n,2),1 ); % node 2 x coord
    y(2) = nodes( fibers(n,2),2 ); % node 2 y coord
    z(2) = nodes( fibers(n,2),3 ); % node 2 z coord
    
    % Plot fiber
    p_fib = plot3(x, y, z, 'Color', 'b');
      
    hold on;
       
end


if 1==1
    A = meshgrid(linspace(0, 2*pi, 20), linspace(0, 2*pi, 20)) ;
    
    switch direction 
        case 1
            Z = R.* cos(A)+l;
            Y = R.* sin(A)+k;
            X = meshgrid(linspace(h-disp, h+L-disp, 20), linspace(h-disp, h+L-disp, 20))';
            
        case 2
            X = R.* cos(A)+h;
            Z = R.* sin(A)+l;
            Y = meshgrid(linspace(k-disp, k+L-disp, 20), linspace(k-disp, k+L-disp, 20))';
            
        case 3
            X = R.* cos(A)+h;
            Y = R.* sin(A)+k;
            Z = meshgrid(linspace(-0.5-disp, -0.5+L-disp, 20), linspace(-0.5-disp, -0.5+L-disp, 20))';
    end
    s = surf(X, Y, Z, 'FaceColor', 'green');
    set(s,'LineWidth', 0.6)
end

set(gcf, 'color', 'white');


axis equal; 

view(-42,22);

switch direction
    case 1
        axis([-1.5 1 -0.25 1 -1 1]);
        view(0, 0)
    case 2
        axis([-1 1 -1.5 1 -1 1]);
        view(90, 0)
    case 3
end

axis off
fig_handle = getframe(gcf);
end
