function plot_network(nodes, fibers, direction, h, k, l, R, L) %, boundary_nodes)

%--------------------------------------------------------------------------
%
% PLOT NETWORK
%
% nodes         N x 3 array of xyz coordinates for N nodes
% fibers        N x 2 array of start-end nodes for N fibers
%
%--------------------------------------------------------------------------



figure;

num_nodes = length(nodes)/3;

for i = 1:num_nodes      
    nodes_x = nodes(3*i-2);
    nodes_y = nodes(3*i-1);
    nodes_z = nodes(3*i);
    
%    if ismember(i, boundary_nodes)
%          plot3( nodes_x , nodes_y, nodes_z,'o', 'LineWidth',1, 'MarkerEdgeColor','r','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6);
%    else
        plot3( nodes_x , nodes_y, nodes_z,'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',6);
%    end

hold on;
end

plot3( nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',1, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',3);


num_fibers = length(fibers)/2;

nodes = cat(2,nodes(1:3:end)', nodes(2:3:end)', nodes(3:3:end)');
fibers = cat(2,fibers(1:2:end)', fibers(2:2:end)'); 

%keyboard
for n = 1 : num_fibers
    
    x(1) = nodes( fibers(n,1),1 ); % node 1 x coord
    y(1) = nodes( fibers(n,1),2 ); % node 1 y coord
    z(1) = nodes( fibers(n,1),3 ); % node 1 z coord
        
    x(2) = nodes( fibers(n,2),1 ); % node 2 x coord
    y(2) = nodes( fibers(n,2),2 ); % node 2 y coord
    z(2) = nodes( fibers(n,2),3 ); % node 2 z coord
     
        
    plot3( x, y, z, 'b' );
       
    hold on;
       
end


if 1==1
    A = meshgrid(linspace(0, 2*pi, 20), linspace(0, 2*pi, 20)) ;
    
    switch direction 
        case 1
            Z = R.* cos(A)+l;
            Y = R.* sin(A)+k;
            X = meshgrid(linspace(h, h+L, 20), linspace(h, h+L, 20))';
            
        case 2
            X = R.* cos(A)+h;
            Z = R.* sin(A)+l;
            Y = meshgrid(linspace(k, k+L, 20), linspace(k, k+L, 20))';
            
        case 3
            X = R.* cos(A)+h;
            Y = R.* sin(A)+k;
            Z = meshgrid(linspace(-0.5, -0.5+L, 20), linspace(-0.5, -0.5+L, 20))';
    end

surf(X, Y, Z,'AlphaData',0.2);

    
end


set(gcf, 'color', 'white');



axis equal; %axis off;

%view( 13, 12 );
view(-42,22);
%axis_val = max( nodes(1:end) ) + 0.1;     % find a max/min for plot
axis_val = 0.5; 
axis( [-axis_val axis_val -axis_val axis_val -axis_val axis_val] );
%axis([-0.1 0.1 -0.5 0.3 -0.1 0.1]);

box on;

xlabel('x')'; ylabel('y'); zlabel('z');

end
