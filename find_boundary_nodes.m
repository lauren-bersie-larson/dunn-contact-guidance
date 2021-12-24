function [bnd_node_nums] = find_boundary_nodes(nodes,boundaries)

% Out
% ===
% bnd_node_nums -- N x 1 col vector of N boundary node numbers


% SETUP VARIABLES

xmin = boundaries(1);
xmax = boundaries(2);
ymin = boundaries(3);
ymax = boundaries(4);
zmin = boundaries(5);
zmax = boundaries(6);

num_nodes = length(nodes) / 3;

nodes_x = nodes(1:3:end);
nodes_y = nodes(2:3:end);
nodes_z = nodes(3:3:end);

bnd_node_nums = []; % will populate with bnd node indices


% FIND BOUNDARY NODES GENERICALLY


for n = 1 : num_nodes
    
    if( nodes_x(n)==xmin || nodes_x(n)==xmax || nodes_y(n)==ymin || nodes_y(n) == ymax || nodes_z(n)==zmin || nodes_z(n) == zmax )
        bnd_node_nums = cat(2,bnd_node_nums, n); %[bnd_node_nums n];
        
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RVE Initial Boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     7-------8   node1 <-0.5,-0.5,0.5>
%    /!      /!   node2 <0.5,-0.5,0.5>
%   / !     / !   node3 <-0.5,-0.5,-0.5>
%  5--!----6  !   node4 <0.5,-0.5,-0.5>
%  !  3----!--4   node5 <-0.5,0.5,0.5>
%  ! /     ! /    node6 <0.5,0.5,0.5>
%  !/      !/     node7 <-0.5,0.5,-0.5>
%  1-------2      node8 <0.5,0.5,-0.5>

% ttrve.x(1)=xmin;
% ttrve.y(1)=ymin;
% ttrve.z(1)=zmax;
% 
% ttrve.x(2)=xmax;
% ttrve.y(2)=ymin;
% ttrve.z(2)=zmax;
% 
% ttrve.x(3)=xmin;
% ttrve.y(3)=ymin;
% ttrve.z(3)=zmin;
% 
% ttrve.x(4)=xmax;
% ttrve.y(4)=ymin;
% ttrve.z(4)=zmin;
% 
% ttrve.x(5)=xmin;
% ttrve.y(5)=ymax;
% ttrve.z(5)=zmax;
% 
% ttrve.x(6)=xmax;
% ttrve.y(6)=ymax;
% ttrve.z(6)=zmax;
% 
% ttrve.x(7)=xmin;
% ttrve.y(7)=ymax;
% ttrve.z(7)=zmin;
% 
% ttrve.x(8)=xmax;
% ttrve.y(8)=ymax;
% ttrve.z(8)=zmin;



end