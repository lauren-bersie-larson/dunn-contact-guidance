function [init_lens] = calc_lens(nodes,fibers)


% IN
% ==
% nodes - 1-by-3N vect of xyz for N nodes  
% fibers - 1-by-2N vect for node numbers 1-2 for N fibers
%
% OUT
% ===
% init_lens - 1 x N vect of N fiber lengths

num_fibers = length( fibers ) / 2;

nodes_x = nodes(1:3:end);
nodes_y = nodes(2:3:end);
nodes_z = nodes(3:3:end);

init_lens = zeros(num_fibers,1);

for n = 1 : num_fibers
    
    node_1_num = fibers(n*2-1);
    node_2_num = fibers(n*2-0);
    
    node_1_vect = [nodes_x(node_1_num) nodes_y(node_1_num) nodes_z(node_1_num)];
    node_2_vect = [nodes_x(node_2_num) nodes_y(node_2_num) nodes_z(node_2_num)];
    
    init_lens(n) = norm( node_1_vect - node_2_vect );

end



end