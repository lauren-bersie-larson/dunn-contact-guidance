function [ dXdt] = dampedFiber(t,int_nodes, int_nodes_idx,...
    nodes, fibers, init_lens, fib_mods, fib_areas, fib_bs, L, R, direction, intfib_nodes, orig_node_coords, int_node_nums)
% using drag equation to find new solution
%   force is the nodal force (Fx,Fy,Fz)
%   dxdt is the nodal positon derivative (dx/dt, dy/dt, dz/dt)
% Ryan Mahutga
% 11-05-18

nodes(int_nodes_idx) = int_nodes;

[forces,~] = calc_forces(nodes, fibers, init_lens, fib_mods, fib_areas, fib_bs, L, R, direction, int_node_nums);

C = 1e-14;  %12; %-9; %1e-8;

dXdt = forces(int_nodes_idx)'./C;
end

