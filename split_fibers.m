function [nodes, fibers, split_fiber_index, cyl_bnd_nodes] = split_fibers(nodes, fibers, fib_num, sectno) %sect_no)

%create new nodes along a single fiber by splitting the fiber into a
%defined number of equal-length segments.
%to create n number of segments, you create (n-1) nodes along the fiber

% Changed 2-23-2018 LMB
% Now instead of creating certain # of fib segments, creates as many as
% needed to make segments of certain length
% sectno is now sect_dist (the distance needed for 

start_node_num = fibers(2*fib_num-1);
end_node_num = fibers(2*fib_num);

fib_start_node = [nodes(3*start_node_num-2) nodes(3*start_node_num-1) nodes(3*start_node_num)];
fib_end_node = [nodes(3*end_node_num-2) nodes(3*end_node_num-1) nodes(3*end_node_num)];

dir_v = fib_end_node - fib_start_node; %directional vector of fiber to be split

fib_len = sqrt(sum((fib_end_node - fib_start_node).^2));

% sectno = floor(fib_len/sect_dist);


% Make sure sectno is even
if mod(sectno, 2) == 0
    % Is even - do nothing
else 
    sectno = sectno + 1; % Make even
end

%create a new variable to track all split fibers. Each row contains the
%connecting fiber numbers along this split fiber.
if sectno > 1
    split_fiber_index = [fib_num];

    for i = 1:sectno-1

        new_node(i,:) = fib_start_node + i*(dir_v/sectno);

        %update overall nodes list - add new node coordinates to the end of
        %nodes list
        nodes = [nodes new_node(i,:)];

        if i == 1
            %set original entire fiber to now become the first segment of split fiber
            fibers(2*fib_num-1) = start_node_num;
            fibers(2*fib_num) = length(nodes)/3;

        else
            %create new fiber, with start and end points corresponding to the last two created nodes
            fibers(end+1) = length(nodes)/3 - 1;
            fibers(end+1) = length(nodes)/3;
            split_fiber_index = [split_fiber_index length(fibers)/2];   %keep track of fiber numbers along this split fiber
        end

    end

    %create last fiber segment
    fibers(end+1) = length(nodes)/3;
    fibers(end+1) = end_node_num;
    split_fiber_index = [split_fiber_index length(fibers)/2];
    cyl_bnd_nodes = [start_node_num; end_node_num];
else 
    split_fiber_index = [];
    cyl_bnd_nodes = [];
end
end


