function [nodes, fibers, intfib_index, num_fibers, cyl_bnd_nodes] = discretize_fibs(num_fibers, fibers, nodes, sect_no, direction, R, prot_node1, prot_node2, rad1, rad2)
    %%% Finds all possible ways a ray (fiber) can intersect a finite cylinder. 
    %%% Once found, fibers that are within the cylinder are discretized. 

    %%% Function originally written by Victor Lai
    %%% Updated to include all cases and generalize for a finite length 
    %%% cylinder by Lauren Bersie-Larson.
   
    intfib_index = []; %keeps track of fiber segment numbers of all intersecting fibers
    cyl_bnd_nodes = [];
    
    for num = 1 : num_fibers
        start_node_num = fibers(2*num-1); % First node of fiber
        end_node_num = fibers(2*num); % Second node of fiber
        node1 = [nodes(3*start_node_num-2) nodes(3*start_node_num-1) nodes(3*start_node_num)];
        node2 = [nodes(3*end_node_num-2) nodes(3*end_node_num-1) nodes(3*end_node_num)];

        dirv=node2-node1; % directional vector of the line

        % Calculate s1 and s2- the possible intersection points of a ray
        % with an infinite cylinder.
        % Values between 0 and 1 mean the intersection is actually on the line.
        % Values outside of 0 and 1 means the the intersection is along the projected extension of the line.
        % Values that are imaginary means the the intersection is not real.        
        s1=(-(2*(node1(rad1)-prot_node1(rad1))*dirv(rad1)+2*(node1(rad2)-prot_node1(rad2))*dirv(rad2))+sqrt((2*(node1(rad1)-prot_node1(rad1))*dirv(rad1)+2*(node1(rad2)-prot_node1(rad2))*dirv(rad2))^2-4*(dirv(rad1)^2+dirv(rad2)^2)*((node1(rad1)-prot_node1(rad1))^2+(node1(rad2)-prot_node1(rad2))^2-R^2)))/(2*(dirv(rad1)^2+dirv(rad2)^2));
        s2=(-(2*(node1(rad1)-prot_node1(rad1))*dirv(rad1)+2*(node1(rad2)-prot_node1(rad2))*dirv(rad2))-sqrt((2*(node1(rad1)-prot_node1(rad1))*dirv(rad1)+2*(node1(rad2)-prot_node1(rad2))*dirv(rad2))^2-4*(dirv(rad1)^2+dirv(rad2)^2)*((node1(rad1)-prot_node1(rad1))^2+(node1(rad2)-prot_node1(rad2))^2-R^2)))/(2*(dirv(rad1)^2+dirv(rad2)^2));
        
        nsargs = isreal(s1) + isreal(s2); % Number of real s results
        
        % Calculate c1 and c2- the possible intersection points of a ray 
        % and the end planes of the cylinder. 
        % Values between 0 and 1 mean the intersection is actually on the line.
        % Values outside of 0 and 1 means the the intersection is along the projected extension of the line.
        % Values that are infinite means the the intersection is a line (parallel to plane).
        n1 = prot_node1(direction); 
        n2 = prot_node2(direction);
        
        c1 = (n1-node1(direction))/dirv(direction);
        c2 = (n2-node1(direction))/dirv(direction);
        
        ncargs = isreal(c1) + isreal(c2); % Number of real c results
        
        % Calculate radial positions of fiber nodes.
        switch direction
            case 1
                d_node1 = norm(node1(2:3));
                d_node2 = norm(node2(2:3));
            case 2
                d_node1 = norm([node1(1),node1(3)]);
                d_node2 = norm([node2(1),node2(3)]);
            case 3
                d_node1 = norm(node1(1:2));
                d_node2 = norm(node2(1:2));
        end
        
        ntotargs = nsargs + ncargs; % Total number of real intersection values
        
        if ntotargs == 4 % If [s1 s2] and [c1 c2] are all real
            s = sort([s1, s2]);
            c = sort([c1, c2]);
            p1 = max(s(1),c(1));
            p2 = min(s(2), c(2)); 
            
            if p1 <= p2 % If the ranges [t1 t2] and [c1 c2] overlap
                if p1 ~= p2 % If the ranges intersect over an interval
                    if (p1 >= 0) && (p1 <= 1) && (p2 >=0) && (p2 <=1) % If the intersection interval is on fiber length
                        % Create nodes at the intersections in the intersection interval
                        int1 = node1+p1*dirv; % Coordinates of first intersection
                        int2 = node1+p2*dirv; % Coordinates of second intersection
                        
                        nodes = [nodes int1]; % Create first new node
                        fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1
                        
                        % Now make new fiber from int1 to int2 (this will
                        % be the one discretized)
                        fibers(end+1) = length(nodes)/3;
                        nodes = [nodes int2]; % Create second new node
                        fibers(end+1) = length(nodes)/3;

                        % Intersecting fiber is the fiber created above
                        int_fiber_num = length(fibers)/2; 

                        % Next, create third fiber from int2 - node 2
                        fibers(end+1) = length(nodes)/3;
                        fibers(end+1) = end_node_num;
                        
                        % Finally, discretize fiber in between two intersections; 
                        [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                        cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                        
                        % Track all intersecting fibers and their associated fiber segments
                        intfib_index = concat_diff(intfib_index, split_fiber_index);  
                        
                    elseif (p1 >= 0) && (p1 <= 1) % Just p1 is on the fiber length
                        if (d_node1 < R) && (node1(direction) <= prot_node2(direction)) && (node1(direction) >= prot_node1(direction)) % If node1 is within cylinder
                            % Make node at intersection and discretize
                            % between point and node1
                            int1 = node1+p1*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Intersecting fiber is the fiber created above
                            int_fiber_num = num; 

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                            
                        else
                            % Make node at intersection and discretize
                            % between point and node2
                            int1 = node1+p1*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Intersecting fiber is the fiber created above
                            int_fiber_num = length(fibers)/2; 

                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                                     
                        end
                    elseif (p2 >=0) && (p2 <=1) % Just p2 is on the fiber length
                        if (d_node1 < R) && (node1(direction) <= prot_node2(direction)) && (node1(direction) >= prot_node1(direction)) % If node1 is within cylinder
                            % Make node at intersection and discretize
                            % between point and node1
                            int1 = node1+p2*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Intersecting fiber is the fiber created above
                            int_fiber_num = num; 

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                            
                        else
                            % Make node at intersection and discretize
                            % between point and node2
                            int1 = node1+p2*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Intersecting fiber is the fiber created above
                            int_fiber_num = length(fibers)/2; 

                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                                     
                        end
                    else % Intersections are outside fiber length
                        % If both nodes are within the cylinder,
                        % discretize. Otherwise don't care.
                        if (d_node1 < R) && (d_node2 < R) && (node1(direction) <= prot_node2(direction)) && (node1(direction) >= prot_node1(direction)) ...
                                && (node2(direction) <= prot_node2(direction)) && (node2(direction) >= prot_node1(direction)) 
                        
                            % Discretize fiber as is
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, num, sect_no); 
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index);          
                        end
                       
                    end
                else  % Ranges intersect at one point (p1 == p2)
                    if (p1 >= 0) && (p2 <= 1) % If intersection is on fiber length
                        if (d_node1 < R) && (node1(direction) <= prot_node2(direction)) && (node1(direction) >= prot_node1(direction)) % If node1 is within cylinder
                            % Make node at intersection and discretize
                            % between point and node1
                            int1 = node1+p1*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Intersecting fiber is the fiber created above
                            int_fiber_num = num; 

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                            
                        else
                            % Make node at intersection and discretize
                            % between point and node2
                            int1 = node1+p1*dirv; % Coordinates of intersection
                            nodes = [nodes int1]; % Create first new node
                            fibers(2*num) = length(nodes)/3; % Edit fiber to now go from node1 - int1

                            % Make new fiber - set intersection as first node, 
                            % old end node of original fiber is second node
                            fibers(end+1) = length(nodes)/3;
                            fibers(end+1) = end_node_num;
                            
                            % Intersecting fiber is the fiber created above
                            int_fiber_num = length(fibers)/2; 

                            % Discretize the fiber going from node1-int1
                            [nodes, fibers, split_fiber_index, c_nodes] = split_fibers(nodes, fibers, int_fiber_num, sect_no);
                            cyl_bnd_nodes = [cyl_bnd_nodes, c_nodes];
                            
                            % Track all intersecting fibers and their associated fiber segments
                            intfib_index = concat_diff(intfib_index, split_fiber_index); 
                                     
                        end
                    end
                end
            end
        else
            %fprintf('Weird case encountered! \n');
        end
    end