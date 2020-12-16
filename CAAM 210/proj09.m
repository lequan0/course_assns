function proj09
% proj09
% plots the sparsity patterns and deformations of n n-section bridges.
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 09
% Last Modified: November 16, 2019

close all


AA = {};

for nos = 1:3
    [adj, xc, yc, len] = build_basic_bridge(nos); % build bridge
    AA(nos) =  {adj};
    [force1, force2] = f(nos); % get force vectors
    % deformations
    [dx1,dy1,~,~,~] = deform_basic_bridge(nos,adj,xc,yc,len,force1);
    [dx2,dy2,~,~,~] = deform_basic_bridge(nos,adj,xc,yc,len,force2);
    % plot
    plotter(adj,xc,yc,dx1,dy1,dx2,dy2,nos)
end

% questions

% 1. If you are driving a light-weight car (0.01 units of weight), 
% what is the maximum length (in nos) 
% the bridge could be while still being considered safe? 
% (What is your criterion for ?safe??)

% Given the ultimate tensile stress, the bridge is safe insofar as the
% force experienced by any given fiber does not exceed the aforementioned
% stress.

% Without this information though, I'd say that 6 sections is as far as I
% would go.

%========================================================================

% 2. What about if you are driving an 18-wheeler truck 
% (0.05 units of weight)?

% The same arguement for question 1 applies to question 2, and I would 
% drive on a bridge with at most 3 sections.

%========================================================================

% 3. How does the spy of the adjacency matrix A 
% for a particular bridge relate to how it deforms?

% A nonzero point in the sparsity matrix represents a point that relates
% the displacement of a point and its elongation. The first fiber, for
% example, is elongated by the same distance the first node is displaces in
% the x-direction.

%========================================================================

% 4. What trends do you notice as nos increases? 
% (Relate your answer to the spy of a particular bridge, 
% how the bridge deforms under different stress conditions, etc.)

% As the number of sections increases, the sparsity pattern of the bridge
% repeates itself along it's own diagonal, demonstrating that the sparsity
% of the adjacency matrix follows the increase in section length.

% The minimum of the deformed bridge, given a constant force, decreases
% with an increase in the number of sections.
end

function [y,z] = f(x)
% [y,z] = f(x)
% gets the force vectors, where the forces are applied on each of the
% nodes.
% Inputs: x, the number of nodes
% Outputs: y, the force vector for a small car
%          z, the force vector for a large car


a = [0 -.01 0 0]'; % a negative force on the bottom node
y = a;

% iterate to make y
for i = 1:x
    y = [y; a];
end
z = 5.*y;
end


function plotter(adj,xc,yc,dx1,dy1,dx2,dy2,nos)
% plotter(adj,xc,yc,dx1,dy1,dx2,dy2,nos)
% plots the sparsity pattern and relative deformation of a bridge
% Inputs: adj, the adjacency matrix
%         xc, the undeformed x-coordinates of each fiber
%         yc, the undeformed y-coordinates of each fiber
%         dx1, the first deformed x-coordinates of each fiber
%         dy1, the first deformed y-coordinates of each fiber
%         dx2, the second deformed x-coordinates of each fiber
%         d21, the second deformed y-coordinates of each fiber
%         nos, the number of sections
% Outputs: none


figure
% plot the sparsity pattern
subplot(1,2,1)
spy(adj)
title('sparsity pattern of adj, and nos = ' + string(nos))

% plot the three forms of each matrix
subplot(1,2,2)
p1=line(xc', yc', 'Color','blue'); % undeformed
p2=line(dx1',dy1', 'Color', 'red'); % deform 1
p3=line(dx2',dy2', 'Color', 'green'); % deform 2
title('deformation of truss bridge')
legend([p1(1) p2(1) p3(1)], 'undeformed','deform 1','deform 2')
axis equal 
set(gca,'XTick',[], 'YTick', []) % gets rid of xy labels
set(gca,'XColor', 'none','YColor','none') % gets rid of xy axis lines
end

function [adj, xc, yc, len] = build_basic_bridge(nos)
% [adj, xc, yc, len] = build_basic_bridge(nos)
% builds a bridge, given the number of sections
% Inputs: nos, the number of sections
% Outputs: adj, the adjacency matrix
%          xc, the undeformed x-coordinates of each fiber
%          yc, the undeformed y-coordinates of each fiber
%          len, the length vector of each fiber


% calculate some helpful numbers
num_nodes = nos*2+2;
num_edges = nos*5+5;
s = 1/sqrt(2);
% initialize the return values
adj = sparse(num_edges, 2*num_nodes);
xc = sparse(num_edges, 2);
yc = sparse(num_edges, 2);
len = ones(num_edges, 1);
% build the left side of bridge (edges 1 and 2)
adj([1 2], [1 3 4]) = [1 0 0;
                       0 s s];

xc([1 2],:) = [0 1;
               0 1];

yc([1 2],:) = [0 0;
               0 1];

len(2) = 1/s;

% build the middle of bridge
for i = 1:nos
    row_frame = (5*i-2):(5*i+2);
    col_frame = (4*i-3):(4*i+4);
    adj(row_frame,col_frame)= [0 -1  0  1  0  0  0  0;
                               0  0 -s  s  s -s  0  0;
                               0  0 -1  0  0  0  1  0;
                              -s -s  0  0  0  0  s  s;
                              -1  0  0  0  1  0  0  0];

    xc(row_frame, :) = [i i;
                        i i+1;
                        i i+1;
                        i i+1;
                        i i+1];
    
    yc(row_frame, :) = [0 1;
                        1 0;
                        1 1;
                        0 1;
                        0 0];
    
    len(row_frame) = [1 1/s 1 1/s 1];
end

% build the right side of bridge
adj((end-2):(end) ,(end-3):(end)) = [0 -1  0  1;
                                     0  0 -s  s;
                                    -1  0  0  0];
xc((end-2):(end), :) = [nos+1 nos+1;
                        nos+1 nos+2;
                        nos+1 nos+2];
yc((end-2):(end), :) = [0 1;
                        1 0;
                        0 0];
len(end-1) = 1/s;
end


function [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force)
% plotter(adj,xc,yc,dx1,dy1,dx2,dy2,nos)
% plots the sparsity pattern and relative deformation of a bridge
% Inputs: nos, the number of sections
%         adj, the adjacency matrix
%         xc, the undeformed x-coordinates of each fiber
%         yc, the undeformed y-coordinates of each fiber
%         len, the length vector of each fiber
%         force, the force applied to each directon on each node

% Outputs: dx, the deformed x-coordinates of each fiber
%          dy, the deformed y-coordinates of each fiber
%          work, the work done by the bridge to resist deformation
%          X, the x displacements for each node
%          Y, the y displacements for each node


% calculate the displacements
stiffness = adj'*diag(1./len)*adj;
displacements = stiffness\force;

% find work done
nos, force
work = displacements'*force

% seperate x and y displacements
X = displacements(1:2:end);
Y = displacements(2:2:end);

% initialize the return values
dx = zeros(size(xc));
dy = zeros(size(yc));

% deform the left side of bridge
dx([1 2],:) = xc([1 2],:) + [0    X(1);
                             0    X(2)];
dy([1 2],:) = yc([1 2],:) + [0    Y(1);
                             0    Y(2)];

% deform the middle of bridge
for i = 1:nos
    row_frame = (5*i-2):(5*i+2);
    dx(row_frame,:) = xc(row_frame,:) + [X(2*i-1) X(2*i);
                                         X(2*i)   X(2*i+1);
                                         X(2*i)   X(2*i+2);
                                         X(2*i-1) X(2*i+2);
                                         X(2*i-1) X(2*i+1)];
    
    dy(row_frame,:) = yc(row_frame,:) + [Y(2*i-1) Y(2*i);
                                         Y(2*i)   Y(2*i+1);
                                         Y(2*i)   Y(2*i+2);
                                         Y(2*i-1) Y(2*i+2);
                                         Y(2*i-1) Y(2*i+1)];
end

% deform the right side of bridge
dx((end-2):end, :) = xc((end-2):end,:) + [X(end-1) X(end);
                                          X(end)   0;
                                          X(end-1) 0];

dy((end-2):end, :) = yc((end-2):end,:) + [Y(end-1) Y(end);
                                          Y(end)   0;
                                          Y(end-1) 0];
end