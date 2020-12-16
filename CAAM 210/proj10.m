function proj10
% proj10
% finds and plots the maximum clique, minimum dominating set, and the
% maximum k-plex
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 10
% Last Modified: November 23, 2019


close all

% read file, initialize relavant data
input = textread("social_network_graph.txt");
num_nodes = input(1,1);
num_edges = input(1,2);
C = input(2:end,:); % connections matrix

% initialize adjacency matrix, plot the graph
adj = adjacency(C, num_nodes, num_edges);
G = graph(adj);
plot(G);
title("Original Social Network Graph")

% find and plot the maximum clique
plex(adj,num_nodes,1);
title("Maximum Clique")

% find and plot the maximum k-plex
plex(adj,num_nodes,5);
title("Maximum 5-plex")

% find and plot the minimum dominating set
domset(adj,num_nodes,G);
title("Minimum Dominating Set")
end


function [adj] = adjacency(C,nodes,edges)
% [adj] = adjacency(C,nodes,edges)
% finds the adjacency matrix
% Inputs:
%     - C, the matrix of connections
%     - nodes, the number of nodes
%     - edges, the number of edges
% Outputs:
%     - adj, the adjacency matrix


% initialize adj
adj = zeros(nodes,nodes);

% symmetrically fill adj
for row = 1:edges
    adj(C(row,1),C(row,2)) = 1;
    adj(C(row,2),C(row,1)) = 1;
end
end



function plex(adj,nodes,k)
% plex(adj,nodes,k)
% finds and plots the maximum k-plex of a graph
% Inputs:
%     - adj, the adjacency matrix
%     - nodes, the number of nodes
%     - k, the parameter for the k-plex
% Outputs: none


% initialize parameters for linear integer programming
f = -ones(nodes,1); % f is negative for maximization
int = 1:nodes;
A = abs(adj - ones(size(adj))) - eye(size(adj));

% let d be the number of nodes that a node is not connected to
d = zeros(nodes,1);
for i = 1:nodes
    d(i) = sum(A(i,:));
end

A2 = A + diag(d-k+1); % finalized A matrix:
b = d;
ub = ones(nodes,1);
lb = zeros(nodes,1);
Aeq = []; beq = [];

% find solution to the linear integer problem
[opt_sol, fval] = intlinprog(f,int,A2,b,Aeq,beq,lb,ub);

% modify the adjacency matrix to plot only the k-plex
[I,~] = find(opt_sol > 1e-3);
A_subset = adj(I,I);
nodenames = string(I);

% plot it
figure
G = graph(A_subset, nodenames);
plot(G)

% display size of k-plex
disp("size of maximum " + string(k) + "-plex:")
disp(-fval)
disp("")
end



function domset(adj,nodes,G)
% domset(adj,nodes,G)
% finds and highlights the dominating set of a graph G
% Inputs:
%     - adj, the adjacency matrix
%     - nodes, the number of nodes
%     - G, the given graph
% Outputs: none


% initialize relevant parameters for the linear integer programming problem
A = adj + eye(size(adj));
f = ones(nodes,1);
b = ones(nodes,1);
int = 1:nodes;
ub = ones(nodes,1);
lb = zeros(nodes,1);
Aeq = []; beq = [];

% solve the problem
[opt_sol, fval] = intlinprog(f,int,-A,-b,Aeq,beq,lb,ub);

% highlight the points in the dominating set
figure
[I,~] = find(opt_sol > 1e-4);
Gplot = plot(G);
highlight(Gplot, I, 'NodeColor','g')

% display the size of the set
disp("size of minimum dominating set:")
disp(fval)
disp("")
end




