function proj11
% proj11
% driver function, runs evo and plots the fraction of cooperators for 2
% sets of initial conditions
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 11
% Last Modified: December 3, 2019


% first set of initial conditions
close all
M = 39;
N = 39;
rnd = 200;
b = 1.9;

% runs evo, plots the final board
[~,C1] = evo(M,N,b,rnd);

% plots the proportion of cooperators
figure
plot(0:rnd, C1)
title("M=39, N=39")
xlabel("rounds")
ylabel("fraction of cooperators")



% second set of initial conditions
M = 35;
N = 35;
rnd = 50;

% runs evo, plots the final board
[~,C2] = evo(M,N,b,rnd);

% plots the proportion of cooperators
figure
plot(0:rnd, C2)
title("M=35, N=35")
xlabel("rounds")
ylabel("fraction of cooperators")


% What is different about the outcome of this game? 

% The game played on a 35x35 size board results in all of the squares
% becoming defectors by round 50. It is different because there exists an
% stable equilibrium state for the game.

% What would happen if we played more rounds?

% There would be no change to the state of each of the players. Because all
% the players are defectors, all players recieve a score of 0, and remain
% defectors.
end


function [A,C] = evo(M,N,b,rnd)
% [A,C] = evo(M,N,b,rnd)
% evolution function, runs the evolution based on the initial conditions
% Inputs: 
%     - M, the number of rows in the board
%     - N, the number of columns in the board
%     - b, the score parameter for a defector
%     - rnd, the number of rounds to run
% Outputs:
%     - C, the vector of the proportion of cooperators


% initialize A, An, matricies
A = ones(M,N);
An = A;
An(ceil(M/2),ceil(N/2)) = 0;
C = [1];

% go through evolutionary rounds
for i = 1:rnd
    % evodisp(A, An)
    A = An;
    [S]  = score(A,b);
    [An] = advance(S,A);
    C(i+1) = sum(sum(A))/(M*N);
end

% plot the final condition of the board
figure
evodisp(A, An)
title("round " + string(rnd))
end

function [N] = neighbors(A,i,j)
% [N] = neighbors(A,i,j)
% finds the neighbors of a given element
% Inputs:
%     - A, the matrix to be operated on
%     - i, the row location of the element
%     - j, the column location of the element
% Outputs:
%     - N, the matrix of neighbors


[M, N] = size(A);

wrap = @(x, S) (1 + mod(x-1, S));
vec_nbr = @(y, S) ([wrap(y-1, S), y, wrap(y+1,S)]);

N = A(vec_nbr(i, M), vec_nbr(j, N));
end


function [S]  = score(A,b)
% [S]  = score(A,b)
% scores a given board based on the roles of each player
% Inputs:
%     - A, the matrix to be operated on
%     - b, the defector score parameter
% Outputs:
%     - S, the matrix of scores


% get matrix size
[M, N] = size(A);

% score each tile
for i = 1:M
    for j = 1:N
        nbrs = neighbors(A, i, j);
        if A(i,j) == 1
            S(i,j) = sum(sum(nbrs));
        else
            S(i,j) = b * sum(sum(nbrs));
        end
    end
end
end



function [An] = advance(S,A)
% [An] = advance(S,A)
% rebuilds and adjusts player strategies based on score matric
% Inputs:
%     - A, the matrix to be operated on
%     - S, the matrix of scores
% Outputs:
%     - An, the next player board matrix


% initialize An, get size A
An = A;
[M,N] = size(A);

% obtain the strategy with the highest score
for i = 1:M
        for j = 1:N
            nbrS = neighbors(S, i, j);
            [imax, jmax] = find(max(max(nbrS)) == nbrS);
            nbrA = neighbors(A, i, j);
            An(i,j) = nbrA(imax(1), jmax(1));
        end
end

% alternate method with augmented matrix

% S_aug = augment(S);
% A_aug = augment(A);
% 
% for i = 1:M
%         for j = 1:N
%             nbrS = neighbors(S_aug, i+1, j+1);
%             [imax, jmax] = find(max(max(nbrS)) == nbrS);
%             nbrA = neighbors(A_aug, i+1, j+1);
%             An(i,j) = nbrA(imax(1), jmax(1));
%         end
% end
end

function evodisp(A,An)
% evodisp(A,An)
% plots the evolution of the player board
% Inputs:
%     - A, the matrix to be operated on
%     - An, the next player board matrix
% Outputs: none


% initials
[M,N] = size(A);
D = zeros(M,N,3);

% evolutionary differences
Asum = A+An;
Adif = An-A;

% based on difference and sum matricies, color plot
D(:,:,3) = (Asum==2);
D(:,:,1) = (Asum==0) + (Adif==-1);
D(:,:,2) = (Adif==1) + (Adif==-1);

% % worse method
% for i = 1:M
%     for j = 1:N
%         D(i,j,:) = color_code(Asum(i,j), Adif(i,j));
%     end
% end

% plot array
image(D)
axis square
end

function cc = color_code(B, C)
% cc = color_code(B, C)
% unused func that determines the color code to be filled in
% Inputs:
%     - B, the sum matrix
%     - C, the difference matrix
% Outputs:
%     - cc, the color code vector


if B == 2
    cc = [0 0 1]; % if a tile remains cooperative
elseif B == 0
    cc = [1 0 0]; % if a tile remains defective
else
    if C == 1
        cc = [0 1 0]; % switch to coop
    else
        cc = [1 1 0]; % defect
    end
end
end

function [F] = augment(A)
% [F] = augment(A)
% unused function to create augmented matric
% Inputs:
%     - A, the start matrix
% Outputs:
%     - F, the wrapped matrix


[M,N] = size(A);
F = zeros(size(A) + [2 2]);
F(2:(M+1), 2:(N+1)) = A;
F([1 M+2], 2:(N+1)) = A([1 M], :);
F(2:(N+1), [1 N+2]) = A(:, [1 N]);
F([1 M+2], [1 N+2]) = A([1 M], [1 N]);
end