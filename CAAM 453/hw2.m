function hw2

warning off

close all

% problem 1

% initialize p values, points a and b
p = [1 2 8 Inf]';
a = [0 0];
b = [1 2];

% initialize titles
tt1 = "||x||_"+p+" = 1";
tt1(4) = '||x||_\infty = 1';
tt2 = "||x-a||_{" + p + "} = ||x-b||_{" + p + "}";
tt2(4) = "||x-a||_{\infty} = ||x-b||_{\infty}";

for i = 1:4
    figure(1)
    subplot(2,2,i)
    fimplicit(@(x,y) norm([x y],p(i))-1); % unit circle
    axis([-2 2 -2 2]) 
    axis square
    title(tt1(i))
    
    figure(2)
    subplot(2,2,i)
    hold on
    % plots the implicit functions representing the bisection
    fimplicit(@(x,y) norm([x y]-a,p(i))-norm([x y]-b,p(i))); 
    plot([0 1],[0 2],'ko') % plots a and b
    hold off
    axis([-.5 2.5 -.5 2.5]) 
    axis square 
    title(tt2(i))
end

% Problem 3 (50 points)

% A = rand(3);
% [V,R] = householder_qr(A);
% b = rand(3,1);
% x = solve_linear_system_with_qr(V,R,b);
% A\b, x
% inv(A)
% invert_using_qr(A)

% Give the number of operations that is requires to compute A−1 
% using that algorithm. 

% You are to test your codes for correctness using the following script:
% The number Diff should be on the order of 10−10 or less.
n = 50;
A = randn(n);
Ainv = invert_using_qr(A); 
Diff = norm(Ainv*A - eye(n));
disp("norm of (Ainv*A - I): " + Diff)

disp("running time of Ainv is O(n^3), for the approximately n or so matrix-matrix operations in QR decomposition, which dominates solve_linear_system")
end


function [V,R] = householder_qr(A)
% Input:
%    A    is a square n x n real matrix
% Output:
% V is a "cell array": vk is (a vector of length n-k+1) is stored in V{k}. 
% R is an upper triangular n x n matrix.

[~, n] = size(A);
V = eye(n);
for j = 1:n 
    aj = A(j:n,j);
    vj = aj + sign(aj(1))*norm(aj, 2)*[1; zeros(n-j,1)];
    Hj = eye(n-j+1) - 2*(vj*vj')/(vj'*vj); 
    Qj = [eye(j-1), zeros(j-1,n-j+1); 
          zeros(n-j+1,j-1), Hj];
    A = Qj*A; 
    V = V*Qj; 
end
R = A;

end

function x = solve_linear_system_with_qr(V,R,b)
% Input:
% V,R is the QR factorization of A that is given by householder_qr(A). 
% b is the right hand side of the equation.
% Output:
% x is the solution of Ax = b.

if det(V*R) == 0
    disp("A is singular") % test for singularity
    return
end

[~, n] = size(R);
b = V'*b;
for j = n:-1:1
    x(j,:) = b(j,:)./R(j,j);
    b = b - x(j,:).*R(:,j);
end
end

function Ainv = invert_using_qr(A)
% Input:
%     A   is a n x n real matrix
% Output:
%     Ainv is the inverse of A

[V,R] = householder_qr(A);
[~, n] = size(A);
Ainv = solve_linear_system_with_qr(V,R,eye(n));
end
