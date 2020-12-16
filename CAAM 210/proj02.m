function qdrive
% qdrive
% qdrive defines the polynomials to be checked, the intervals of x and y,
% and the number of iterations. qdrive then runs qnewt for each function
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 02
% Last Modified: September 28, 2019

% Initals
q1 = [1 0 -0.84 0 -0.16];
q2 = [1 0 -0.84 -0.1 -0.16];
q3 = [1 -0.1 -0.84 0 -0.16];
q4 = [1 -0.1i -0.84 0 -0.16];
xt = [.45 .0001 .55];
yt = [-.05 .0001 .05];
maxiter = 20;
% Iterates over the initial functions
testqs = [q1; q2; q3; q4];
[rows,~] = size(testqs);
for idx1 = 1:rows
    figure(idx1)
    qnewt(testqs(idx1,:), xt, yt, maxiter)
end
end

function [dq] = mypolyder(q)
% [dq] = mypolyder(q)
% mypolyder finds the derivative of the polynomial q
% Inputs: q, the vector of the coefficients of the polynomial
% Outputs: dq, the vector of the coefficients of the derivative

len = length(q);
for i= 1:(len-1)
    dq(i) = q(i)*(len-i);
end
end

function qnewt(q,xt,yt,maxiter)
% qnewt(q,xt,yt,maxiter)
% qnewt takes the initial parameters and plots a graph of the coordinates
% that converge to a root, one color for each root
% Inputs: q, the vector of the coefficients of the polynomial
%         xt, the x interval of the points to check
%         yt, the y interval of the points to check
%         maxiter, the number of times to run newtons method
% Outputs: a plot of which coordinates converge

% Create field of the complex points to check
xrng = xt(1):xt(2):xt(3);
yrng = yt(1):yt(2):yt(3);
[X,Y] = meshgrid(xrng,yrng);
Z = X+1i*Y;
% Runs newton's method on the field
dq = mypolyder(q);
for k=1:maxiter
    f = polyval(q,Z);
    df = polyval(dq,Z);
    Z = Z - f./df;
end
% Check which points converge to which root, plot
qroots = roots(q);
hold on
colors = ['y','r','g','b']; % I like these colors more
for idx2 = 1:length(qroots)
    [i1,j1] = find(abs(Z-qroots(idx2))<0.1);
    plot(j1,i1,strcat('.', colors(idx2)),'MarkerSize',1)
    title(qlab(q))
    axis tight
end
hold off
end

function [qlab] = qlab(q)
% [qlab] = qlab(q)
% qlab generates a title for the plots, assumes that the first element of 
% 1 is nonzero
% Inputs: q, the polynomial
% Outputs: qlab, the title

% Initalizes qlab
len = length(q);
if q(1) == 1
    qlab = strcat('z^',num2str(len-1));
else
    qlab = strcat(num2str(q(1)),'z^',num2str(len-1));
end
% Creates the rest of qlab
for k = 2:(len)
    if q(k) ~= 0
        qlab = strcat(qlab,'+',num2str(q(k)),'z^',num2str(len-k));
    end
end
% Cleans up title
qlab =strrep(qlab, '+-', '-');
qlab =strrep(qlab, '+0', '');
qlab =strrep(qlab, '+1z', '+z');
qlab =strrep(qlab, 'z^1', 'z');
qlab =strrep(qlab, 'z^0', '');
end