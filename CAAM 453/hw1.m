%% driver

function hw1
% function driver
% inputs: none
% outputs: none
% quan le, CAAM 453, fall 2020, HW 1
% last modified: sept 23, 2020
close all

%% problem 1

% initialize constants
Ti = 20; Ts = -15;
a = 0.183e-6;
t = 60*24*60*60;

% initialize funtions
syms x
f = Ts + (Ti-Ts)*erf(x/(2*sqrt(a*t)));
fx = diff(f);

% plot f
figure 
fplot(f,[0 5])
title("Plot of f(x) on the interval [0, 5]")

% determine roots
f = matlabFunction(f);
fx = matlabFunction(fx);
t = 1e-9;
r1 = bisect(f,0,5,t);
r2 = newton(f,fx,.01,t);
% r2 = newton(f,fx,5,t);
disp("root from bisection:       " + r1)
disp("root from newton's method: " + r2)

drawnow

%% problem 2

% initialize function, find roots
p = [1 0 -2 2];
r = roots(p);
disp("roots of z^3-2z+2:")
disp(r)

figure
qnewt(p,1001,5,20)
title("Newton Fractal of z^3-2z+2")

% checking for nonconvergent guesses
% to view the behaviour, unsuppress the output in newton()
g = @(x) x^3 - 2*x +2;
dg = @(x) 3*x^2 - 2;
% newton(g,dg,-0.5+.85i,t)
% newton(g,dg,-0.8-1.75i,t)
% newton(g,dg,1,t)
% newton(g,dg,1+.1i,t)


%% problem 3

% finding x*

% inefficient method using arctan x
f = @(x) atan(x);
df = @(x) 1./(x.^2+1);

x = linspace(0,2,20000);
Z = x;
for k=1:20
    Z = Z - f(Z)./df(Z);
end
x(find(abs(Z)<0.001, 1, 'last'));


% more effective method using g
syms x
g = 2*x - (x.^2+1).*atan(x);
dg = diff(g);
g = matlabFunction(g);
dg = matlabFunction(dg);

r1 = bisect(g,-2,0,t);
r2 = newton(g,dg,2,t);

% numerical tests for behavior of x*
for i = 10.^-(1:10)
    x = [-r2-i, -r2+i, -i, i, r2-i, r2+i];
    disp("epsilon = "+ i)
    for xi = x
        y = newton(f,df,xi,t);
        if isnan(y)
            disp("newton method for x=" + xi + ": diverges")
        else
            disp("newton method for x=" + xi + ": " + y)
        end
    end
    disp(" ")
end


end

%% supporting functions

function [x] = bisect(f,a,b,tol)
% The function returns a root of f in [a, b], if it exists
% Inputs: f, function to be tested
%         a, b, bounds of the testing interval
%         tol, the tolerance of the root and the result
% Output: x, root

if f(a)*f(b) > 0 % check for root in the interval
    disp("there may not be a root in the interval, please try again")
    return
else
    x = (a+b)/2;
    while abs(f(x)) > tol
        if f(a)*f(x) < 0
            b=x;
        else
            a=x;
        end
        x = (a+b)/2;
    end
end
end

function [x] = newton(f,df,x,tol)
% implementation of newton's method to find a root of f
% inputs: f, function of interest
%         x, initial guess
%         df, derivative of f
%         tol, the tolerance
% output: x, the root 
while abs(f(x)) > tol 
    x = x -f(x)/df(x);
end
end

function qnewt(q,n,l,maxiter)
% newton's method for complex valued polynomials
% Inputs: q, the vector of the coefficients of the polynomial
%         n, number of points for each axis
%         l, length of the testing area
%         maxiter, the number of times to run newtons method
% outputs: plot of newton basins

% Create field of the complex points to check
x = linspace(-l,l,n);
y = linspace(-l,l,n);
[X,Y] = meshgrid(x,y);
Z = X-1i*Y;
% Runs newton's method on the field
dq = polyder(q);
for k=1:maxiter
    f = polyval(q,Z);
    df = polyval(dq,Z);
    Z = Z - f./df;
end
% Check which points converge to which root, plot
r = roots(q);
M = zeros(size(Z));
for i = 1:length(r)
    M = M + i*(abs(Z-r(i)) < 0.01);
end
imagesc(M)
colorbar

end