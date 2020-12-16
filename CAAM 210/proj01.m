% Driver

function cooldrive
% cooldrive
% cooldrive is a driver function that first determines and plots how the cooling
% of a bar is affected by the bar length, then shows the relationships
% between the tolerance and number of iterations for both Newton's method
% and for the bisection function
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 01
% Last Modified: September 23, 2019

% Generate bisection roots over a range of bar lengths
a = 0.1;
b = 3;
L = 1:0.1:4;
tol = 0.01;

for lengthindex = 1:length(L)
    x(lengthindex) = bisect(a,b,L(lengthindex),tol)
    b = 1.1*x(lengthindex);
end

% First plot, bar length versus decay rate
figure(1)
plot(L, x.^2, 'x-')
title('Cooling Rate vs. Bar Length')
ylabel('Length, L')
xlabel('Decay Rate, x^2')

% Generate iteration numbers for Newton's and the bisection method over a
% series of tolerances
L=1;
a = 0.1;
b = 3;
ex = (a+b)/2;

for j = 1:8
    tol(j)=10.^(-j);
end

for tolidx1 = 1:length(tol)
    [x(tolidx1), iterx(tolidx1)] = bisect(a,b,L,tol(tolidx1));
end

for tolidx2 = 1:length(tol)
    [y(tolidx2), itery(tolidx2)] = newton(ex,L,tol(tolidx2));
end

% Second plot, tolerance versus iterations for both methods
figure(2)
semilogx(tol, iterx, 'x-', tol, itery, 'o-')
title('Bisection vs. Newton')
ylabel('# of Iterations')
xlabel('Tolerance')
legend('Bisection','Newton')


end

% Bisection Method

function [x,iter] = bisect(a,b,L,tol)
% [x,iter] = bisect(a,b,L,tol)
% The function returns a root of the bar cooling function using the
% bisection method
% Inputs: a, lower bound of the testing interval
%         b, upper bound of the testing interval
%         L, length of bar
%         tol, the tolerance of the root and the result
% Output: x, the root of the bar cooling function within [a,b]
%         iter, the number of iterations to get the a root within the given
%         tolerance

% x = [];
if bcooln(a,L)*bcooln(b,L) > 0
    disp("there may not be a root in the interval, please try again")
    return
else
    iter = 0;
    x = (a+b)/2;
    while abs(bcooln(x,L)) > tol
        if bcooln(a,L)*bcooln(x,L) < 0
            b=x;
        else
            a=x;
        end
        x = (a+b)/2;
        iter = iter+1;
    end
end
end

% Newton's Method 

function [x, iter] = newton(x,L,tol)
% [x, iter] = newton(x,L,tol)
% The function returns a root of the bar cooling function using Newton's
% method
% Inputs: x, initial guess
%         L, length of bar
%         tol, the tolerance of the root and the result
% Output: x, the root of the bar cooling function near the initial x
%         iter, the number of iterations to get the a root within the given
%         tolerance

iter = 0;
while abs(bcooln(x,L)) > tol
    x = x -bcooln(x,L)/bcoolndx(x,L);
    iter = iter +1;
end
end

% Bar cooling function and its derivative

function [y] = bcooln(x,L)
% [y] = bcooln(x,L)
% The least positive root is the square root of the rate of heat decay
% Inputs: x, x value on the x axis
%         L, length of bar
% Output: the resultant value of sin(Lx)+xcos(Lx)

y = sin(x.*L) + x.*cos(x.*L);
end

function [y] = bcoolndx(x,L)
% [y] = bcooln(x,L)
% The derivative of the above function
% Inputs: x, x value on the x axis
%         L, length of bar
% Output: the derivative the original bar cooling method, evaluated for x
% and L.

y = (L+1).*cos(x.*L) - L.*x.*sin(x.*L);
end


