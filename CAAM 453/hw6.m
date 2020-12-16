function hw6
% quan le
% caam 453: Numerical Analysis II
% problem set 6

close all

problem1
problem2
problem3
end

function problem1
eps2 = abs(2*(3/2 - 1) - 1)
eps3 = abs(3*(4/3 - 1) - 1)
eps4 = abs(4*(5/4 - 1) - 1)
eps5 = abs(5*(6/5 - 1) - 1)
eps % This is the internal machine precision.

% verifies internal machine precision is 2^-52
eps*2^52;

% verifies that only fractional powers of 2 are exact
M = [];
for i=1:50
    M(i) = abs(i*((i+1)/i - 1) - 1);
end
M';
end

function problem2
% problem 2.1: newtoncotes() and clenshawcurtis()

% problem 2.2

f1 = @(x) exp(-x.^2);
f2 = @(x) (1+25*x.^2).^-1;
f3 = @(x) abs(x);

% show convergence
dataf1 = []; dataf2 = []; dataf3 = [];
for N = 1:50
    dataf1(N,:) = [newtoncotes(f1,-1,1,N);
                   clenshawcurtis(f1,-1,1,N);
                   comptrap(f1,-1,1,N)]';
    dataf2(N,:) = [newtoncotes(f2,-1,1,N)
              clenshawcurtis(f2,-1,1,N)
              comptrap(f2,-1,1,N)];
    dataf3(N,:) = [newtoncotes(f3,-1,1,N)
              clenshawcurtis(f3,-1,1,N)
              comptrap(f3,-1,1,N)];
end

% plot values of integral
figure
semilogy(dataf1)
title("integral of f_1 from -1 to 1")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")

figure
semilogy(dataf2)
title("integral of f_2 from -1 to 1")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")

figure
semilogy(dataf3)
title("integral of f_3 from -1 to 1")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")

% NOTE that the negative values for the integral are omitted


% problem 2.3

% on the interval -1 to 1, max|f1'| is small compared 
% to |f2'|, and |f3'| has a discontinuity. Thus, we see 
% terrible errors for newton-cotes, and better errors for 
% clenshaw-curtis. Since none of the functions are highly
% oscillitory, we have standard convergence for composite trapezoid

% problem 2.4

tic;
clenshawcurtis(f1,-1,1,20);
t1 = toc; tic;
comptrap(f1,-1,1,20);
t2 = toc; tic;
quad(f1,-1,1,1e-15);
t3 = toc;

disp("time to compute clenshaw-curtis")
disp(t1)
disp("time to compute composite trapezoid")
disp(t2)
disp("time to compute matlab quad")
disp(t3)

% % error for the above time calculations
% format long
% clenshawcurtis(f1,-1,1,20)-sqrt(pi)*erf(1)
% comptrap(f1,-1,1,20)-sqrt(pi)*erf(1)
% format short

% bonus error plots
E1 = abs(dataf1-sqrt(pi)*erf(1));
E2 = abs(dataf2-2*atan(5)/5);
E3 = abs(dataf3-1);

figure
subplot(1,3,1)
semilogy(E1)
title("error for approximation of f_1")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")
subplot(1,3,2)
semilogy(E2)
title("error for approximation of f_2")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")
subplot(1,3,3)
semilogy(E3)
title("error for approximation of f_3")
legend("newton-cotes", "clenshaw-curtis", "composite trapezoid")

set(gcf,'position',100+0.6*get(0,'ScreenSize'))
end

function I = newtoncotes(f,a,b,N)
% Newton-Cotes quadrature 
% that approximates the integral of a function 
% by the  integral of the degree-N
% polynomial that interpolates it at
% N+1 equispaced points.

% input: f function
%        a,b domain of integration
%        N order of polynomial

x = linspace(a,b,N+1)'; y = f(x);
p = polynomial_interp(x,f)';
I = polyval(polyint(p),b)-polyval(polyint(p),a);
end

function I = clenshawcurtis(f,a,b,N)
% Clenshaw-Curtis quadrature that 
% approximates the integral of a function by the integral 
% of the degree-N polynomial that interpolates it at
% N+1 Chebyshev points.

% inputs: see ncquad

x = cos(pi*(N:-1:0)/N)'; % chebyshev points
x = x*(b-a)/2 + (b+a)/2; % scale points
p = polynomial_interp(x,f)';
I = polyval(polyint(p),b)-polyval(polyint(p),a);
end

function p = polynomial_interp(x,f)
% problem 3.1, hw 5

% Input :
%    x is a vector of size n+1 that contains interpolation points. 
%    f is a handle to the function.
%   Output:
%    p  is a vector of size n+1 that contains the 
%       coefficients of the nth order polynomial 
%       that interpolates f at x.
A = vander(x);
p = A\f(x);
end

function I = comptrap(f,a,b,N)
% computes the composite trapezoid method
% usind N+1 points and N intervals

x = linspace(a,b,N+1)';

% I = 0;
% for i = 1:N
%     I = I + newtoncotes(f,x(i),x(i+1),1);
% end

y = f(x);
h = (b-a)/N;
I = h*(f(a)+f(b)+2*sum(y(2:N)))/2;
end

function problem3

% [val1, N1] = myquad(@(x) 1+atan(30.*x), -1, 1, 1e-6)
% [val1, N1] = myquad(@(x) x+sin(x.^4), 0, 3, 1e-12)
% [val1, N1] = myquad(@(x) exp(-x.^2), -100, 100, 1e-12)

% initialize functions and parameters

f1 = @(x) exp(-x.^2);
f2 = @(x) (1+25*x.^2).^-1;
f3 = @(x) abs(x);
f4 = @(x) 1+atan(30.*x);
f5 = @(x) x+sin(x.^4);
f6 = @(x) x+sin(x.^4);

F = {f1, f2, f3, f4, f5, f6};
A = [-1,1,1e-12]; A = [A; A; A];
A(4:6,:) = [-1, 1, 1e-6; 0, 3, 1e-12; -100, 100, 1e-12];

% compute values of integral, determine number of evaluations
QQ = []; NN = [];
for i = 1:6
    [Q1,N1] = quad(F{i},A(i,1),A(i,2),A(i,3));
    [Q2,N2] = myquad(F{i},A(i,1),A(i,2),A(i,3));
    QQ(i,:) = [Q1 Q2];
    NN(i,:) = [N1 N2];
end

format long
disp("val of quad vs myquad")
disp(QQ)
disp("evals of quad vs myquad")
disp(NN)
format short

% % visualizing the evaluation points
% [~,~,pts] = myquadvis(f1,-1,1,1e-8);
% plot(pts(:,1),pts(:,2),".")

% brief discussion on the myquad algorithm

% we wish to minimize the number of iterations, while also attaining
% a high degree of accuracy. Since we are using a composite method,
% let us choose the simplest closed form method given: simpson's method

% simpson's method is fourth order accurate, but to increase accuracy,
% we also use richardson extrapolation to get a romberg integration 
% method. A single step yeilds an eighth order accurate approximation:

% let ø(h) be simpsons method

% R(j,0) := ø(h/2^j)    j≥0
% R(j,k) := 2^rk R(j,k-1) - R(j-1,k-1)    k≥j≥0
%           ¯¯¯¯¯¯¯ 2^rk - 1 ¯¯¯¯¯¯¯¯¯¯

% R(1,1) := (16R(1,0)-R(0,0))/15

% Then, |R11-R10| is O(h^4). If |R11-R10|<eps, then
% R11 is very likely to be within tolerance, 
% since its error is O(h^8).
% We cannot say for certain though, since if the derivatives of 
% ø(h) near zero are very large, Ω(h^-4), then 
% the estimate may not be in tolerance. However, this is in
% the extreme case. Well behaved functions should converge as normal.
end

function [val,N] = myquad(f,a,b,eps)
% adaptive quadrature

% Input: f   is a handle to a function f(x)
%        a   is the lower bound of integration
%        b   is the upper bound of integration
%        eps is the user-defined precision

% Output: val is the computed value of the integral
%         N   is the number of evaluations of f

c = (a + b)/2;
y = f([a, b, c]);
fa = y(1); fb = y(2); fc = y(3);
[val,N] = iterquad(f,a,b,fa,fc,fb,eps,3);

if N > 1e4
    disp("maximum number of iterations attained")
end
end


function [R11,N] = iterquad(f,a,b,fa,fc,fb,eps,N)
% recursive step

% Input: f   is a handle to a function f(x)
%        a   is the lower bound of integration
%        b   is the upper bound of integration
%        fa  = f(a)
%        fc  = f(c)
%        fb  = f(b)
%        eps is the user-defined precision
%        N   is the number of evaluations of f

% Output: R11 is the computed value of the integral
%         N   is the number of evaluations of f

h = b - a;
c = (a + b)/2;
fd = f((a + c)/2); fe = f((c + b)/2);
N = N + 2;

% simpson's rule one interval
R00 = h*(fa + 4*fc + fb)/6;

% simpson's rule two intervals
R10 = h*(fa + 4*fd + 2*fc + 4*fe + fb)/12;

% romberg integration step
R11 = R10 + (R10 - R00)/15;

if abs(R11-R10) < eps
    % satisfying bound
    return
end

if N > 1e4
    % max iterations
    return
end

% recursive step
[Qac,N] = iterquad(f,a,c,fa,fd,fc,eps,N);
[Qcb,N] = iterquad(f,c,b,fc,fe,fb,eps,N);
R11 = Qac + Qcb;
end


% function [val,N,xy] = myquadvis(f,a,b,eps)
% 
% c = (a + b)/2;
% y = f([a, b, c]);
% fa = y(1); fb = y(2); fc = y(3);
% [val,N,xy] = iterquad2(f,a,b,fa,fc,fb,eps,3,[a fa; b fb; c fc]);
% if N > 1e4
%     disp("maximum number of iterations attained")
% end
% xy = sortrows(xy);
% end
% 
% function [R11,N,xy] = iterquad2(f,a,b,fa,fc,fb,eps,N,xy)
% % recursive step
% 
% h = b - a;
% c = (a + b)/2;
% N = N + 2;
% d = (a + c)/2;
% e = (c + b)/2;
% fd = f(d);
% fe = f(e);
% 
% xy = [xy; d fd; e fe];
% 
% % simpson's rule one interval
% R00 = h*(fa + 4*fc + fb)/6;
% 
% % simpson's rule two intervals
% R10 = h*(fa + 4*fd + 2*fc + 4*fe + fb)/12;
% 
% % romberg integration step
% R11 = R10 + (R10 - R00)/15;
% 
% if abs(R11-R10) < eps
%     % satisfying bound
%     return
% end
% 
% if N > 1e4
%     % max iterations
%     return
% end
% 
% [Qac,N,xy] = iterquad2(f,a,c,fa,fd,fc,eps,N,xy);
% [Qcb,N,xy] = iterquad2(f,c,b,fc,fe,fb,eps,N,xy);
% R11 = Qac + Qcb;
% end
