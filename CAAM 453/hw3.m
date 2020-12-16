function hw3
close all
% problem 1: newton polynomial interpolation

% % test cases

% f = @(x) abs(x);
% % f = @(x) cos(x);
% xn = linspace(-5,5,10);
% [p] = newtonpoly(xn',f(xn)');
% xp = linspace(-5, 5, 200);
% xpp = linspace(-5, 6, 200);
% % plot(xp, f(xp), xp, polyval(p,xp))
% pnew = newtoncont(p,xn,6,f(6));
% plot(xp,f(xp),xp,hornersmethod(p,xp),xpp,hornersmethod(pnew,xpp))


figure(1)
xp = linspace(.1, 4, 1000);     % plotting points
hold on
plot(xp, gamma(xp),'k')
n = [5 10 15];                  % number of interpolation points
shape = [":"; "--"; "-."];
for i = 1:3
    x = linspace(.1, 4, n(i))';
    p = newtonpoly(x,gamma(x));
    plot(xp, hornersmethod(p, xp),shape(i))
end
title("polynomial interpolation: newton basis")
legend("\gamma function","p_5 interpolant","p_{10} interpolant","p_{15} interpolant")
hold off


figure(2)
hold on
plot(xp, gamma(xp),'k')
x = linspace(.5, 3.5, 10)';     % original set of interpolation pts
p = newtonpoly(x,gamma(x));
plot(xp, hornersmethod(p, xp),'--')     % plot first newton poly
pnew = newtoncont(p,x,.1,gamma(.1));    % adds (.1,gamma(.1))
x = [x; .1];
pnew = newtoncont(pnew,x,4,gamma(4));   % adds (4,gamma(4))
plot(xp, hornersmethod(pnew, xp),'-.')  % plots updated newton poly
title("polynomial interpolation: updating newton basis")
legend("\gamma function","p_{10} interpolant","p_{12} interpolant")
hold off


% problem 3: piecewise cubic Hermite polynomials
% computes f, f', f''
syms x
f = 1/(1+x^2); df = diff(f); ddf = diff(df);
f = matlabFunction(f); 
df = matlabFunction(df);
ddf = matlabFunction(ddf);
F = {f, df, ddf};           % collects f, f', f''

xp = linspace(-4, 4, 1000); % plotting points
N = [3 4 8];                % number of intervals
C = [0.4 1.0 0.4;           % plotting colors
     1.0 0.4 0.4;
     0.4 0.4 1.0];
T = ["piecewise interpolation of f=1/(1+x^2)";
     "first derivative of piecewise interpolant";
     "second derivative of piecewise interpolant"];
L = ["3 intervals", "5 intervals", "8 intervals"]; % titles
L = [L; L; L];
L(:,4) = ["f","f_{x}","f_{xx}"]'; % legend



for i = 1:3                 % for each number of intervals
    Fp = [];
    for j = 1:N(i)          % for each interval
        a=(j-1)*8/N(i) - 4; b=j*8/N(i) - 4;     % endpoints
        p = hermite(f,df,a,b);
        dp = polyder(p);
        ddp = polyder(dp);
        x = linspace(a, b, 8000/N(i));      % plotting domain
        A = [x; hornersmethod(p, x); hornersmethod(dp, x); hornersmethod(ddp, x)];
        Fp = [Fp A];                        % plotting points
    end
    for k = 1:3             % for each derivative of p
        figure(2+k)
        hold on
        plot(Fp(1,:),Fp(k+1,:),'color',C(i,:)) 
        hold off
    end
end
for i = 1:3                 % plot f,f',f''
    figure(2+i)
    g = cell2mat(F(i));
    hold on
    plot(xp, g(xp),'k')
    hold off
    title(T(i))
    legend(L(i,:))
end

% computes max error for each N
y = [];
for N = 1:100
    ymax = [];
    for j = 1:N
        a=(j-1)*8/N - 4; b=j*8/N - 4;
        p = hermite(f,df,a,b);
        x = linspace(a, b, 8000/N);
        ymax = [ymax max(abs(polyval(p,x)-f(x)))]; % max of each interval
    end
    y = [y max(ymax)]; % max over all the intervals
end
figure(6)
semilogy(y)
title("error bound on |H_N(x)-f(x)|")

% comments: we see convergence as expected, with slight variation
% between even and odd numbered intervals. 
end


function p = hermite(f,df,a,b)
% f, df: function and its derivative
% a, b: endpoints of interval
p = f(a)*Aj(a,b) + f(b)*Aj(b,a);        % match endpoints
p = df(a)*Bj(a,b) + df(b)*Bj(b,a) + p;  % match derivatives
end

function [p] = Aj(xj,xk)
% applies cubic hermite polynomials for continuity of p
p = conv([1 -xk]',[1 -xk]');
p = conv(p, [-2; 3*xj-xk]);
p = p./((xj-xk)^3);
end

function [p] = Bj(xj,xk)
% applies cubic hermite polynomials for continuity of p'
p = conv([1 -xk]',[1 -xk]');
p = conv(p,[1 -xj]');
p = p./((xj-xk)^2);
end

function [p] = newtonpoly(x,y)
% generates interpolating polynomial using newton basis
% input:    x, y: vectors of interpolation points
% output:   p: interpolating polynomial

n = length(x)-1;

% construct system to find coefficients of newton basis
b = ones(n+1,1);
A(:,1) = b;
for i = 1:n
     b = b .* (x-x(i));
     A(:,i+1) = b;
end
c = A\y;     % solve for coefficients

% iteratively create interpolating polynomial
q = 1;
p = q*c(1);
for j = 1:n
    q = conv(q,[1 -x(j)]);      % newton basis
    p = polyadd(p,c(j+1)*q');   % construct the polynomial
end
end

function w = polyadd(u, v)
% adds coefficient vectors of polynomials
% u, v are column vectors, w is their sum

n1 = length(u);
n2 = length(v);

% matches length
if n1 < n2
    u = [zeros(n2-n1,1); u];
elseif n1 > n2
    v = [zeros(n1-n2,1); v];
end
w = u+v;  
end

function b = hornersmethod(p, x)
% evaluates polynomial
% input:    p: vec, polynomial coefficients
%           x: plotting domain
% output:   b: p evaluated at x

n = length(p)-1;
b = p(1);
for i = 1:n
    b = p(1+i) + b.*x; % horner's method step
end
end

function pnew = newtoncont(p,x,xnew,ynew)
% updates newton poly with new point
% input:    p: old newton poly
%           x: set of interpolation points
%           xnew, ynew: new data point
% output:   pnew: new newton poly

q = 1;
for j = 1:length(p)
    q = conv(q,[1 -x(j)]);  % newton basis function q_{n+1}
end
c = ynew-polyval(p,xnew);   % evaluates c_{n+1}  
c = c/polyval(q,xnew);
pnew = polyadd(p,c*q');     % updates p
end