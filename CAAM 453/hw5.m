function hw5
% driver, no inputs/outputs
close all
% problem1
% problem2
problem3
end

function problem1
% problem 1.9
% plots the legendre polynmials p_0 to p_5

N = 6; % number of polynomials to generate
p = zeros(N,N);
p(1,end) = 1; % p0 = 1
p(2,end-1) = 1; 

% generating legendre polynomials using 
% bonnet's recursion formula
for n = 1:(N-2) 
    q = (2*n+1)*conv([1 0],p(n+1,:))/(n+1); % recursive step
    q = q(2:end);
    p(n+2,:) = q-(n)*p(n,:)/(n+1);
end

% plotting legendre polynomials for N=0 to 5
x = -1:0.01:1;
figure
hold on 
for n = 1:min(6,N)
    plot(x,polyval(p(n,:),x))
end
hold off
title("legendre polynomials of order n=0,..., 5")

% % verify result for problem 1.7
% I = [];
% for i = 1:N
%     for j = 1:N
%         r = polyder(p(i,:));
%         t = conv(r,p(j,:));
%         f = @(x) polyval(t,x);
%         I(i,j) = integral(f,-1,1);
%     end
% end
% int8(I) 


% % plotting P4'*P3
% r = polyder(p(4,:));
% t = conv(r,p(3,:));
% f = @(x) polyval(t,x);
% I = integral(f,-1,1);
% plot(x, f(x),x,polyval(r,x),x,polyval(p(3,:),x),x,zeros)
end

function problem2
% radius
a = 1; 
[X, Y] = meshgrid(-linspace(-a,a,600)); % initialize plotting points
X = flip(X')';
[Q,R] = cart2pol(X,Y); % convert to polar coords

figure
N = [2 10 100];
for i = 1:length(N)
    T = harmonic(R,Q,N(i),a);
    T = T - T.*(R > a); % zero out terms outside domain
    
    % color plot
    subplot(length(N),3,3*i-2)
    imagesc(flip(T'))
    title("color plot of T, N = "+ N(i))
    colormap jet
    colorbar
    axis square
    
    % verify boundary conditions
    subplot(length(N),3,3*i-1)
    q = 0:0.1:pi;
    plot(q,harmonic(1,q,N(i),a))
    title("T(1,\theta)  N = "+ N(i))
    axis square
    
    subplot(length(N),3,3*i)
%     T = T - (T).*(T > 1.5);
%     T = T - (T).*(T < -.5);
    mesh(T)
    title("mesh plot of T, N = "+ N(i))
    
    Tc = harmonic(0,0,N(i),a);
    disp("N = "+N(i)+", T(0,0) = "+Tc)
end

% determine figure size
set(gcf,'position',100+0.6*get(0,'ScreenSize'))

% % testing validity of intP
% N = 10;
% p = zeros(N,N);
% p(1,end) = 1;
% p(2,end-1) = 1;
% % generate legendre polynimials
% for n = 1:(N-2) 
%     q = (2*n+1)*conv([1 0],p(n+1,:))/(n+1); % recursive step
%     q = q(2:end);
%     p(n+2,:) = q-(n)*p(n,:)/(n+1);
% end
% I = zeros(N+1,2);
% for n = 1:N
%     f = @(x) polyval(p(n,:),x);
%     I(n,:) = [integral(f,0,1) intP(n-1)];
% end
% I
end

function T = harmonic(r,q,N,a)
% computes an N-th order approimation of the 
% temperature at radius r, degree q

T = 0;

P{1} = 1; % P{k} is the k-1 degree legendre polynomial
P{2} = [1 0];

for n = 0:N
    m = n+1;
    if m > 2 % recursively generates P_n
        P{m} = polyadd((2*n-1)*conv([1 0], P{m-1}), -(n-1)*P{m-2})./(n);
    end
    An = intP(n); 
%     f = @(x) (1-x).*polyval(P{m},x);
%     An = integral(f,-1,1);  % alternate boundary conditions
    An = An*(2*n+1)/(2*a^n);
    T = T + An.*r.^(n).*polyval(P{m}, cos(q));
end


end

function y = intP(n)
% computes the integral of an n-th degree
% Legendre polynomial on the interval (-1, 1)

if n == 0
    y = 1;
elseif mod(n,2) == 0
    y = 0;
else
    y = 0;
    for j = 0:(n-1)/2
        y = y + (4*j+1)*(-1)^j*factorial(2*j)/(2^(2*j)*factorial(j)^2);
    end
    y = y/(n^2+n);
end
end

function w = polyadd(u, v)
% adds coefficient vectors of polynomials
% u, v are row vectors, w is their sum

n1 = length(u);
n2 = length(v);

% matches length
if n1 < n2
    u = [zeros(n2-n1,1)' u];
elseif n1 > n2
    v = [zeros(n1-n2,1)' v];
end
w = u+v;  
end


function problem3
% =============================================
% problem 3.1
% =============================================

for n = [11 21]
    f = @(x) 1./(1+x.^2);
    
    % equally spaced points
    x = linspace(-5,5,n+1)';
    p = polynomial_interp(x,f);
    xp = linspace(-5,5,500);
    yp1 = polyval(p,xp);
    
    % chebyshev roots
    x = 5*cos((2*(1:(n+1))-1)*pi/(2*(n+1)))';
    p = polynomial_interp(x,f);
    yp2 = polyval(p,xp);
    
    % legendre roots
    syms x
    x = vpasolve(legendreP(n-1,x) == 0);
    x = 5*[-1; x; 1];
    p = polynomial_interp(x,f);
    p = double(p);
    yp3 = polyval(p,xp);
    
    % plotting
    figure
    plot(xp,f(xp),xp,yp1,xp,yp2,xp,yp3)
    title("order "+n+" polynomial interpolation")
    legend("f(x)", "equally spaced", "Chebyshev", "Legendre")
end

% =============================================
% problem 3.2
% =============================================

E = zeros(21,3);
for n = 1:21
    
    % equally spaced points
    x = linspace(-5,5,n+1)';
    p = polynomial_interp(x,f);
    xp = linspace(-5,5,500);
    yp1 = polyval(p,xp);
    
    % chebyshev roots
    x = 5*cos((2*(1:(n+1))-1)*pi/(2*(n+1)))';
    p = polynomial_interp(x,f);
    yp2 = polyval(p,xp);
    
    % legendre roots
    syms x
    if n ~= 1
        x = vpasolve(legendreP(n-1,x) == 0);
    else x = [];
    end
    x = 5*[-1; x; 1];
    p = polynomial_interp(x,f);
    p = double(p);
    yp3 = polyval(p,xp);
    
    E(n,:) = [norm(f(xp)-yp1,Inf), 
              norm(f(xp)-yp2,Inf), 
              norm(f(xp)-yp3,Inf)]; 
end
figure
semilogy(E)
title("error for order "+n+" polynomial interpolation")
legend("equally spaced", "Chebyshev", "Legendre")

% =============================================
% problem 3.4
% =============================================

% If n goes to infinity, it's likely that 
% the interpolations on the Chebyshev and Legendre roots
% will continue to converge, while that on the equally 
% spaced points will diverge. The behavior of the Legendre 
% and Chebyshev point interpolations have been shown to converge,
% while that of the equally spaced points seems that way
% from the generated error bounds. 

% This hints that the choice of interpolation points is quite 
% significant when trying to produce viable interpolation polynomials; 
% as the interpolation is strongly variant based off of such a choice.

end

function p = polynomial_interp(x,f)
% problem 3.1

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
