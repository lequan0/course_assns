function hw4

close all

disp("please note that the error is given in terms of the standard euclidean norm, as opposed to its square")

% % test cases for LLS_SVD
% m = 4;
% n = 3;
% A = rand(m,n);
% b = rand(m,1);
% A\b
% LLS_SVD(A,b)


% % import data from website
data = [62 41 3.68;	 
        67 44 2.98;	 
        73 51 3.36;	 
        79 57 3.60;	 
        86 66 5.15;	 
        91 72 5.35;	 
        94 74 3.18;	 
        94 73 3.83;	 
        89 68 4.33;	 
        82 59 4.50;	 
        72 50 4.19;	 
        65 43 3.69];
datat = data(:,1); % we consider monthly average high
datap = data(:,3);

% % problem 5.2

% % find c0 and c1, and find least squares regression 
% % for 1 and cos
a = 2*pi/365; % angular frequency 
days = [31 28 31 30 31 30 31 31 30 31 30 31];
yp = [];
b = 0;
for i = 1:12
    c = b + days(i);
    yp((b+1):c) = datat(i); % construct daily data from monthly
    b = c;
end
xp = 1:365; % plotting points
T = [cos(a.*xp)];
T = [ones(365,1) T'];
c = LLS_SVD(T,yp');
f = @(d) c(1) + c(2)*cos(a.*d);
figure
plot(xp,f(xp),xp,yp);
title("two term temperature LLS")
disp("Error for two term temp LLS: " + norm(T*c-yp'))

% % problem 5.3
% we attempt to use an additional sin function to
% improve regression

% % find least squares regression for 1, cos, sin
T = [T sin(a.*xp')];
c = LLS_SVD(T,yp');
g = @(d) c(1) + c(2)*cos(a.*d) + c(3)*sin(a.*d);
figure
plot(xp,g(xp),xp,yp);
title("three term temperature LLS")
disp("Error for three term temp LLS: " + norm(T*c-yp'))

% % problem 5.4

% % repeat for the precipitation data
a = 2*a; % lets try double the angular frequency
days = [31 28 31 30 31 30 31 31 30 31 30 31];
yp = [];
b = 0;
for i = 1:12
    c = b + days(i);
    yp((b+1):c) = datap(i); % construct daily data
    b = c;
end
T = [sin(a.*xp)];
T = [ones(365,1) T'];
c = LLS_SVD(T,yp');
h = @(d) c(1) + c(2)*sin(a.*d);
figure
plot(xp,h(xp),xp,yp);
title("two precipitation temp LLS")
disp("Error for two term precipitation LLS: " + norm(T*c-yp'))
end


function x = LLS_SVD(A,b)
% linear least squares solver
% Input:
%       A   is a m x n matrix
%       b   is a m vector
% Output:
%       x   is the solution to the LLS problem


[U, S, V] = svd(A);
I = find(S~= 0);
S(I) = 1./S(I);
x = V*S'*U'*b;
end

% % using a monthly averave scheme
% yp = [];
% a = 2*pi/365;
% days = [31 28 31 30 31 30 31 31 30 31 30 31];
% b = 0;
% for i = 1:12
%     c = b + days(i);
%     T(i) = mean(cos(a.*((b+1):c)));
%     yp((b+1):c) = datat(i);
%     b = c;
% end
% T = [ones(12,1) T'];
% c = LLS_SVD(T,datat);
% f = @(d) c(1) + c(2)*cos(a.*d);
% 
% xp = 1:365;
% plot(xp,f(xp),xp,yp);
