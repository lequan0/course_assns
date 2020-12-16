function hw7
close all

problem1
problem2
problem3
problem4
end


function problem1
% simple, nonquarantine Kermack-McKentrick epidemic model

c = 1; d= 5;
y0 = [95; 5; 0];

f = @(t,y) [-c*y(1)*y(2); 
             c*y(1)*y(2)-d*y(2);
             d*y(2)];
       
[T,Y] = ode45(f,[0 1],y0);

figure
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3))
title("Kermack-McKentrick epidemic model")
y1 = "healthy people who can get sick";
y2 = "people who are currently sick";
y3 = "people who are either recovered or dead";
legend(y1,y2,y3)
end


function problem2

% problem 2.3
% no quarantine

c1=12.5; c2=4.5; d=.2; r=4.8;
y0 = [0.9999; 0.0001; 0; 0];

Q = @(t) 0;
f = @(t,y) [-(c1*(1-Q(t))+c2*Q(t))*y(1)*y(2); 
    (c1*(1-Q(t))+c2*Q(t))*y(1)*y(2)-(d+r)*y(2);
    r*y(2);      d*y(2)];

[T,Y] = ode45(f,[0 4],y0);

figure
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3),T,Y(:,4))
title("Kermack-McKentrick epidemic model with recovery")
y1 = "susceptible";
y2 = "infected";
y3 = "recovered";
y4 = "dead";
legend(y1,y2,y3,y4)

rmort = Y(end,4)/(Y(end,2)+Y(end,3)+Y(end,4));
disp("mortality rate")
disp(rmort)

% problem 2.4

% death toll, no quarantine, 2024
dtoll = Y(end,4)*33e7;
disp("death toll by 2024, no quarantine")
disp(round(dtoll))

% enforced qurantine
Q = @(t) 1;
f = @(t,y) [-(c1*(1-Q(t))+c2*Q(t))*y(1)*y(2); 
    (c1*(1-Q(t))+c2*Q(t))*y(1)*y(2)-(d+r)*y(2);
    r*y(2);      d*y(2)];

[~,Y] = ode45(f,[0 4],y0);

% death toll, yes quarantine, 2024
dtoll = Y(end,4)*33e7;
disp("death toll by 2024, quarantine")
disp(round(dtoll))
end


function problem3
% problem 3.1
% no quarantine, yes vaccine

c1=12.5; c2=4.5; d=.2; r=4.8;
y0 = [.9999; .0001; 0; 0; 0];

V = @(t) (t<1)*0 + (t>=1)/3;
Q = @(t) 0;
f = @(t,y) [-(c1*(1-Q(t))+c2*Q(t))*y(1)*y(2) - min(V(t),y(1)); 
             (c1*(1-Q(t))+c2*Q(t))*y(1)*y(2)-(d+r)*y(2);
             r*y(2);
             d*y(2);
             min(V(t),y(1))];

[T,Y] = ode45(f,[0 4],y0);

figure
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3),T,Y(:,4),T,Y(:,5))
title("Kermack-McKentrick model with vaccine")
y1 = "susceptible";
y2 = "infected";
y3 = "recovered";
y4 = "dead";
y5 = "vaccinated";
legend(y1,y2,y3,y4,y5)

% death toll, no quarantine, 2024, vaccine
dtoll = Y(end,4)*33e7;
disp("death toll by 2024, Q = 0, vaccine present after 2021")
disp(round(dtoll))

% problem 3.2
f = @D;
a = 1.4; b = 1.5; x = (a+b)/2;
iter = 0;

while abs(f(x)) > 1e-4/6
    if f(a)*f(x) < 0
        b = x;
    else
        a = x;
    end
    x = (a+b)/2; iter = iter+1;
    if iter > 1000
        disp("maximum number of iterations attained")
        x = b;
        break
    end
end
disp(" ")
dtoll = round(f(x)*33e7+100000);
disp("minimum quarantine "+x+" years")
disp("for a death toll of "+dtoll)

xm = floor((x-1)*12);
xd = (x-1)*365-31-28-31-30-31;
disp("or 1 year, "+xm+" months, and "+xd+" days")
disp(" ")

% yp = [];
% xp = 1.4:0.001:1.5;
% for x = xp
%     yp = [yp; f(x)];
% end
% figure
% plot(xp,yp)
end


function dtoll = D(tstar)

% outputs the death rate - 100000/33e7 at t = 2024,
% given that quarantine is enforced before tstar

c1=12.5; c2=4.5; d=.2; r=4.8;
y0 = [.9999; .0001; 0];

V = @(t) (t<1)*0 + (t>=1)/3;
Q = @(t) (t<=tstar) + 0*(t>tstar);
f = @(t,y) [-(c1*(1-Q(t))+c2*Q(t))*y(1)*y(2) - min(V(t),y(1)); 
             (c1*(1-Q(t))+c2*Q(t))*y(1)*y(2)-(d+r)*y(2);
             d*y(2)];

options = odeset('RelTol', 1e-6);
[~,Y] = ode45(f,[0 4],y0,options);

% conditional quarantine, yes vaccine
dtoll = Y(end,3) - 100000/33e7;
end


function problem4
c1=12.5; c2=4.5; d=.2; r=4.8; 
y0 = [.9999; .0001; 0; 0; 0];
% y0 = [y0; .9999*.17; .0001*.17; 0; 0; 0];

V = @(t) (t<1)*0 + (t>=1)/4;
R = @(t,y) 1.50*c1*(y(2) < 0.04);
Q = @(t,y) 0.05 + 0.75*(y(2)/0.2);
f = @(t,y) [-(c1*(1-Q(t,y))+c2*Q(t,y))*y(1)*y(2) - min(V(t),y(1))- R(t,y)*y(1); 
             (c1*(1-Q(t,y))+c2*Q(t,y))*y(1)*y(2)-(d+r)*y(2) + R(t,y)*y(1);
             r*y(2);
             d*y(2);
             min(V(t),y(1))];

[T,Y] = ode45(f,[0 4],y0);
% Y = Y(:,1:5) + Y(:,6:10);
figure
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3),T,Y(:,4),T,Y(:,5))
title("Kermack-McKentrick model, travel and conditional quarantine")
y1 = "susceptible";
y2 = "infected";
y3 = "recovered";
y4 = "dead";
y5 = "vaccinated";
legend(y1,y2,y3,y4,y5)


% death toll, conditional quarantine, 2024, vaccine
dtoll = Y(end,4)*33e7;
disp("death toll by 2024, conditional Q, travel")
disp(round(dtoll))

end