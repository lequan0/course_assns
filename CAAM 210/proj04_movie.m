function moviedriver
% moviedriver
% creates a movie showing the time dependence of the predator prey model
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 04
% Last Modified: October 14, 2019


[timerange,r,f] = ode_method; % generates data
length(timerange)

% initaliizes video
vid = VideoWriter('proj04 video','MPEG-4'); 
open(vid) % begins writing
for index = 1:(length(timerange)) % for each data point
    clf % ensures graphs don't layer into frames
    hold on
    time = timerange(index);
    plot(r,f) % graph
    plot(r(index),f(index),'o') % point
    
    title("time: " + string(time))
    xlabel("Rabbits")
    ylabel("Foxes")
    
    hold off
    Frame = getframe(gcf);
    writeVideo(vid,Frame); 
end

close(vid)

end



function [t,r,f] = ode_method % predator-prey system, ode45
% [t,r,f] = ode_method
% Uses ode45 to find the solution to the predator-prey system, based on
% pre-specified initial conditions
% Inputs: None
% Outputs: t, the vector domain of the solution generated, correspods to 
%             the time sample
%          r, the vector of rabbit populations, with respect to time
%          f, the vector of fox populations, with respect to time


% initialize constants
k1=3; k2=3e-3; k3=6e-4; k4=0.5; 
tspan = [0 15];
% initial conditions
r0 = 1000;
f0 = 500;
i0 = [r0;f0];
% define the system of equations
odefun = @(t,Y) [k1*Y(1)-k2*Y(1)*Y(2);
                k3*Y(1)*Y(2)-k4*Y(2)];

% note that we will let Y(1) equal to the rabbit pop, with the Y(2) being
% the fox pop

options = odeset('RelTol', 1e-6); % precision of model
[t, rf] = ode45(odefun,tspan,i0,options); % generates the actual data
% splits the data set to r and f
r = rf(:,1);    f = rf(:,2);
size(t)
end