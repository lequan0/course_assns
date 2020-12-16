function kaprekar
% kaprekar
% Generates all kaprekar paths for two sets of data, plots them
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 03
% Last Modified: October 6, 2019

% generate kaprekar paths for all valid four digit integers
x = 4; maxiter = 8;
data1 = quest(x,maxiter);
% plot kaprekar paths
figure(1)
plot(0:8, data1, 'o-')
title('kaprekar paths for all 4-digit integers')
xlabel('iterations')
ylabel('value')

% % generate kaprekar paths for all valid five digit integers
x = 5; maxiter = 12;
data2 = quest(x,maxiter);
% plot kaprekar paths
figure(2)
plot(0:12, data2, 'o-')
title('kaprekar paths for all 5-digit integers')
xlabel('iterations')
ylabel('value')
end

function d = dis(x)
% d = dis(x)
% Returns the first invalid integer for a given digit length
% Inputs: x, the number of digits
% Outputs: d, the first invalid
str = '';
for digit = 1:x
    str = strcat(str, '1');
end
d = str2double(str);
end

function data = quest(x,maxiter)
% data = quest(x,maxiter)
% Generates all kaprekar paths for a given integer length
% Inputs: x, the number of digits
%         maxiter, maximum iterations
% Outputs: data, a matrix of all the kaprekar paths

% Set a range of integers to operate on
low = 10^(x-1);
range = setdiff( low:9*dis(x), (1:9)*dis(x) );
% Generate the kaprekar paths
data = [];
for num = range
    % Initialize a kaprekar path for a given num
    y(1) = num;
    % Generate a kaprekar path
    for n = 1:maxiter
        if y(n) > low-1
            str1 = sort(num2str(y(n)));
        else % Must be done for integers like 0999
            str1= strcat('0'+ sort(num2str(y(n))));
        end
        y(n+1) = str2double(flip(str1))-str2double(str1);
    end
    % Concatenate kaprekar path
    data = [data; y];
end
end

% Remarks on the data from the 4-digit kaprekar paths:
% All converge to 6174

% Remarks on the data from the 5-digit kaprekar paths:
% There seems to be no kaprekar constants, and there are 3 kaprekar holes:
% Two with four elelments each, and one with two elements