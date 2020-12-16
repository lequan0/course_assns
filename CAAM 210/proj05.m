function moneymaker
% moneymaker
% generates the outcomes for 2000 trials of a 40yr investment strategy
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 05
% Last Modified: October 18, 2019


% initials
A = [];
trials = 2000;
% implementing strategy, generating data
for trial = 1:trials
    balance = 0; % initialzes balance for a new year
    yeartrial = []; % initializes balance vector
    for year = 1:40 % note on this in the strategy comments **** 
        balance = balance + 10000;
        distro = strategy(balance,year);
        arg = (balance).*distro; % prepares arguements for oneyear
        balance = oneyear(arg(1), arg(2), arg(3), arg(4), arg(5));
        % where arg(1...5) are: stocks, bonds, bank, inverted, casino
        yeartrial = [yeartrial; balance];
    end
    A = [A yeartrial];
end
plotter(A) % plotting first 20 trials
summary(A,trials) % statistical characteristics of summary

% explaination of summary

% The primary goal of this investment strategy was to maximize the amount
% of final returns. To this end, different investment strategies were 
% implemented based on the expected returns for each investment option. The
% following expected values were calculated:

% stocks: 0.0540
% bonds: 0.0311
% bank: 0.9994
% inverted: 0.0217
% casino: -0.05

% Clearly, stocks have the highest expected value, and thus, all money was
% invested in stocks became the default strategy. However, the proportion
% above 1.5m was only 0.33. To increase this number, we must consider the
% two extreme cases, banks, and casinos.

% If the balance can, in fact, reach 1.5m in 1% intervals for the years
% remaining, then we can simply place the entirety of the balance into the
% bank to reduce loss.

% As for the casino, if a value can only reach the target in 60% intervals,
% it is incerdibly unlikely to have such a jump for each year in stocks.
% Thus, we can instead place it all in the casino to increase the number of
% trials that do reach target.

% The same logic applies to switching to bonds, this is just done to reduce
% risk in later years.

% **** afternote
% note that 1:40 is used instead of 0:39 to compensate for not including 
% the 10000$ increases in the balance each year for the thresholds 

% tbh this was originally an error, but it ended up helping lol
end

function [distro]=strategy(balance,year)
% [distro]=strategy(balance,idx)
% determines the strategy to be used
% Inputs: balance, the amount of money to be invested
%         year, the year in which the investment is done
% Outputs: distro, the strategy to be used


if balance >= (1500000/(1.01^(40-year))) % safest strategy, for bank
    distro = [0 0 1 0 0];
elseif balance >= (1500000/(1.025^(40-year))) % safe strategy, for bonds
    distro = [0 1 0 0 0];
elseif balance <= (1500000/(1.80^(40-year))) % risky strategy, for casino
    distro = [0 0 0 0 1];
else
    distro = [1 0 0 0 0]; % default strategy
end
end

function plotter(A)
% plotter(A)
% plots the first 20 columns in A
% Inputs: A, a 40x2000 matrix, where the columns are the trials and the
%            rows are time dependent behavior
% Outputs: None


clf % clears the old graph
hold on
B = A(:, 1:20); % grabs the first 20 columns
plot(B)
title("balances over time")
xlabel("years")
ylabel("monetary balance")
% plot(0:40, ones(41,1).*1.5e6) % plots goal line
axis tight
hold off
end

function summary(A,trials)
% summary(A,trials)
% prints a statistical summary of the trials run
% Inputs: A, a 40x2000 matrix, where the columns are the trials and the
%            rows are time dependent behavior
%         trials, the number of trials run for each 40-year interval
% Outputs: None


E = A(40,:); % grabs the final values
% finds the proportion of final values over 1500000
madeit = length(find(E>1500000))/trials; 
% actual summary
disp("Trials: " + string(trials))
disp("Maximum Ending Value: " + string(max(E)))
disp("Median Ending Value: " + string(median(E)))
disp("Mean Ending Value: " + string(mean(E)))
disp("Minimum Ending Value: " + string(min(E)))
disp("Proportion above 1500000: " + string(madeit))
disp(" ")
end

function [balance] = oneyear(stocks, bonds, bank, inverted, casino)
% [balance] = oneyear(stocks, bonds, bank, inverted, casino)
% generates a new balance after a year of investments
% Inputs: stocks, the amount of money invested into stocks
%         bonds, the amount of money invested into bonds
%         bank, the amount of money invested into the bank
%         inverted, the amount of money invested into inverted stocks
%         casino, the amount of money invested into the casino
% Outputs: balance, balance generated from weighted probabilities, sum of
%                   all the resultant investments 


% generate randoms
yrtype = rand(); 
outcometype = rand(1,5); % assuming independent investment options
rowp = [.15 .35 .45 .05]; % row probabilities
colp = [.2 .4 .2 .2]; % column probabilities
% determines type of year based on cumulative probabilities
cumprob1 = 0;
for i = 1:4
    cumprob1 = cumprob1 + rowp(i);
    if yrtype < cumprob1
        rowtype = i; % rowtype
        break
    end
end
% determines type of result per stock based on cumulative probabilities
% note that this type must be independent for each stock
coltype = [];
for num = outcometype
    cumprob2 = 0;
    for j = 1:4
        cumprob2 = cumprob2 + colp(j);
        if num < cumprob2
            newtype = j;
            break
        end
    end
    coltype = [coltype newtype]; % column type
end
% probabilities for stocks, bonds, inverted
% note that the rows are the returns for a given year
stocksM = [.08  .16  .24  .50;
          -.04  .04  .08  .12;
          -.04  .00  .04  .10;
          -.40 -.08 -.04  .00];

bondsM =  [.04  .06  .07  .10;
          -.02  .03  .06  .09;
          -.02  .03  .05  .08;
          -.25 -.01 -.50  .03];

istocksM = [-.15  .05  .10  .15;
            -.15  .05  .10  .15;
            -.05 -.02  .00  .02;
            -.06  .00  .06  .8];
% generates returned values for the three investment options above
stocks = stocks * (1 + stocksM(rowtype,coltype(1)));
bonds = bonds * (1 + bondsM(rowtype,coltype(2)));
inverted = inverted * (1 + istocksM(rowtype,coltype(3)));
% generates returned value for bank
if rand() <= 0.9999
    bank = bank * (1.01);
else
    bank = bank * (.5);
end
% generates returned value for casino
if rand() <= .5 
    casino = casino * (1.9);
else
    casino = 0;
end
% generates final balance
balance = stocks + bonds + inverted + bank + casino;
end
