function proj07
% proj07
% uses the metropolis algorithm to decode a simple cypher
% Inputs: None
% Outputs: None
% Quan Le, CAAM 210, Fall 2019, Project 07
% Last Modified: October 29, 2019


% import strings
string1 = fileread('encodedtext1.txt');
string2 = fileread('encodedtext2.txt');
string3 = fileread('encodedtext3.txt');

% run the decoder on each encoded message
disp("decoded data" + newline)
disp("decoded encodedtext1.txt")
decoder(string1, 10000)
disp("decoded encodedtext2.txt")
decoder(string2, 10000)
disp("decoded encodedtext3.txt")
decoder(string3, 10000)
disp(newline)
end



function [message] = decoder(string,iter)
% [message] = decoder(string,iter)
% decodes an encoded message using a key mapping from the metropolis
% algorithm
% Inputs: string, an encoded message string
%         iter, the number of times to run the metropolis alg
% Outputs: message, the decoded message string


% generate the key using the metropolis algorithm
y = Metropolis(string,iter);

% construct the message
message = "";
len = length(string);
for char = string(1:(len-1)) % note that len-1 is used to exclude " ' "
    yval = downlow(double(char));
    add = downlowinv(y(yval));
    message = message + add;
end

% replace ` with spaces
message = strrep(message,'`',' ');
end



function [y] = Metropolis(string,trials)
% [y] = Metropolis(string,trials)
% generates a key for a simple cypher encoded message
% Inputs: string, an encoded message string
%         trials, the number of times to run the metropolis algorithm
% Outputs: y, the cypher key


% generate a completely random y
y = randperm(27);

% metropolis algotithm
for trial = 1:trials
    ymaybe = guess(y); % guess is a helper function, swaps two values in y
    proby = loglike(y, string);
    probymaybe = loglike(ymaybe, string);
    if probymaybe > proby
        y = ymaybe;
    else
        r = rand;
        threshold = exp(-proby+probymaybe);
        if r < threshold
            y = ymaybe;
        end
    end
end
end



function [ymaybe]=guess(y)
% [ymaybe]=guess(y)
% randomly swaps two elements in a cypher key, y
% Inputs: y, the initial cypher key
% Outputs: ymaybe, the edited cypher key


% swaps
ymaybe = y;
swap = randsample(27,2);
swap1 = y(swap(1));
swap2 = y(swap(2));
ymaybe(swap(1)) = swap2;
ymaybe(swap(2)) = swap1;
end

function [prob] = loglike(y, string)
% [prob] = loglike(y, string)
% finds the log-likelihood for a given y and string
% Inputs: y, the prospective cypher key
%         string, the encoded message
% Outputs: prob, the log-likelihood of y and string


% load the letter transition probabilities
load("letterprob.mat")

% generate total probabilites
prob = 0;
for idx = 1:(length(string)-2) % iteration over character pairs
    char1 = string(idx);
    i = y(downlow(double(char1))); % conversion to code27
    char2 = string(idx+1);
    j = y(downlow(double(char2)));
    prob = prob+ log(letterprob(i,j));
end
end



function [code27] = downlow(ascii_num)
% [code27] = downlow(ascii_num)
% given an ascii code, it returns its equivalent in 1:27
% Inputs: ascii_num, the ascii representation of the character
% Outputs: code27, the representation of the character as an index in 1:27

code27 = ascii_num-95;
end



function [char_want] = downlowinv(code27)
% [char_want] = downlowinv(code27)
% given the code27 representation, it returns the character
% Inputs: code27, the representation of the character as an index in 1:27
% Outputs: char_want, the actual character

char_want = char(code27+95);
end