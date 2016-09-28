function s = randstr(N)
%s = RANDSTR(N)
%
%This function generates a random string of length N.

symbols = ['a':'z' 'A':'Z' '0':'9'];
nums = randi(numel(symbols), 1, N);
s = symbols (nums);
