function [rows, cols] = findsubplotsize(n)
% function [rows, cols] = findsubplotsize(n)
% finds the optimal number of rows and columsn for a given number of plot
% elements

rows = 1;
cols = 1;
switcher = 0;

while (rows*cols < n)
    switcher = switcher + 1;
    if mod(switcher,2)==0
        rows = rows+1;
    else
        cols = cols+1;
    end
end