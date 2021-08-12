function [means, stds]=errorplot(x,y, uniquex, varargin)
% function errorplot(x,y,varargin)
% X and Y are vectors of the same length; This function plots error bars
% on Y, binned by the values in X. 
% X is a categorical variable expressed as integers.  
% Y is a numerical measurement conditional on X. 
% This function bins X into different integers, calculates a mean and
% standard deviation for each bin, and plots an error bar at that location.
% 

means = [];
stds = [];

if nargin==2
    uniquex = unique(x,'stable');
end
    
for i = 1:length(uniquex)
    currx = uniquex(i);
    means = [means,mean(y(x==currx))];
    stds = [stds,std(y(x==currx))];  
end

confplot(1:length(uniquex), means, stds, stds, varargin{:})
%errorbar(uniquex, means, stds, varargin{:}) % 'k','LineWidth',1.5,'linestyle','none'); 