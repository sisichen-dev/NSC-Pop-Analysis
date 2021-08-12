function [dist] = model_dist(beta, input)
% compute the distribution of the log linear model given 
% beta: current value of stacks of beta vectors 
% input: input condition 
    score = beta * input;
    dist = exp(score) ./ sum(exp(score), 1);
end