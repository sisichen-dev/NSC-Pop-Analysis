function [grad] = grad_beta(dist, input, beta)
% compute gradient of log likelihood given 
% dist: sample distribution
% input: input condition 
% beta: current value of stacks of beta vectors 
    cur_dist = model_dist(beta, input);
    baye_scale = sum(dist, 1) / sum(dist(:)) * size(dist, 2);
    grad = (size(dist, 2) * dist / sum(dist(:))  - cur_dist .* baye_scale)  * input';
    grad(end, :) = 0; % the last row is set as 0 for uniquely determine model parameters
end

