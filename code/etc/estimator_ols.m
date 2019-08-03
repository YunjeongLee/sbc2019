function [theta, cost_val, cal_time] = estimator_ols(data, cost_ols, params, theta0, lb, ub)
if(nargin < 5)
    lb = [];
    ub = [];
end

[theta, cost_val, cal_time] = minimizer(data, cost_ols, params, theta0, lb, ub);