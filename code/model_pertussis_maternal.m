function dydt = model_pertussis_maternal(t, y, params)
%% Assign Parameter
names = params(:,1);
for k = 1:size(params,1)
    cmd_string = sprintf('%s_ = params{%d, 2};', names{k}, k);
    eval(cmd_string);
end

%% Generate ODE system
% S : y(1:num_grps_)
% P1_AI : y(num_grps_+1:2*num_grps_)
% P2_AI : y(2*num_grps_+1:3*num_grps_)
% P3_AI : y(3*num_grps_+1:4*num_grps_)
% C_AI : y(4*num_grps_+1:5*num_grps_)
% I1 : y(5*num_grps_+1:6*num_grps_)
% I2 : y(6*num_grps_+1:7*num_grps_)
% I3 : y(7*num_grps_+1:8*num_grps_)
% R : y(8*num_grps_+1:9*num_grps_)
dydt = zeros(9 * num_grps_ + 1, 1);

Lambda = contact_matrix_ * (y(5*num_grps_+1:6*num_grps_) ...
    + rho1_ * y(6*num_grps_+1:7*num_grps_) + rho2_ * y(7*num_grps_+1:8*num_grps_));

%% Generate susceptible
dydt(1:num_grps_) = ...
    - Lambda .* y(1:num_grps_) ...
    + sigma0_ * y(num_grps_+1:2*num_grps_) ...
    - mu_ .* y(1:num_grps_) ...
    + [0; c_] .* (1 - vacc_) .* [0; y(1:num_grps_-1)] ...
    - [c_; 0] .* y(1:num_grps_) + [B_; zeros(num_grps_-1, 1)];

%% Generate partially vaccinated group 1 (P1_AI)
dydt(num_grps_+1:2*num_grps_) = ...
    - Lambda .* y(num_grps_+1:2*num_grps_) ...
    - sigma0_ * y(num_grps_+1:2*num_grps_) ...
    + tauP_ * y(2*num_grps_+1:3*num_grps_) ...
    - mu_ .* y(num_grps_+1:2*num_grps_) ...
    + [0; c_] .* (1 - vacc_) .* [0; y(num_grps_+1:2*num_grps_-1)] ...
    - [c_; 0] .* y(num_grps_+1:2*num_grps_) ...
    + [0; c_] .* vacc_ .* [0; y(1:num_grps_-1)];

%% Generate partially vaccinated group 2 (P2_AI)
dydt(2*num_grps_+1:3*num_grps_) = ...
    - Lambda .* y(2*num_grps_+1:3*num_grps_) ...
    - tauP_ * y(2*num_grps_+1:3*num_grps_) ...
    + tauP_ * y(3*num_grps_+1:4*num_grps_) ...
    - mu_ .* y(2*num_grps_+1:3*num_grps_) ...
    + [0; c_] .* (1 - vacc_) .* [0; y(2*num_grps_+1:3*num_grps_-1)] ...
    - [c_; 0] .* y(2*num_grps_+1:3*num_grps_) ...
    + [0; c_] .* vacc_ .* [0; y(num_grps_+1:2*num_grps_-1)];

%% Generate partially vaccinated group 3 (P3_AI)
dydt(3*num_grps_+1:4*num_grps_) = ...
    - Lambda .* y(3*num_grps_+1:4*num_grps_) ...
    - tauP_ * y(3*num_grps_+1:4*num_grps_) ...
    + tau_ * y(4*num_grps_+1:5*num_grps_) ...
    - mu_ .* y(3*num_grps_+1:4*num_grps_) ...
    + sigma_ * y(8*num_grps_+1:9*num_grps_) ...
    + [0; c_] .* (1 - vacc_) .* [0; y(3*num_grps_+1:4*num_grps_-1)] ...
    - [c_; 0] .* y(3*num_grps_+1:4*num_grps_) ...
    + [0; c_] .* vacc_ .* [0; y(2*num_grps_+1:3*num_grps_-1)];

%% Generate fully vaccinated group (C_AI)
dydt(4*num_grps_+1:5*num_grps_) = ...
    - tau_ * y(4*num_grps_+1:5*num_grps_) ...
    - mu_ .* y(4*num_grps_+1:5*num_grps_) ...
    + [0; c_] .* [0; y(4*num_grps_+1:5*num_grps_-1)] ...
    - [c_; 0] .* y(4*num_grps_+1:5*num_grps_) ...
    + [0; c_] .* vacc_ .* [0; y(3*num_grps_+1:4*num_grps_-1)];

%% Generate primed infectious group (I1)
dydt(5*num_grps_+1:6*num_grps_) = ...
    Lambda .* y(1:num_grps_) ...
    - gamma_ .* y(5*num_grps_+1:6*num_grps_) ...
    - mu_ .* y(5*num_grps_+1:6*num_grps_) ...
    + [0; c_] .* [0; y(5*num_grps_+1:6*num_grps_-1)] ...
    - [c_; 0] .* y(5*num_grps_+1:6*num_grps_);

%% Generate less infectious group (I2)
dydt(6*num_grps_+1:7*num_grps_) = ...
    Lambda .* y(num_grps_+1:2*num_grps_) ...
    - gamma_ .* y(6*num_grps_+1:7*num_grps_) ...
    - mu_ .* y(6*num_grps_+1:7*num_grps_) ...
    + [0; c_] .* [0; y(6*num_grps_+1:7*num_grps_-1)] ...
    - [c_; 0] .* y(6*num_grps_+1:7*num_grps_);

%% Generate mild infectious group
dydt(7*num_grps_+1:8*num_grps_) = ...
    Lambda .* y(2*num_grps_+1:3*num_grps_) ...
    - gamma_ .* y(7*num_grps_+1:8*num_grps_) ...
    - mu_ .* y(7*num_grps_+1:8*num_grps_) ...
    + [0; c_] .* [0; y(7*num_grps_+1:8*num_grps_-1)] ...
    - [c_; 0] .* y(7*num_grps_+1:8*num_grps_);

%% Generate recovered group
dydt(8*num_grps_+1:9*num_grps_) = ...
    Lambda .* y(3*num_grps_+1:4*num_grps_) ...
    + gamma_ .* (y(5*num_grps_+1:6*num_grps_) ...
    + y(6*num_grps_+1:7*num_grps_) + y(7*num_grps_+1:8*num_grps_)) ...
    - sigma_ * y(8*num_grps_+1:9*num_grps_) ...
    - mu_ .* y(8*num_grps_+1:9*num_grps_) ...
    + [0; c_] .* [0; y(8*num_grps_+1:9*num_grps_-1)] ...
    - [c_; 0] .* y(8*num_grps_+1:9*num_grps_);

dydt(9*num_grps_ + 1) = ...
    - Lambda(1) * y(9*num_grps_ + 1) - mu_(1) * y(9*num_grps_ + 1) ...
    - c_(1) * y(9*num_grps_ + 1) + B_ * f_A_;
