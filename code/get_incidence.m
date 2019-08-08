function [incdI1, incdI2] = get_incidence(sol, params)
%% Assign Parameter
names = params(:,1);
for k = 1:size(params,1)
    cmd_string = sprintf('%s_ = params{%d, 2};', names{k}, k);
    eval(cmd_string);
end

%% Divide into compartments
if ~exist('f_A_', 'var')
    S = sol(:, 1:num_grps_);
    P1 = sol(:, num_grps_+1:2*num_grps_);
    I1 = sol(:, 5*num_grps_+1:6*num_grps_);
    I2 = sol(:, 6*num_grps_+1:7*num_grps_);
    I3 = sol(:, 7*num_grps_+1:8*num_grps_);
    
    %% Incidence for I1
    incdI1 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* S';
    
    %% Incidence for I2
    incdI2 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* P1';
    
    %% Merge two incidences
    incd = incdI1' + incdI2';
    
else
    S = sol(:, 1:num_grps_);
    P1 = sol(:, num_grps_+1:2*num_grps_);
    I1 = sol(:, 5*num_grps_+1:6*num_grps_);
    I2 = sol(:, 6*num_grps_+1:7*num_grps_);
    I3 = sol(:, 7*num_grps_+1:8*num_grps_);
    X0 = sol(:, 9*num_grps_+1);
    
    %% Incidence for I1
    incdI1 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* (S + X0)';
    
    %% Incidence for I2
    incdI2 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* P1';
    
    %% Merge two incidences
    incd = incdI1' + incdI2';
    
end
end