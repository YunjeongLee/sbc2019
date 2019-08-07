clear; close all; clc;

%% Parameter values
N = 40091359;
age_grps = [0, 2/12, 4/12, 6/12, 1, 18/12, 2:14, 15:5:50, 55, 65];
vaccine_age = [2/12, 4/12, 6/12, 18/12, 6];
num_grps_val = length(age_grps);

% Initial condition
S0_val = ones(num_grps_val, 1) - 1/N;
P10_val = zeros(num_grps_val, 1);
P20_val = zeros(num_grps_val, 1);
P30_val = zeros(num_grps_val, 1);
C0_val = zeros(num_grps_val, 1);
I10_val = 1 - S0_val;
I20_val = zeros(num_grps_val, 1);
I30_val = zeros(num_grps_val, 1);
R0_val = zeros(num_grps_val, 1);

% Parameters
sigma_val = 1/11;
tau_val = 1/2;
tauP_val = 1/2;
sigma0_val = 1/100;
VE = 0.9;
p_val = [0.95 * ones(3, 1); 0.85; 0.95];
mu_val = [zeros(num_grps_val-1, 1); 1];
rho1_val = 0.5;
rho2_val = 0.25;
gamma_val = 365/21;
c_val = 1./(diff(age_grps))';
B_val = 1/75;

%% Generate vaccine vector
vacc_val = zeros(num_grps_val, 1);
for i = 1:num_grps_val
    f = p_val * VE;
    vacc_val(i) = sum(f * ismember(age_grps(i), vaccine_age));
end

%% Load contact matrix
b = [0.02, 0.12, 0.23, 0.73, 2.47, 0.95, 0.54, 0.16, 0.00];
contact_matrix_val = generate_contact_matrix(b, num_grps_val);

%% Generate params
y0 = [S0_val; P10_val; P20_val; P30_val; C0_val; I10_val; I20_val; I30_val; R0_val];
params = {'contact_matrix', contact_matrix_val; ...
          'sigma', sigma_val; ...
          'sigma0', sigma0_val; ...
          'tau', tau_val; ...
          'tauP', tauP_val; ...
          'mu', mu_val; ...
          'rho1', rho1_val; ...
          'rho2', rho2_val; ...
          'gamma', gamma_val; ...
          'c', c_val; ...
          'B', B_val; ...
          'vacc', vacc_val; ...
          'num_grps', num_grps_val};

%% Solve ODE
fode = @(t, y) model_pertussis(t, y, params);
time_stamp = linspace(0, 20, 1e3);
[~, sol] = ode15s(fode, time_stamp, y0);

