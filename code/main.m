clear; close all; clc;

%% Parameter values
N = 40091359;
S0_val = 0.99;
I0_val = 1 - S0_val;
R0_val = 0;

sigma_val = 1/11;
tau_val = 1/2;
tauP_val = 1/2;
sigma0_val = 1/100;
VE = 0.9;
% p_val = [0; 0.95 * ones(3, 1); 0; 0.85; zeros(4, 1); 0.95; zeros(19, 1)];
p_val = [0.95 * ones(3, 1); 0.85; 0.95];
mu_val = 1/75.5;
rho1_val = 0.5;
rho2_val = 0.25;
gamma = 365/21;

%% Load contact matrix


%% Generate params

