clear; close all; clc;

%% Parameter values
N = 1;
age_grps = [0, 2/12, 4/12, 6/12, 1, 18/12, 2:14, 15:5:50, 55, 65]*365;
num_grps_val = length(age_grps);
pop_age = diff([age_grps, 75*365]')/75/365*N;

% Initial condition
S0_val = 0.999 * pop_age;
P10_val = zeros(num_grps_val, 1);
P20_val = zeros(num_grps_val, 1);
P30_val = zeros(num_grps_val, 1);
C0_val = zeros(num_grps_val, 1);
I10_val = pop_age - S0_val;
I20_val = zeros(num_grps_val, 1);
I30_val = zeros(num_grps_val, 1);
R0_val = zeros(num_grps_val, 1);
X0_val = 0;

% Parameters
sigma_val = 1/11/365;
tau_val = 1/2/365;
tauP_val = 1/2/365;
sigma0_val = 1/100/365;
mu_val = [zeros(num_grps_val-1, 1); 1/10/365];
rho1_val = 0.5;
rho2_val = 0.25;
gamma_val = 1/21;
c_val = 1./(diff(age_grps))';
B_val = 1/75/365;

% Parameters for maternal immunity
f_A_val = nan;

%% Parameters w.r.t. vaccine
VE = 1;
p_val = [0.95 * ones(3, 1); 0.85; 0.1; 0.1];
vaccine_age = [2/12, 4/12, 6/12, 18/12, 6, 11]*365;

%% Generate vaccine vector
vacc_val = zeros(num_grps_val, 1);
for i = 1:num_grps_val
    f = p_val * VE;
    for j = 1:length(vaccine_age)
        vacc_val(i) = vacc_val(i) + f(j) * kDelta(age_grps(i), vaccine_age(j));
    end
end

%% Load contact matrix
b = [0.02, 0.12, 0.23, 0.73, 2.47, 0.95, 0.54, 0.16, 0.00];
contact_matrix_val = generate_contact_matrix(b, num_grps_val);

%% Generate params
y0 = [S0_val; P10_val; P20_val; P30_val; C0_val; I10_val; I20_val; I30_val; R0_val; X0_val];
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
          'p1', 0.5; ...
          'p2', 0.6; ...
          'p3', 0.7; ...
          'f_A', f_A_val; ...
          'vacc', vacc_val; ...
          'num_grps', num_grps_val};

%% Change maternal immunity values
time_stamp = 0:50*365;
tauA = 2.5*365; % unit: days
frac_immune_mother_range = [0, 0.65, 0.8, 0.99];
incd_aggregate_baby = zeros(length(time_stamp), length(frac_immune_mother_range));
incd_aggregate_all = zeros(length(time_stamp), length(frac_immune_mother_range));
for i = 1:length(frac_immune_mother_range)
    % Assign fractions of immuned mothers
    frac_immune_mother = frac_immune_mother_range(i);
    f_A_val = sigma_val * tauA * frac_immune_mother;
    
    % Parameters w.r.t. vaccine
    p_val = [0.95 * ones(3, 1); 0.85; 0.1; 0.1; frac_immune_mother; frac_immune_mother; frac_immune_mother];
    vaccine_age = [2/12, 4/12, 6/12, 18/12, 6, 11, 20, 25, 30]*365;

    % Generate vaccine vector
    vacc_val = zeros(num_grps_val, 1);
    for k = 1:num_grps_val
        f = p_val * VE;
        for j = 1:length(vaccine_age)
            vacc_val(k) = vacc_val(k) + f(j) * kDelta(age_grps(k), vaccine_age(j));
        end
    end
    
    % Update params
    params_temp = params;
    params_temp{end-2, 2} = f_A_val;
    params_temp{end-1, 2} = vacc_val;
    
    % Solve ode
    fode = @(t, y) model_pertussis_maternal(t, y, params_temp);
    options = odeset('NonNegative', 1:num_grps_val*9+1);
    [~, sol] = ode45(fode, time_stamp, y0, options);

    % Get incidence
    [incd, ~] = get_incidence(sol, params_temp);
    incd_aggregate_baby(:,i) = sum(incd(:, 1), 2);
    incd_aggregate_all(:,i) = sum(incd, 2);
end

%% Visualize
mkdir('results/maternal_immune')
ttl_baby = '0-2 months';
ttl_all = 'All age groups';
lgd = {};
for i = 1:length(frac_immune_mother_range)
    lgd{end+1} = sprintf('prop: %d%%', frac_immune_mother_range(i) * 100);
end
% xlims = [10*365+1, time_stamp(end)]/365;
xlims = [time_stamp(1), time_stamp(end)]/365;
ylims = [0, 1];
visualize(incd_aggregate_baby, time_stamp, ttl_baby, lgd, xlims, ylims)
saveas(gca, 'results/maternal_immune/incd_below2months.png', 'png');
ylims = [0, 50];
visualize(incd_aggregate_all, time_stamp, ttl_all, lgd, xlims, ylims)
saveas(gca, 'results/maternal_immune/incd_all', 'png');

%% Visualize as bar graph
wts_period = 10*365+1:time_stamp(end);
total_incd_baby = sum(incd_aggregate_baby(wts_period, :));
total_incd_all = sum(incd_aggregate_all(wts_period, :));

reduced_ratio_baby = 1 - total_incd_baby./(total_incd_baby(1));
reduced_ratio_all = 1 - total_incd_all./(total_incd_all(1));

figure('pos', [10 10 1600 700]);
subplot(1,2,1)
bar(reduced_ratio_baby(2:end)*100);
ylabel('reduced level (%)')
yticks(0:5:50);
x_labels = lgd(2:end);
title('0-2 months babies')
set(gca, 'fontsize', 18);
set(gca,'XTickLabel', {'','',''});
format_ticks(gca, x_labels);
ylim([0 50])

subplot(1,2,2)
bar(reduced_ratio_all(2:end)*100);
ylabel('reduced level (%)')
yticks(0:5:50);
x_labels = lgd(2:end);
title('All age groups')
set(gca, 'fontsize', 18);
set(gca,'XTickLabel', {'','',''});
format_ticks(gca, x_labels);
ylim([0 50])

saveas(gca, 'results/maternal_immune/effect_change_maternal_prop.png', 'png');
