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
          'p1', 0.5; ...
          'p2', 0.6; ...
          'p3', 0.7; ...
          'vacc', vacc_val; ...
          'num_grps', num_grps_val};

%% Change 6 years old and 11 years old vaccination proportion
time_stamp = 0:50*365;
year6_vaccine = [0.1, 0.7, 0.7];
year11_vaccine = [0.1, 0.1, 0.7];
incd_aggregate_baby = zeros(length(time_stamp), length(year6_vaccine));
incd_aggregate_all = zeros(length(time_stamp), length(year6_vaccine));
for i = 1:length(year6_vaccine)
    % Copy vaccine proportion
    p_val_temp = p_val;
    
    % Update vaccine proportion
    p_val_temp(end-1) = year6_vaccine(i);
    p_val_temp(end) = year11_vaccine(i);
    
    % Vaccine update
    vacc_val = zeros(num_grps_val, 1);
    for k = 1:num_grps_val
        f = p_val_temp * VE;
        for j = 1:length(vaccine_age)
            vacc_val(k) = vacc_val(k) + f(j) * kDelta(age_grps(k), vaccine_age(j));
        end
    end
    
    % Update params
    params_temp = params;
    params_temp{end-1,2} = vacc_val;

    % Solve ODE
    fode = @(t, y) model_pertussis(t, y, params_temp);
    options = odeset('NonNegative', 1:num_grps_val*9);
    [~, sol] = ode45(fode, time_stamp, y0, options);
    
    incd = get_incidence(sol, params_temp, time_stamp);
    incd_aggregate_baby(:,i) = sum(incd(:, 1), 2);
    incd_aggregate_all(:,i) = sum(incd, 2);
end

%% Visualization
mkdir('results/change_vacc_prop');
group_for_title_baby = '0-1 year';
group_for_title_all = '20+ year';
text_for_legend = {'6 years: 10%, 11 years: 10%', '6 years: 70%, 11 years: 10%', ...
    '6 years: 70%, 11 years: 70%'};
xlims = [time_stamp(1), time_stamp(end)]/365;
ylims = [0, 2];
visualize(incd_aggregate_baby, time_stamp, group_for_title_baby, text_for_legend, xlims, ylims)
saveas(gca, 'results/change_vacc_prop/incd_below1year.png', 'png');
ylims = [0, 300];
visualize(incd_aggregate_all, time_stamp, group_for_title_all, text_for_legend, xlims, ylims)
saveas(gca, 'results/change_vacc_prop/incd_all.png', 'png');

%% Draw reduced proportion as bar graph
total_incd_baby = sum(incd_aggregate_baby);
total_incd_all = sum(incd_aggregate_all);

reduced_ratio_baby = 1 - total_incd_baby./max(total_incd_baby);
reduced_ratio_all = 1 - total_incd_all./max(total_incd_all);

figure('pos', [10 10 1600 700]);
subplot(121)
bar(reduced_ratio_baby(2:end)*100);
ylabel('reduced level (%)')
yticks(1:10);
x_labels = {sprintf('6 years: 10%% \n 11 years: 10%%'), ...
    sprintf('6 years: 70%% \n 11 years: 10%%'), ...
    sprintf('6 years: 70%% \n 11 years: 70%%')};
title('0-2 months babies')
set(gca, 'fontsize', 15);
set(gca,'XTickLabel', {'','',''});
format_ticks(gca, x_labels);
ylim([0 10])

subplot(122)
bar(reduced_ratio_all(2:end)*100);
ylabel('reduced level (%)')
yticks(1:10);
x_labels = {sprintf('6 years: 10%% \n 11 years: 10%%'), ...
    sprintf('6 years: 70%% \n 11 years: 10%%'), ...
    sprintf('6 years: 70%% \n 11 years: 70%%')};
title('All age groups')
set(gca, 'fontsize', 15);
set(gca,'XTickLabel', {'','',''});
format_ticks(gca, x_labels);
ylim([0 10])
