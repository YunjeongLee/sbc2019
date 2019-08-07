function [] = visualize(sol, params, time_stamp)
%% Assign Parameter
names = params(:,1);
for k = 1:size(params,1)
    cmd_string = sprintf('%s_ = params{%d, 2};', names{k}, k);
    eval(cmd_string);
end

%% Divide into compartments
S = sol(:, 1:num_grps_);
P1 = sol(:, num_grps_+1:2*num_grps_);
P2 = sol(:, 2*num_grps_+1:3*num_grps_);
P3 = sol(:, 3*num_grps_+1:4*num_grps_);
C = sol(:, 4*num_grps_+1:5*num_grps_);
I1 = sol(:, 5*num_grps_+1:6*num_grps_);
I2 = sol(:, 6*num_grps_+1:7*num_grps_);
I3 = sol(:, 7*num_grps_+1:8*num_grps_);
R = sol(:, 8*num_grps_+1:9*num_grps_);

dt = time_stamp(2) - time_stamp(1);

previous_lambda = contact_matrix_ ...
    * (I1(1:end-1, :)' + rho1_ * I2(1:end-1, :)' + rho2_ * I3(1:end-1, :)');
previous = previous_lambda .* S(1:end-1, :)';
next_lambda = contact_matrix_ ...
    * (I1(2:end, :)' + rho1_ * I2(2:end, :)' + rho2_ * I3(2:end, :)');
next = next_lambda .* S(2:end, :)';
incdI1 = ((previous + next)/2 * dt)';

% incdI1 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* S';

previous_lambda = contact_matrix_ ...
    * (I1(1:end-1, :)' + rho1_ * I2(1:end-1, :)' + rho2_ * I3(1:end-1, :)');
previous = previous_lambda .* P1(1:end-1, :)';
next_lambda = contact_matrix_ ...
    * (I1(2:end, :)' + rho1_ * I2(2:end, :)' + rho2_ * I3(2:end, :)');
next = next_lambda .* P1(2:end, :)';
incdI2 = ((previous + next)/2 * dt)';
% incdI2 = contact_matrix_ * (I1' + rho1_ * I2' + rho2_ * I3') .* P1';

incd = incdI1 + incdI2;

figure('pos', [10 10 1600 900]);
plot(time_stamp(2:end)/365, sum(incdI1 + incdI2, 2) * 1e6, 'linewidth', 2)
xlabel('time (years)')
ylabel('Incidence (cases/year per 1,000,000 population)')
ylim([0 500])

figure('pos', [10 10 1600 900]);
plot(time_stamp(2:end)/365, (incdI1 + incdI2) * 1e6, 'linewidth', 2)
xlabel('time (years)')
ylabel('Incidence (cases/year per 1,000,000 population)')
ylim([0 100])
     