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

incdI1 = zeros(size(sol, 1)-1, num_grps_);
for i = 1:size(sol, 1)-1    % time
    previous_lambda = contact_matrix_ ...
            * (I1(i, :)' + rho1_ * I2(i, :)' + rho2_ * I3(i, :)');
    previous = previous_lambda .* S(i, :)';
    next_lambda = contact_matrix_ ...
            * (I1(i+1, :)' + rho1_ * I2(i+1, :)' + rho2_ * I3(i+1, :)');
    next = next_lambda .* S(i+1, :)';
    incdI1(i, :) = ((previous + next)/2 * dt)';
end

incdI2 = zeros(size(sol, 1)-1, num_grps_);
for i = 1:size(sol, 1)-1    % time
    previous_lambda = contact_matrix_ ...
            * (I1(i, :)' + rho1_ * I2(i, :)' + rho2_ * I3(i, :)');
    previous = previous_lambda .* P1(i, :)';
    next_lambda = contact_matrix_ ...
            * (I1(i+1, :)' + rho1_ * I2(i+1, :)' + rho2_ * I3(i+1, :)');
    next = next_lambda .* P1(i+1, :)';
    incdI2(i, :) = ((previous + next)/2 * dt)';
end

figure('pos', [10 10 1600 900]);
plot(time_stamp(2:end), sum(incdI1 + incdI2, 2) * 1e5, 'linewidth', 2)
xlabel('time (years)')
ylabel('Incidence (cases/year per 100 000 population)')
xlim([0 10])
ylims = get(gca, 'Ylim');
ylims(1) = 0;
set(gca, 'Ylim', ylims, 'FontSize', 20)
     