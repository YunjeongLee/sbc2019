function [] = visualize(incd, time_stamp, group_for_label, text_for_legend, ylims)
%% Plot the graphs
figure('pos', [10 10 1600 900]);
plot(time_stamp(2:end)/365, incd * 1e6, 'linewidth', 2)
xlabel('time (years)')
ylabel(sprintf( ...
    'Incidence of %s month age group (cases/year per 1,000,000 population)', ...
    group_for_label))
legend(text_for_legend)
ylim(ylims)
