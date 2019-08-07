function [] = visualize(incd, time_stamp, group_for_ttl, text_for_legend, xlims, ylims)
%% Plot the graphs
figure('pos', [10 10 1600 900]);
plot(time_stamp(2:end)/365, incd * 1e6, 'linewidth', 2)
grid on; grid minor;
xlabel('time (years)')
ylabel(sprintf( ...
    'Incidence of %s age group (cases/year per 1,000,000 population)', ...
    group_for_label))
legend(text_for_legend)
xlim(xlims)
ylim(ylims)
set(gca, 'fontsize', 20);
