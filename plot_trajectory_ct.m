function plot_trajectory_ct(A_cl, B, K, h, x0, plotTitle)
    T_sim = 15; 
    time = 0:h:T_sim;
    n = length(x0);
    m = size(B,2);
    x_hist = zeros(n, length(time));
    u_hist = zeros(m, length(time));
    
    for k = 1:length(time)
        x_hist(:, k) = expm(A_cl * time(k)) * x0;
        u_hist(:, k) = K * x_hist(:, k);  
    end
    
    % 1. States (H1, H2, H3, XB3) - Bu fonksiyonun içinde figür açıldığı varsayılıyor
    plot_chemical_states(time, x_hist, plotTitle, 'Continuous-Time');
    saveas(gcf, sprintf('States_%s.jpg', plotTitle)); % En son aktif figürü kaydeder

    % 2. Figure: Feed Flow Rates
    flow_indices = [1, 3, 5];
    flow_titles = {'F_{f1} (Feed Flow 1)', 'F_{f2} (Feed Flow 2)', 'F_R (Recycle Flow)'};
    f1 = figure('Name', [plotTitle ' - Flows'], 'Color', 'w', 'Position', [100, 100, 600, 800]);
    for i = 1:3
        subplot(3, 1, i);
        plot(time, u_hist(flow_indices(i),:), 'LineWidth', 1.5, 'Color', [0.1 0.4 0.1]);
        grid on; xlabel('Time (s)'); ylabel('\Delta Flow [kg/s]');
        title(flow_titles{i}); hold on; line([time(1) time(end)], [0 0], 'Color', 'k', 'LineStyle', ':');
    end
    sgtitle([plotTitle ' | Flow Deviations']);
    saveas(f1, sprintf('Flows_%s.jpg', plotTitle));

    % 3. Figure: Heat Duties
    heat_indices = [2, 4, 6];
    heat_titles = {'Q_1 (Heat 1)', 'Q_2 (Heat 2)', 'Q_3 (Heat 3)'};
    f2 = figure('Name', [plotTitle ' - Heat'], 'Color', 'w', 'Position', [750, 100, 600, 800]);
    for i = 1:3
        subplot(3, 1, i);
        plot(time, u_hist(heat_indices(i),:), 'LineWidth', 1.5, 'Color', [0.8 0.2 0.2]);
        grid on; xlabel('Time (s)'); ylabel('\Delta Heat [kJ/s]');
        title(heat_titles{i}); hold on; line([time(1) time(end)], [0 0], 'Color', 'k', 'LineStyle', ':');
    end
    sgtitle([plotTitle ' | Heat Deviations']);
    saveas(f2, sprintf('Heat_%s.jpg', plotTitle));
end