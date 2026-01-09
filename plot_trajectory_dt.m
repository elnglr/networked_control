%% Discrete-Time Trajectory + Control Effort (ZOH-aware)
function plot_trajectory_dt(F_cl, G, K, h, x0, plotTitle)
    % Ayrık zamanlı kapalı çevrim yörüngeyi ve kontrol çabalarını hesaplar ve çizer.
    T_sim = 10; % 10 saniye
    n_steps = T_sim / h;
    n = length(x0);
    m = size(G,2);

    x_hist = zeros(n, n_steps+1);
    u_hist = zeros(m, n_steps+1);

    x_hist(:,1) = x0;
    u_hist(:,1) = K * x0;

    for k = 1:n_steps
        x_hist(:,k+1) = F_cl * x_hist(:,k);
        u_hist(:,k+1) = K * x_hist(:,k+1);
    end

    time = 0:h:(h*n_steps);

    % Plot states with ZOH
    plot_chemical_states(time, x_hist, [plotTitle ' - States'], 'Discrete');

    % Plot control efforts with ZOH
    figure('Name',[plotTitle ' - Control Efforts'],'Color','w');
    for i = 1:m
        subplot(m,1,i);
        stairs(time, u_hist(i,:), 'LineWidth',1.5); % ZOH effect
        grid on;
        xlabel('Time (s)');
        ylabel(['u_' num2str(i)]);
        title(['Control Input u_' num2str(i)]);
    end
    sgtitle([plotTitle ' | Discrete-Time Control Efforts']);
end