function plot_eigenvalues(K, sys_type, structure_name, A_or_F, B_or_G, target_val, zeta)
    % K: Kazanç matrisi
    % sys_type: 'CT' (Continuous) veya 'DT' (Discrete)
    % structure_name: Yapı ismi (Centralized, String, vb.)
    % target_val: CT için Alpha, DT için Radius
    
    if nargin < 7, zeta = 0; end 

    % Kapalı döngü matrisi ve özdeğer hesabı
    CL_mat = A_or_F + B_or_G * K;
    evs = eig(CL_mat);
    re_evs = real(evs);
    im_evs = imag(evs);
    
    % --- Dinamik Metrik Hesaplamaları ---
    if strcmpi(sys_type, 'CT')
        spectral_metric = max(re_evs); % Spectral Abscissa (En sağdaki kutup)
        metric_label = 'Spectral Abscissa (\alpha)';
    else
        spectral_metric = max(abs(evs)); % Spectral Radius (En dıştaki kutup)
        metric_label = 'Spectral Radius (\rho)';
    end

    % --- Figure Hazırlığı ---
    figure('Color', 'w', 'Position', [200, 200, 700, 600]);
    hold on; grid on; box on; axis equal;
    
    if strcmpi(sys_type, 'CT')
        % --- S-Plane (Continuous Time) ---
        % Eksen limitleri (Tüm kutupları kapsayacak geniş bir alan)
        x_min = min([-target_val*1.5, min(re_evs)*1.2, -5]);
        x_max = max([2, max(re_evs)*1.2]);
        y_max = max([abs(im_evs)*1.5; 5]);
        
        % Bölge Sınırları
        line([-target_val -target_val], [-y_max y_max], 'Color', [0.8 0 0], ...
            'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Alpha Limit');
        
        if zeta > 0 && zeta < 1
            theta_rad = acos(zeta);
            x_line = linspace(x_min, 0, 100);
            plot(x_line, -x_line * tan(theta_rad), 'k:', 'LineWidth', 1.5, 'DisplayName', 'Damping Sector');
            plot(x_line, x_line * tan(theta_rad), 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        
        % Ana eksenler
        line([0 0], [-y_max y_max], 'Color', 'k', 'LineWidth', 1);
        line([x_min x_max], [0 0], 'Color', 'k', 'LineWidth', 1);
        
        xlim([x_min x_max]); ylim([-y_max y_max]);
        xlabel('Real Part'); ylabel('Imaginary Part');
        
    else
        % --- Z-Plane (Discrete Time) ---
        % Limitler: Birim çember ve kutupları kapsayacak şekilde
        view_limit = max([1.2, max(abs(evs))*1.2, target_val*1.2]);
        
        % Birim Çember ve Hedef Çember
        t = linspace(0, 2*pi, 200);
        plot(cos(t), sin(t), 'k', 'LineWidth', 1, 'DisplayName', 'Unit Circle');
        plot(target_val*cos(t), target_val*sin(t), 'r--', 'LineWidth', 2, 'DisplayName', 'Target Radius');
        
        % Eksenler
        line([-view_limit view_limit], [0 0], 'Color', [0.5 0.5 0.5]);
        line([0 0], [-view_limit view_limit], 'Color', [0.5 0.5 0.5]);
        
        xlim([-view_limit view_limit]); ylim([-view_limit view_limit]);
        xlabel('Real (z)'); ylabel('Imaginary (z)');
    end
    
    % --- Özdeğerleri Çiz (Görsel olarak belirgin Blue X) ---
    plot(re_evs, im_evs, 'bx', 'MarkerSize', 12, 'LineWidth', 2.5, 'DisplayName', 'CL Eigenvalues');
    
    % --- Başlık (İstediğin Eklemelerle) ---
    title_str = sprintf('%s Structure\n%s = %.4f', structure_name, metric_label, spectral_metric);
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    
    legend('Location', 'northeastoutside');
    hold off;
end