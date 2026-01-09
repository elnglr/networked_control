function perform_centralized_dt_analysis(F, G, h, x0, rho_target, Q_cost, R_cost)
% Merkezi Ayrık Zaman Kontrolü Analiz Fonksiyonu

disp('####################################################');
disp('# AYRIK ZAMAN (DT) MERKEZİ KONTROL ANALİZİ #');
disp('####################################################');

% Açık Çevrim Analizi
rho_open_dt = max(abs(eig(F)));
fprintf('AÇIK ÇEVRİM: Spektral Yarıçap (rho) = %.4f. Stabilite: %s\n', rho_open_dt, ternary(rho_open_dt < 1, 'Kararlı', 'KARARSIZ'));

results = struct();

% --- 1. Stabilizasyon LMI ---
[K_stab, ~] = solve_discrete_lmi_centralized(F, G);
if ~isempty(K_stab)
    F_cl = F + G * K_stab;
    results.Stabilization = struct('K', K_stab, 'Rho', max(abs(eig(F_cl))));
    fprintf('\n1. STABİLİZASYON: max(|lambda|) = %.4f. Kutsal Koşul Sağlandı.\n', results.Stabilization.Rho);
    plot_trajectory_dt(F_cl, h, x0, 'Ayrık Zaman Stabilizasyon (|lambda|<1)');
end

% --- 2. Bölge Kontrolü LMI (Hız Kısıtı) ---
[K_region, ~] = solve_discrete_lmi_region(F, G, rho_target);
if ~isempty(K_region)
    F_cl = F + G * K_region;
    results.Region = struct('K', K_region, 'Rho', max(abs(eig(F_cl))), 'Rho_Target', rho_target);
    fprintf('\n2. HIZ KONTROLÜ: max(|lambda|) = %.4f. Hedef bölge (rho < %.4f) içinde.\n', results.Region.Rho, rho_target);
    plot_trajectory_dt(F_cl, h, x0, sprintf('Ayrık Zaman Hız Kontrolü (\\rho<%.3f)', rho_target));
end

% --- 3. H2 Optimal Kontrol ---
[K_H2, trace_S] = solve_discrete_lmi_H2(F, G, Q_cost, R_cost);
if ~isempty(K_H2)
    F_cl = F + G * K_H2;
    results.H2 = struct('K', K_H2, 'Trace_S', trace_S, 'Rho', max(abs(eig(F_cl))));
    fprintf('\n3. H2 OPTİMAL KONTROL: Trace(S) = %.4f. max(|lambda|) = %.4f\n', trace_S, results.H2.Rho);
    plot_trajectory_dt(F_cl, h, x0, 'Ayrık Zaman H2 Optimal Kontrol');
end

% --- 4. H-Infinity Robust Kontrol ---
[K_Hinf, gamma_opt] = solve_discrete_lmi_Hinf(F, G, Q_cost, R_cost);
if ~isempty(K_Hinf)
    F_cl = F + G * K_Hinf;
    results.Hinf = struct('K', K_Hinf, 'Gamma', gamma_opt, 'Rho', max(abs(eig(F_cl))));
    fprintf('\n4. H-INFINITY KONTROL: Optimal Gamma = %.4f. max(|lambda|) = %.4f\n', gamma_opt, results.Hinf.Rho);
    plot_trajectory_dt(F_cl, h, x0, 'Ayrık Zaman H-Infinity Robust Kontrol');
end

fprintf('\nANALİZ SONUÇLARI ÖZETİ (DT):\n');
fprintf('----------------------------------------------------\n');
if isfield(results, 'Stabilization'), fprintf('1. Stabilizasyon: max(|lambda|) = %.4f\n', results.Stabilization.Rho); end
if isfield(results, 'Region'), fprintf('2. Hız Kontrolü: max(|lambda|) = %.4f\n', results.Region.Rho); end
if isfield(results, 'H2'), fprintf('3. H2 Optimal: Trace(S) = %.4f, max(|lambda|) = %.4f\n', results.H2.Trace_S, results.H2.Rho); end
if isfield(results, 'Hinf'), fprintf('4. H-Infinity: Gamma = %.4f, max(|lambda|) = %.4f\n', results.Hinf.Gamma, results.Hinf.Rho); end
fprintf('----------------------------------------------------\n');

end