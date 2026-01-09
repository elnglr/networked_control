function perform_centralized_ct_analysis(A, B, h, x0, alpha, zeta, Q_cost, R_cost)

disp('####################################################');
disp('# Continuous-Time (CT) Centralized Control #');
disp('####################################################');

% Open Loop Anaylsis
rho_open_ct = max(real(eig(A)));
fprintf('Open Loop : Spektral Apsis (rho) = %.4f. Stability: %s\n', rho_open_ct, ternary(rho_open_ct < 0, 'Stable', 'Unstable'));

results = struct();

% --- 1. Stabilization LMI ---
[K_stab] = solve_continuous_lmi_centralized(A, B); % Sizin ilk yazdığınız fonksiyon
if ~isempty(K_stab)
    A_cl = A + B * K_stab;
    results.Stabilization = struct('K', K_stab, 'Rho', max(real(eig(A_cl))));
    fprintf('\n1. Stabilization: max(Re(lambda)) = %.4f. Kutsal Koşul Sağlandı.\n', results.Stabilization.Rho);
    plot_trajectory_ct(A_cl, h, x0, 'Sürekli Zaman Stabilizasyon (Re(lambda)<0)');
end


% --- 2. Bölge Kontrolü LMI (Hız ve Sönümleme) ---
[K_region] = solve_continuous_lmi_region(A, B, alpha, zeta); % Sizin bölge kontrolü fonksiyonunuz
if ~isempty(K_region)
    A_cl = A + B * K_region;
    results.Region = struct('K', K_region, 'Rho', max(real(eig(A_cl))), 'Alpha_Target', alpha, 'Zeta_Target', zeta);
    fprintf('\n2. BÖLGE KONTROLÜ: max(Re(lambda)) = %.4f. Hedef bölge (-%.4f, zeta>%.4f) içinde.\n', results.Region.Rho, alpha, zeta);
    plot_trajectory_ct(A_cl, h, x0, sprintf('Sürekli Zaman Bölge Kontrolü (\\alpha=%.3f, \\zeta=%.3f)', alpha, zeta));
end

% --- 3. H2 Optimal Kontrol ---
[K_H2, trace_S] = solve_continuous_lmi_H2(A, B, Q_cost, R_cost);
if ~isempty(K_H2)
    A_cl = A + B * K_H2;
    results.H2 = struct('K', K_H2, 'Trace_S', trace_S, 'Rho', max(real(eig(A_cl))));
    fprintf('\n3. H2 OPTİMAL KONTROL: Trace(S) = %.4f. max(Re(lambda)) = %.4f\n', trace_S, results.H2.Rho);
    plot_trajectory_ct(A_cl, h, x0, 'Sürekli Zaman H2 Optimal Kontrol');
end

% --- 4. H-Infinity Robust Kontrol ---
[K_Hinf, gamma_opt] = solve_continuous_lmi_Hinf(A, B, Q_cost, R_cost);
if ~isempty(K_Hinf)
    A_cl = A + B * K_Hinf;
    results.Hinf = struct('K', K_Hinf, 'Gamma', gamma_opt, 'Rho', max(real(eig(A_cl))));
    fprintf('\n4. H-INFINITY KONTROL: Optimal Gamma = %.4f. max(Re(lambda)) = %.4f\n', gamma_opt, results.Hinf.Rho);
    plot_trajectory_ct(A_cl, h, x0, 'Sürekli Zaman H-Infinity Robust Kontrol');
end

fprintf('\nANALİZ SONUÇLARI ÖZETİ (CT):\n');
fprintf('----------------------------------------------------\n');
if isfield(results, 'Stabilization'), fprintf('1. Stabilizasyon: max(Re(lambda)) = %.4f\n', results.Stabilization.Rho); end
if isfield(results, 'Region'), fprintf('2. Bölge Kontrolü: max(Re(lambda)) = %.4f\n', results.Region.Rho); end
if isfield(results, 'H2'), fprintf('3. H2 Optimal: Trace(S) = %.4f, max(Re(lambda)) = %.4f\n', results.H2.Trace_S, results.H2.Rho); end
if isfield(results, 'Hinf'), fprintf('4. H-Infinity: Gamma = %.4f, max(Re(lambda)) = %.4f\n', results.Hinf.Gamma, results.Hinf.Rho); end
fprintf('----------------------------------------------------\n');

end