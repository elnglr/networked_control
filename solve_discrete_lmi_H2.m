function [K_H2_disc, trace_S] = solve_discrete_lmi_H2(F, G, Q_cost, R_cost)
% Ayrık Zamanlı Merkezi H2 Optimal Kontrol LMI'sını çözer.
% Amaç: İz(S)'yi minimize etmek (Çıktı enerji sınırının üst sınırı).

n = size(F, 1); % 12
m = size(G, 2); % 6

% Genelleştirilmiş Sistem Matrisleri (H2 için yaygın varsayımlar)
Bw = eye(n); % Gürültü girdi matrisi (w, tüm durumlara etki eden gürültü)
Cz = [sqrtm(Q_cost); zeros(m, n)]; % Çıktı durum ağırlığı (z = Cz*x + Du*u)
Du = [zeros(n, m); sqrtm(R_cost)]; % Çıktı girdi ağırlığı

p_cost = size(Cz, 1); % Maliyet çıktısı boyutu (n + m)

disp('--- Ayrık Zaman Merkezi H2 Optimal Kontrol Çözülüyor ---');

%% 1. YALMIP Değişkenlerini Tanımlama
Q = sdpvar(n, n, 'symmetric');    % Lyapunov Matrisi (Q=P)
L = sdpvar(m, n, 'full');         % Yardımcı Değişken (K*Q)
S = sdpvar(p_cost, p_cost, 'symmetric'); % Çıktı Enerjisi Üst Sınırı

%% 2. LMI Kısıtlamaları

epsilon = 1e-6; 

% LMI 1: Lyapunov Kararlılık/Kontrol Kısıtı (Schur Komplementi)
% [Q - F*Q*F' - F*L'*G' - G*L*F' - Bw*Bw', G*L; L'*G', Q] > 0
LMI_1 = [[Q - F*Q*F' - F*L'*G' - G*L*F' - Bw*Bw',  G*L;
          L'*G',                                   Q]] >= epsilon*eye(2*n);

% LMI 2: H2 Norm Kısıtı (Schur Komplementi)
% [S, (Cz*Q + Du*L); (Cz*Q + Du*L)', Q] > 0
LMI_2 = [ S,                 Cz*Q + Du*L;
          (Cz*Q + Du*L)',    Q ] >= epsilon*eye(p_cost + n);

constraints = [LMI_1, LMI_2, Q >= epsilon*eye(n)];

%% 3. Optimizasyon (Amaç: H2 Normunu Minimize Etmek => Trace(S)'yi Minimize Etmek)
options = sdpsettings('verbose', 0, 'solver', 'sdpt3'); 
sol = optimize(constraints, trace(S), options);

%% 4. Sonuçlar
if sol.problem == 0
    disp('Başarılı: Ayrık H2 LMI Çözüldü.');
    
    Q_val = double(Q);
    L_val = double(L);
    
    K_H2_disc = L_val / Q_val; 
    trace_S = double(trace(S));
    
    fprintf('Ayrık H2 Normu Üst Sınırı (Trace(S)): %.4f\n', trace_S);
else
    disp('HATA: Ayrık H2 LMI Çözümü Başarısız Oldu.');
    disp(yalmiperror(sol.problem));
    K_H2_disc = [];
    trace_S = NaN;
end
end