function [K_Hinf_disc, gamma_opt] = solve_discrete_lmi_Hinf(F, G, Q_cost, R_cost)
% Ayrık Zamanlı Merkezi H-Infinity Robust Kontrol LMI'sını çözer.
% Amaç: En kötü durum kazancı (gamma)'yı minimize etmek.

n = size(F, 1); % 12
m = size(G, 2); % 6

% Genelleştirilmiş Sistem Matrisleri (Hinf için yaygın varsayımlar)
Bw = eye(n); % Gürültü/Disturbance girdi matrisi
Cz = [sqrtm(Q_cost); zeros(m, n)]; % Çıktı durum ağırlığı
Du = [zeros(n, m); sqrtm(R_cost)]; % Çıktı girdi ağırlığı
Dw = zeros(size(Cz, 1), n); % Disturbance'den çıktıya doğrudan besleme yok
p_cost = size(Cz, 1); % Maliyet çıktısı boyutu (n + m)

disp('--- Ayrık Zaman Merkezi H-Infinity Robust Kontrol Çözülüyor ---');

%% 1. YALMIP Değişkenlerini Tanımlama
Q = sdpvar(n, n, 'symmetric');    % Lyapunov Matrisi (Q=P)
L = sdpvar(m, n, 'full');         % Yardımcı Değişken (K*Q)
gamma = sdpvar(1);                % H-Infinity Normu (Minimize edilecek)

%% 2. LMI Kısıtlamaları (Bounded Real Lemma - Q formunda)

epsilon = 1e-6; 

% Matris Blokları (Ayrık Zaman için Q=P kullanılarak lineerleştirilmiş BRL)
% [F*Q*F' + G*L*F' + F*L'*G' + G*L*G' - Q, F*Q*Cz' + G*L*Cz' + Bw*Dw'; 
%  ... simetrik ..., Q - Bw*Bw', Q*Cz' + L'*Du' + Bw*Dw';
%  ... simetrik ..., Q*Cz' + L'*Du' + Bw*Dw', Q - Du*Du' - Dw*Dw'] > 0
% Bu, 4x4 matris formunun Schur Komplementi alınmış halidir,
% ancak daha basit 3x3 Schur Komplementi formunu kullanacağız:

% Matris Blokları (Bu form 3x3 bloklu BRL matrisini kullanır)
M11 = -Q;
M12 = Q*F' + L'*G';
M13 = Q*Cz' + L'*Du';
M14 = Bw;

M22 = -Q; % F'in katsayısı
M23 = Q*Cz' + L'*Du';
M24 = Bw;

M33 = -gamma*eye(p_cost); % Çıktı (z)
M34 = Dw;

M44 = -gamma*eye(n); % Disturbance (w)

% Ana Bounded Real Lemma Matrisi (Simetrik 3x3 Blok Matris, w ve z kısıtlı)
% Daha yaygın Schur Komplementi formatını kullanalım:
LMI_Hinf_schur = [ 
    Q,             F*Q + G*L,       Bw,         zeros(n, p_cost);
    (F*Q + G*L)',  Q,             zeros(n, n),   Cz'*Q + Du'*L;
    Bw',           zeros(n, n),   gamma*eye(n), Dw';
    zeros(p_cost, n), Cz*Q + Du*L, Dw,         gamma*eye(p_cost)
];

constraints = [LMI_Hinf_schur >= epsilon*eye(3*n + p_cost), Q >= epsilon*eye(n), gamma >= epsilon];

%% 3. Optimizasyon (Amaç: H-Infinity Normu gamma'yı Minimize Etmek)
options = sdpsettings('verbose', 0, 'solver', 'sdpt3'); 
sol = optimize(constraints, gamma, options);

%% 4. Sonuçlar
if sol.problem == 0
    disp('Başarılı: Ayrık H-Infinity LMI Çözüldü.');
    
    Q_val = double(Q);
    L_val = double(L);
    
    K_Hinf_disc = L_val / Q_val; 
    gamma_opt = double(gamma);
    
    fprintf('Optimal Ayrık H-Infinity Normu (gamma): %.4f\n', gamma_opt);
else
    disp('HATA: Ayrık H-Infinity LMI Çözümü Başarısız Oldu.');
    disp(yalmiperror(sol.problem));
    K_Hinf_disc = [];
    gamma_opt = NaN;
end
end