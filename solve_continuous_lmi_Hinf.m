function [K_Hinf, gamma_opt] = solve_continuous_lmi_Hinf(A, B, Q_cost, R_cost)
% Sürekli Zamanlı Merkezi H-Infinity Robust Kontrol LMI'sını çözer.
% Amaç: En kötü durum kazancı (gamma)'yı minimize etmek.

n = size(A, 1); % 12
m = size(B, 2); % 6

% Genelleştirilmiş Sistem Matrisleri (Hinf için yaygın varsayımlar)
Bw = eye(n); % Gürültü/Disturbance girdi matrisi
Cz = [sqrtm(Q_cost); zeros(m, n)]; % Çıktı durum ağırlığı
Du = [zeros(n, m); sqrtm(R_cost)]; % Çıktı girdi ağırlığı
Dw = zeros(size(Cz, 1), n); % Disturbance'den çıktıya doğrudan besleme yok
p_cost = size(Cz, 1); % Maliyet çıktısı boyutu (n + m)

disp('--- Sürekli Zaman Merkezi H-Infinity Robust Kontrol Çözülüyor ---');

%% 1. YALMIP Değişkenlerini Tanımlama
Q = sdpvar(n, n, 'symmetric');    % Lyapunov Matrisi (P^-1)
L = sdpvar(m, n, 'full');         % Yardımcı Değişken (K*Q)
gamma = sdpvar(1);                % H-Infinity Normu (Minimize edilecek)

%% 2. LMI Kısıtlamaları (Bounded Real Lemma - Q formunda)

epsilon = 1e-6; 

% Matris Blokları
M11 = A*Q + Q*A' + B*L + L'*B';
M12 = Bw;
M13 = Q*Cz' + L'*Du';
M22 = -gamma*eye(n);
M23 = Dw';
M33 = -gamma*eye(p_cost);

% Ana Bounded Real Lemma Matrisi (Simetrik 3x3 Blok Matris)
LMI_Hinf = [ M11, M12, M13;
             M12', M22, M23;
             M13', M23', M33 ];
             
constraints = [LMI_Hinf <= -epsilon*eye(n + n + p_cost), Q >= epsilon*eye(n), gamma >= epsilon];

%% 3. Optimizasyon (Amaç: H-Infinity Normu gamma'yı Minimize Etmek)
options = sdpsettings('verbose', 0, 'solver', 'sdpt3'); 
sol = optimize(constraints, gamma, options);

%% 4. Sonuçlar
if sol.problem == 0
    disp('Başarılı: H-Infinity LMI Çözüldü.');
    
    Q_val = double(Q);
    L_val = double(L);
    
    K_Hinf = L_val / Q_val; 
    gamma_opt = double(gamma);
    
    fprintf('Optimal H-Infinity Normu (gamma): %.4f\n', gamma_opt);
else
    disp('HATA: H-Infinity LMI Çözümü Başarısız Oldu.');
    disp(yalmiperror(sol.problem));
    K_Hinf = [];
    gamma_opt = NaN;
end
end