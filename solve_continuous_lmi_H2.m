function [K, trace_S, feas] = solve_continuous_lmi_H2(A, Bdec, Cdec, N, ContStruc, Q_cost, R_cost, alpha, zeta, u_min, u_ss, u_max, x0_norm)
% solve_continuous_lmi_H2_comprehensive: 
% H2 Optimizasyonu + Bilgi Yapısı + Bölgesel Yerleşim + Giriş Kısıtları

[Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
theta = acos(zeta);
yalmip clear
epsilon = 1e-7;

% H2 Sistemi için Ağırlık Matrislerini Hazırla
Bw = eye(ntot); % Tüm durumlara etki eden dış gürültü
Cz = [sqrtm(Q_cost); zeros(mtot, ntot)]; 
Du = [zeros(ntot, mtot); sqrtm(R_cost)];
p_cost = size(Cz, 1);

%% 1. Karar Değişkenleri
Q = [];
for i=1:N
    Q = blkdiag(Q, sdpvar(n_sub(i), n_sub(i), 'symmetric')); % Blok-diyagonal Lyapunov
end
L = sdpvar(mtot, ntot, 'full');
S = sdpvar(p_cost, p_cost, 'symmetric'); % H2 performans değişkeni

%% 2. Kısıtlar (Constraints)
% A. Yapısal Bilgi Kısıtı (ContStruc)
L_constraints = apply_structure_constraints(L, N, m_sub, n_sub, ContStruc); %

% B. H2 Performans LMI'ları
% LMI 1: Enerji Sınırı (AQ + QA' + BL + L'B' + Bw Bw' < 0)
AQ_BL = A*Q + Btot*L;
LMI_H2_1 = AQ_BL + AQ_BL' + Bw*Bw' <= -epsilon*eye(ntot);

% LMI 2: Çıktı Enerjisi Schur Complement [S > (CzQ + DuL)Q^-1(CzQ + DuL)']
LMI_H2_2 = [ S,                 Cz*Q + Du*L;
             (Cz*Q + Du*L)',    Q ] >= epsilon*eye(p_cost + ntot);

% C. Bölgesel Özdeğer Yerleşimi (Alpha & Zeta - Şekildeki yeşil bölge)
LMI_speed = (AQ_BL + AQ_BL') + 2*alpha*Q <= -epsilon*eye(ntot); %
LMI_damping = [ sin(theta)*(AQ_BL + AQ_BL'),  cos(theta)*(AQ_BL - AQ_BL');
               -cos(theta)*(AQ_BL - AQ_BL'),  sin(theta)*(AQ_BL + AQ_BL') ] <= -epsilon*eye(2*ntot); %

% D. Giriş Kısıtları (Table 2 - Saturation Avoidance)
LMI_input = [];
if nargin >= 13 && ~isempty(u_min)
    delta_u_upper = u_max - u_ss;
    delta_u_lower = u_ss - u_min;
    mu = min(delta_u_upper, delta_u_lower); 
    for i = 1:mtot
        LMI_input = [LMI_input, [mu(i)^2, L(i,:); L(i,:)', Q] >= 0]; %
    end
    % Başlangıç Koşulu Güvenliği
    % LMI_input = [LMI_input, Q >= (x0_norm^2) * eye(ntot)];
end

constraints = [LMI_H2_1, LMI_H2_2, LMI_input, L_constraints, Q >= epsilon*eye(ntot)];

%% 3. Optimizasyon (Trace(S) Minimize Edilir)
options = sdpsettings('verbose', 0, 'solver', 'sedumi'); 
sol = optimize(constraints, trace(S), options);

%% 4. Sonuçlar
feas = sol.problem;
if feas == 0 || feas == 4
    Q_val = double(Q);
    L_val = double(L);
    K = L_val / Q_val; 
    trace_S = double(trace(S));
    fprintf('H2 Success! Trace(S): %.4f | Actual Norm(K): %.4f\n', trace_S, norm(K,2));
else
    K = []; trace_S = NaN;
    fprintf('H2 Infeasible! Check constraints or region.\n');
end
end