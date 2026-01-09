function [K, rho, feas] = solve_continuous_lmi_constrained(Atot, Bdec, Cdec, N, ContStruc, alpha, zeta, u_min, u_ss, u_max, x0_norm)
% u_min, u_ss, u_max: Tablo 2'deki değerler (vektör olarak)
% x0_norm: Beklenen maksimum başlangıç koşulu normu (örn: 0.1 veya 1)

[Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
theta = acos(zeta); 

% --- Giriş Sınırlarını Hesapla (Shifted) ---
% Lineer model için izin verilen sapma limitleri
delta_u_upper = u_max - u_ss;
delta_u_lower = u_ss - u_min;
% Simetrik sınır matrisi (Mu)
mu = min(delta_u_upper, delta_u_lower); 

yalmip clear
epsilon = 1e-7;

%% 1. Lyapunov ve Yardımcı Matrisler
Y = [];
for i=1:N
    Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
end
L = sdpvar(mtot, ntot, 'full');

%% 2. Bilgi Yapısı Kısıtları (ContStruc)
L_constraints = [];
row_ptr = 0;
for i = 1:N
    col_ptr = 0;
    for j = 1:N
        if ContStruc(i,j) == 0
            L_constraints = [L_constraints, L(row_ptr+1:row_ptr+m_sub(i), col_ptr+1:col_ptr+n_sub(j)) == 0];
        end
        col_ptr = col_ptr + n_sub(j);
    end
    row_ptr = row_ptr + m_sub(i);
end

%% 3. LMI Bölge Kısıtları
AY_BL = Atot*Y + Btot*L;
LMI_speed = (AY_BL + AY_BL') + 2*alpha*Y <= -epsilon*eye(ntot);
LMI_damping = [ sin(theta)*(AY_BL + AY_BL'),  cos(theta)*(AY_BL - AY_BL');
               -cos(theta)*(AY_BL - AY_BL'),  sin(theta)*(AY_BL + AY_BL') ] <= -epsilon*eye(2*ntot);

%% 4. Giriş Kısıtlaması (Input Constraints)
% Invariant Ellipsoid: u_i^2 <= mu_i^2 koşulu için LMI
% Schur Complement ile: [mu_i^2 * I , L_i ; L_i' , Y] >= 0
LMI_input = [];
for i = 1:mtot
    LMI_input = [LMI_input, [mu(i)^2, L(i,:); L(i,:)', Y] >= 0];
end

% Başlangıç koşulu elipsoidi (x0' * Y^-1 * x0 <= 1 olması için)
% Eğer x0 normu biliniyorsa (örn ||x0|| <= rho_x), Y >= rho_x^2 * I olmalı
LMI_init = Y >= (x0_norm^2) * eye(ntot);

constraints = [LMI_speed, LMI_damping, LMI_input, LMI_init, L_constraints];

%% 5. Çözüm
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
sol = optimize(constraints, [], options);
feas = sol.problem;

if feas == 0 || feas == 4
    K = double(L) / double(Y);
    rho = max(real(eig(Atot + Btot*K)));
else
    K = []; rho = [];
end
end