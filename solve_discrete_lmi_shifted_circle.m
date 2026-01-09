function [K, rho, feas] = solve_discrete_lmi_shifted_circle(Ftot, Gdec, Hdec, N, ContStruc, alpha_shift, r_damp)
% solve_discrete_lmi_shifted_circle: Kutupları z-düzleminde kaydırılmış bir 
% çember içine hapseder (LMI Region: Disk).
%
% alpha_shift: Çemberin merkezi (Reel eksende, örn: 0.3)
% r_damp: Çemberin yarıçapı (örn: 0.3)

[Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);

yalmip clear
epsilon = 1e-8;

%% 1. Lyapunov Matrisi (Y) ve Karar Değişkeni (L)
Y = [];
for i=1:N
    Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
end

L = sdpvar(mtot, ntot, 'full');

% Bilgi Yapısı Kısıtları (ContStruc)
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

%% 2. Shifted Circle LMI (Schur Complement uygulanmış hali)
% (F_cl - alpha*I)*Y*(F_cl - alpha*I)' - r^2*Y < 0
% F_tilde = F - alpha*I;
F_tilde = Ftot - alpha_shift * eye(ntot);

% Görseldeki matris yapısı: [r^2*Y , F_tilde*Y + G*L ; (F_tilde*Y + G*L)' , Y] > 0
LMI_shifted = [ (r_damp^2)*Y ,  F_tilde*Y + Gtot*L ;
                (F_tilde*Y + Gtot*L)', Y ] >= epsilon * eye(2*ntot);

constraints = [LMI_shifted, Y >= epsilon*eye(ntot), L_constraints];

%% 3. Çözücü Ayarları ve Çözüm
options = sdpsettings('solver', 'sedumi', 'verbose', 0);
sol = optimize(constraints, [], options);
feas = sol.problem;

if feas == 0 || feas == 4
    K = double(L) / double(Y);
    rho = max(abs(eig(Ftot + Gtot*K)));
    fprintf('Success! Shifted Circle: Max |z| = %.4f\n', rho);
else
    K = []; rho = [];
    fprintf('Infeasible! Bölgeyi genişletmeyi deneyin.\n');
end
end