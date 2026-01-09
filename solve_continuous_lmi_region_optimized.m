function [K, rho, feas] = solve_continuous_lmi_region_optimized(Atot, Bdec, Cdec, N, ContStruc, alpha, zeta)
    % [cite: 76-117] Bölgesel Özdeğer Yerleşimi + [cite: 79-95] Kontrol Çabası Minimizasyonu
    [Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
    theta = acos(zeta); 
    yalmip clear

    %% 1. Değişkenler
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric')); % [cite: 260]
    end
    L = sdpvar(mtot, ntot, 'full');
    kL = sdpvar(1,1); % [cite: 85]
    kY = sdpvar(1,1); % [cite: 87]

    %% 2. Kısıtlar (Constraints)
    epsilon = 1e-7;
    AY_BL = Atot*Y + Btot*L;
    
    % Yapısal Kısıt (ContStruc)
    L_constraints = apply_structure_constraints(L, N, m_sub, n_sub, ContStruc); % [cite: 259]

    % Bölgesel Kısıtlar (Performance)
    LMI_speed = (AY_BL + AY_BL') + 2*alpha*Y <= -epsilon*eye(ntot); % [cite: 102]
    LMI_damping = [ sin(theta)*(AY_BL + AY_BL'),  cos(theta)*(AY_BL - AY_BL');
                   -cos(theta)*(AY_BL - AY_BL'),  sin(theta)*(AY_BL + AY_BL') ] <= -epsilon*eye(2*ntot); % [cite: 149]

    % Kontrol Çabası Kısıtları (Effort Reduction)
    LMI_kL = [kL*eye(mtot), L; L', eye(ntot)] >= epsilon*eye(mtot + ntot); % [cite: 86]
    LMI_kY = [kY*eye(ntot), eye(ntot); eye(ntot), Y] >= epsilon*eye(2*ntot); % [cite: 89]

    constraints = [LMI_speed, LMI_damping, LMI_kL, LMI_kY, Y >= epsilon*eye(ntot), L_constraints];

    %% 3. Optimizasyon (Maliyet Fonksiyonu ile)
    % Doküman sayfa 79'daki önerilen ağırlıklar: aL=0.01, aY=10 
    Cost = 10*kL + 1*kY;
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, Cost, options);

    feas = sol.problem;
    if feas == 0 || feas == 4
        K = double(L) / double(Y); % [cite: 27]
        rho = max(real(eig(Atot + Btot*K)));
        fprintf('Success! Optimal Effort Bound (||Kx||): %.4f\n', sqrt(double(kL))*double(kY));
    else
        K = []; rho = [];
    end
end