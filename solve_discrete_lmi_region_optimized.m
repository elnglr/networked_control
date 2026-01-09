function [K, rho, feas] = solve_discrete_lmi_region_optimized(Ftot, Gdec, Hdec, N, ContStruc, r)
    % [cite: 107] Spektral Yarıçap + [cite: 79-95] Kontrol Çabası Minimizasyonu
    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    yalmip clear

    %% 1. Değişkenler
    P = [];
    for i=1:N
        P = blkdiag(P, sdpvar(n_sub(i), n_sub(i), 'symmetric')); % [cite: 260]
    end
    L = sdpvar(mtot, ntot, 'full');
    kL = sdpvar(1,1); 
    kP = sdpvar(1,1);

    %% 2. Kısıtlar
    epsilon = 1e-7;
    FY_GL = Ftot*P + Gtot*L;
    L_constraints = apply_structure_constraints(L, N, m_sub, n_sub, ContStruc);

    % Spektral Yarıçap Kısıtı (Radius)
    LMI_radius = [ (r^2)*P , (FY_GL); 
                    (FY_GL)', P ] >= epsilon*eye(2*ntot); % [cite: 107]

    % Kontrol Çabası Kısıtları [cite: 86, 89]
    LMI_kL = [kL*eye(mtot), L; L', eye(ntot)] >= epsilon*eye(mtot + ntot);
    LMI_kP = [kP*eye(ntot), eye(ntot); eye(ntot), P] >= epsilon*eye(2*ntot);

    constraints = [LMI_radius, LMI_kL, LMI_kP, P >= epsilon*eye(ntot), L_constraints];

    %% 3. Optimizasyon
    Cost = 0.01*kL + 10*kP; ,
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, Cost, options);

    feas = sol.problem;
    if feas == 0 || feas == 4
        K = double(L) / double(P); % [cite: 67, 82]
        rho = max(abs(eig(Ftot + Gtot*K)));
        fprintf('Success! Discrete Optimal Effort Bound: %.4f\n', sqrt(double(kL))*double(kP));
    else
        K = []; rho = [];
    end
end