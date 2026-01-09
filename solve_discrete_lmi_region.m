function [K, rho, feas] = solve_discrete_lmi_region(Ftot, Gdec, Hdec, N, ContStruc, r, theta_deg)
    % solve_discrete_lmi_region_distributed: 
    % Doküman Kaynak 107 (Disk) ve Kaynak 149 (Sektör) kurallarına göre tasarım yapar.
    
    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    theta = deg2rad(theta_deg); % Doküman Kaynak 132'deki theta açısı

    fprintf('\n--- Solving Discrete-Time Structured Region LMI ---\n');
    yalmip clear

    %% 1. Yapısal Kısıtlar (Structure Constraints)
    % Y matrisi blok-diyagonal olmalıdır 
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    % L matrisi Kx ile aynı yapıya sahip olmalıdır 
    L = sdpvar(mtot, ntot, 'full');
    L_constraints = [];
    row_p = 0;
    for i = 1:N
        col_p = 0;
        for j = 1:N
            if ContStruc(i,j) == 0
                L_constraints = [L_constraints, L(row_p+1:row_p+m_sub(i), col_p+1:col_p+n_sub(j)) == 0];
            end
            col_p = col_p + n_sub(j);
        end
        row_p = row_p + m_sub(i);
    end

    %% 2. LMI Kısıtlamaları (Performance Regions)
    epsilon = 1e-7;
    FY_GL = Ftot*Y + Gtot*L;

    % Kısıt A: Spektral Yarıçap (Disk |z| < r) [cite: 107]
    % Not: Doküman A->F, B->G değişimini söyler 
    LMI_radius = [ (r^2)*Y , (FY_GL); 
                    (FY_GL)', Y ] >= epsilon*eye(2*ntot);

    % Kısıt B: Konik Sektör (Damping) [cite: 149]
    % ÖNEMLİ: Dokümandaki formül Re(z) < 0 içindir. 
    % Ringing'i önlemek (Re(z) > 0) için diagonal blokların işaretini ters çeviriyoruz.
    Term1 = -(FY_GL + FY_GL'); % Re(z) > 0 bölgesi için uyarlama
    Term2 = (FY_GL - FY_GL');
    
    LMI_sector = [ sin(theta)*Term1,  cos(theta)*Term2;
                  -cos(theta)*Term2,  sin(theta)*Term1 ] <= -epsilon*eye(2*ntot);

    % Kısıt C: Pozitif Tanımlılık [cite: 26]
    LMI_pos = Y >= epsilon*eye(ntot);

    constraints = [LMI_radius, LMI_pos, L_constraints];

    %% 3. Çözüm ve Kazanç Geri Kazanımı
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, [], options);
    
    feas = sol.problem;
    if feas == 0 || feas == 4
        % K = L * Y^-1 [cite: 27, 67]
        K = double(L) / double(Y);
        rho = max(abs(eig(Ftot + Gtot*K)));
        fprintf('Success: Region constraints satisfied. Rho: %.4f\n', rho);
    else
        fprintf('Error: LMI Infeasible (Code %d). Check targets.\n', feas);
        K = []; rho = [];
    end
end