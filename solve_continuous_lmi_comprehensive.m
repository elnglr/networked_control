function [K, rho, feas] = solve_continuous_lmi_comprehensive(Atot, Bdec, Cdec, N, ContStruc, alpha, zeta, u_min, u_ss, u_max, x0_norm, aL, aY)
    % Atot, Bdec, Cdec, N, ContStruc: Sistem tanımları ve bilgi yapısı
    % alpha, zeta: Bölgesel özdeğer yerleşimi hedefleri
    % u_min, u_ss, u_max, x0_norm: Giriş kısıtları (Opsiyonel, boş [] bırakılabilir)
    % aL, aY: Kontrol çabası minimizasyonu ağırlıkları (Opsiyonel)

    [Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
    theta = acos(zeta); 
    yalmip clear
    epsilon = 1e-7;

    %% 1. Karar Değişkenleri
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric')); % Blok-diyagonal Lyapunov matrisi
    end
    L = sdpvar(mtot, ntot, 'full'); % Yardımcı matris (L = KY)
    kL = sdpvar(1,1); % Kontrol çabası üst sınırı
    kY = sdpvar(1,1); % Durum uzayı üst sınırı

    %% 2. Temel Performans ve Yapısal Kısıtlar
    AY_BL = Atot*Y + Btot*L;
    
    % Yapısal Bilgi Kısıtı (ContStruc)
    L_constraints = apply_structure_constraints(L, N, m_sub, n_sub, ContStruc); %
    
    % Bölgesel Özdeğer Yerleşimi (Performance)
    LMI_speed = (AY_BL + AY_BL') + 2*alpha*Y <= -epsilon*eye(ntot); % Hız kısıtı
    LMI_damping = [ sin(theta)*(AY_BL + AY_BL'),  cos(theta)*(AY_BL - AY_BL');
                   -cos(theta)*(AY_BL - AY_BL'),  sin(theta)*(AY_BL + AY_BL') ] <= -epsilon*eye(2*ntot); % Sönümleme
    
    constraints = [LMI_speed, LMI_damping, Y >= epsilon*eye(ntot), L_constraints];

    %% 3. Giriş Kısıtlamaları (Hard Constraints - Opsiyonel)
    if nargin >= 11 && ~isempty(u_min) && ~isempty(u_max) && ~isempty(x0_norm)
        delta_u_upper = u_max - u_ss;
        delta_u_lower = u_ss - u_min;
        mu = min(delta_u_upper, delta_u_lower); 
        
        for i = 1:mtot
            constraints = [constraints, [mu(i)^2, L(i,:); L(i,:)', Y] >= 0]; % Giriş elipsoidi
        end
        % constraints = [constraints, Y >= (x0_norm^2) * eye(ntot)]; % Başlangıç koşulu güvenliği
    end

    %% 4. Kontrol Çabası Minimizasyonu (Effort Reduction - Opsiyonel)
    LMI_effort = [];
    if nargin >= 13 && ~isempty(aL) && ~isempty(aY)
        LMI_effort = [[kL*eye(mtot), L; L', eye(ntot)] >= epsilon*eye(mtot + ntot), ... %
                      [kY*eye(ntot), eye(ntot); eye(ntot), Y] >= epsilon*eye(2*ntot)]; %
        Cost = aL*kL + aY*kY; %
    else
        Cost = []; % Sadece fizibilite problemi
    end
    
    constraints = [constraints, LMI_effort];

    %% 5. Optimizasyon ve Çözüm
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, Cost, options);
    feas = sol.problem;

    if sol.problem == 0 || sol.problem == 4
    K = double(L) / double(Y); % K matrisini hesapla
    rho = max(real(eig(Atot + Btot*K)));
    
    % Gerçek Norm ve Teorik Sınır
    actual_norm = norm(K, 2); 
    theoretical_bound = sqrt(double(kL)) * double(kY); %
    
    fprintf('Success!\n');
    fprintf('  -> Actual norm(K,2): %.4f\n', actual_norm);
    fprintf('  -> Theoretical Bound: %.4f\n', theoretical_bound);
    else
        K = []; rho = [];
        fprintf('Infeasible problem for the given constraints.\n');
    end

end