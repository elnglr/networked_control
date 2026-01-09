function [K, rho, feas] = solve_discrete_lmi_comprehensive(Ftot, Gdec, Hdec, N, ContStruc, alpha_shift, r_damp, u_min, u_ss, u_max, x0_norm, aL, aP)
    % Ftot, Gdec, Hdec, N, ContStruc: System and structure definitions
    % alpha_shift, r_damp: Shifted Circle parameters (Center and Radius)
    % u_min, u_ss, u_max, x0_norm: Input constraints (Table 2 data)
    % aL, aP: Effort minimization weights (e.g., aL=0.01, aP=10)

    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    yalmip clear
    epsilon = 1e-8;

    %% 1. Decision Variables
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric')); % Block-diagonal Lyapunov
    end
    L = sdpvar(mtot, ntot, 'full');
    kL = sdpvar(1,1); % Control gain bound
    kP = sdpvar(1,1); % State bound

    %% 2. Structural & Performance Constraints
    % Information Structure Constraints
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

    % Shifted Circle Region Constraint (Damping & Speed)
    F_tilde = Ftot - alpha_shift * eye(ntot);
    % [r^2*Y, F_tilde*Y + G*L; (F_tilde*Y + G*L)', Y] > 0
    LMI_shifted = [ (r_damp^2)*Y ,  F_tilde*Y + Gtot*L ;
                    (F_tilde*Y + Gtot*L)', Y ] >= epsilon * eye(2*ntot);

    constraints = [LMI_shifted, Y >= epsilon*eye(ntot), L_constraints];

    %% 3. Input Constraints (Hard Saturation Limits)
    if nargin >= 11 && ~isempty(u_min)
        delta_u_upper = u_max - u_ss;
        delta_u_lower = u_ss - u_min;
        mu = min(delta_u_upper, delta_u_lower); 
        
        for i = 1:mtot
            % [mu_i^2, L_i; L_i', Y] >= 0
            constraints = [constraints, [mu(i)^2, L(i,:); L(i,:)', Y] >= 0];
        end
        % Initial condition safety: Y >= x0_norm^2 * I
        constraints = [constraints, Y >= (x0_norm^2) * eye(ntot)];
    end

    %% 4. Effort Minimization (Optional Optimization)
    Cost = [];
    if nargin >= 13 && ~isempty(aL)
        % [kL*I, L; L', I] >= 0
        % [kP*I, I; I, Y] >= 0
        LMI_effort = [[kL*eye(mtot), L; L', eye(ntot)] >= epsilon*eye(mtot + ntot), ...
                      [kP*eye(ntot), eye(ntot); eye(ntot), Y] >= epsilon*eye(2*ntot)];
        constraints = [constraints, LMI_effort];
        Cost = aL*kL + aP*kP; % J = aL*kL + aP*kP
    end

    %% 5. Optimization & Results
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, Cost, options);
    feas = sol.problem;
    if sol.problem == 0 || sol.problem == 4
        K = double(L) / double(Y);
        rho = max(abs(eig(Ftot + Gtot*K)));
        actual_norm = norm(K, 2);
        theoretical_bound = sqrt(double(kL)) * double(kP);
        fprintf('Success! Max |z| = %.4f | Actual ||K||: %.4f | Bound: %.4f\n', ...
                rho, actual_norm, theoretical_bound);
    else
        K = []; rho = [];
        fprintf('Infeasible! Try relaxing r_damp or increasing x0_norm.\n');
    end
end