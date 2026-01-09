function [K, rho, feas] = solve_continuous_lmi_distributed(Atot, Bdec, Cdec, N, ContStruc)
    % solve_continuous_lmi_distributed: Computes a state-feedback gain K 
    % subject to a specific communication topology (Distributed Control).
    %
    % ContStruc: NxN binary matrix where ContStruc(i,j) = 1 means subsystem i 
    %            receives information from subsystem j.

    [Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
    
    disp('--- Solving Continuous-Time Distributed LMI ---');
    fprintf('Applied Communication Structure (ContStruc):\n');
    disp(ContStruc);

    yalmip clear
    
    %% 1. Define Block-Diagonal Lyapunov Matrix (Y)
    % For structured control, Y remains block-diagonal to ensure 
    % localized stability certificates for each subsystem.
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    %% 2. Define Structured Auxiliary Matrix (L)
    % L = K*Y. The structure of K (determined by ContStruc) is imposed on L.
    L = sdpvar(mtot, ntot, 'full');
    L_constraints = [];
    row_ptr = 0;
    
    for i = 1:N
        col_ptr = 0;
        for j = 1:N
            % If there is no communication from subsystem j to i, 
            % the corresponding block in the gain matrix must be zero.
            if ContStruc(i,j) == 0
                L_constraints = [L_constraints, L(row_ptr+1:row_ptr+m_sub(i), col_ptr+1:col_ptr+n_sub(j)) == 0];
            end
            col_ptr = col_ptr + n_sub(j);
        end
        row_ptr = row_ptr + m_sub(i);
    end

    %% 3. LMI Stability Constraints
    % Continuous-time Lyapunov Inequality: A*Y + Y*A' + B*L + L'*B' < 0
    tol = 1e-2; % Numerical tolerance for strict inequality
    
    LMI_stability = (Atot*Y + Y*Atot' + Btot*L + L'*Btot') <= -tol*eye(ntot);
    LMI_positivity = Y >= tol*eye(ntot);
    
    constraints = [LMI_stability, LMI_positivity, L_constraints];
    
    %% 4. Optimization
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, [], options);
    
    feas = sol.problem; % 0 indicates Success
    
    %% 5. Results Processing and Gain Recovery
    if feas == 0
        disp('Success: Continuous-Time Distributed LMI satisfied.');
        
        % Recover Structured Gain: K = L * Y^-1
        K = double(L) / double(Y);
        
        % Calculate Closed-Loop Eigenvalues
        cl_eigenvalues = eig(Atot + Btot*K);
        rho = max(real(cl_eigenvalues));
        
        fprintf('Closed-Loop Spectral Abscissa (rho): %.4f\n', rho);
    else
        disp('ERROR: Continuous-Time Distributed LMI Optimization Failed.');
        fprintf('Error Code: %d\n', feas);
        K = []; rho = [];
    end
end