function [K, rho, feas] = solve_continuous_lmi_decentralized(Atot, Bdec, Cdec, N)
    % solve_continuous_lmi_decentralized: Computes a decentralized state-feedback 
    % gain K where each subsystem is controlled independently.
    
    [Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
    
    disp('--- Solving Continuous-Time Decentralized LMI ---');
    fprintf('Target: Spectral Abscissa (max(Re(lambda))) < 0\n');
    fprintf('Information Structure: Strictly Block-Diagonal (No Communication)\n');
    
    yalmip clear
    
    %% 1. Define Block-Diagonal Lyapunov Matrix (Y)
    % For decentralized stability, Y must be block-diagonal to ensure 
    % a separate Lyapunov function for each subsystem.
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    %% 2. Define Block-Diagonal Auxiliary Matrix (L)
    % L = K*Y. To ensure K is decentralized (block-diagonal), L must also
    % be restricted to a block-diagonal structure.
    L = sdpvar(mtot, ntot, 'full');
    L_constraints = [];
    row_ptr = 0; col_ptr = 0;
    
    for i = 1:N
        for j = 1:N
            % If it is not a diagonal block, set the elements to zero
            if i ~= j
                L_constraints = [L_constraints, L(row_ptr+1:row_ptr+m_sub(i), col_ptr+1:col_ptr+n_sub(j)) == 0];
            end
            col_ptr = col_ptr + n_sub(j);
        end
        row_ptr = row_ptr + m_sub(i); 
        col_ptr = 0;
    end
    
    %% 3. LMI Stability Constraints
    % Continuous-time Lyapunov Stability: A*Y + Y*A' + B*L + L'*B' < 0
    tol = 1e-2; % Numerical tolerance for strict inequality
    
    LMI_stability = (Atot*Y + Y*Atot' + Btot*L + L'*Btot') <= -tol*eye(ntot);
    LMI_positivity = Y >= tol*eye(ntot);
    
    constraints = [LMI_stability, LMI_positivity, L_constraints];
    
    %% 4. Optimization
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    diagnostics = optimize(constraints, [], options);
    
    feas = diagnostics.problem; % 0 indicates Success
    
    %% 5. Result Processing and Gain Recovery
    if feas == 0
        disp('Success: Continuous-Time Decentralized LMI satisfied.');
        
        % Recover Decentralized Gain: K = L * Y^-1
        K = double(L) / double(Y);
        
        % Calculate Closed-Loop Eigenvalues
        cl_eigenvalues = eig(Atot + Btot*K);
        rho = max(real(cl_eigenvalues));
        
        fprintf('Closed-Loop Spectral Abscissa (rho): %.4f\n', rho);
        
        if rho < 0
            disp('Congratulations! The system is stabilized under a decentralized structure.');
        end
    else
        disp('ERROR: Continuous-Time Decentralized LMI Optimization Failed.');
        fprintf('Error Code: %d (%s)\n', feas, yalmiperror(feas));
        K = []; rho = [];
    end
end