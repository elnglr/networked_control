function [K, rho, feas] = solve_discrete_lmi_decentralized(Ftot, Gdec, Hdec, N)
    % solve_discrete_lmi_decentralized: Computes a decentralized state-feedback 
    % gain K in discrete-time where each subsystem is controlled independently.
    %
    % Ftot: Discrete-time global system matrix (F)
    % Gdec: Cell array of discrete-time input matrices {G1, G2, G3}
    % Hdec: Cell array of discrete-time output/state matrices {H1, H2, H3}
    
    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    
    disp('--- Solving Discrete-Time Decentralized LMI ---');
    fprintf('Target: Closed-loop spectral radius (rho) < 1.0000\n');
    
    yalmip clear

    %% 1. Define Block-Diagonal Lyapunov Matrix (Y)
    % For decentralized stability, Y must be block-diagonal to ensure 
    % localized stability certificates for each subsystem.
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    %% 2. Define Block-Diagonal Auxiliary Matrix (L)
    % L = K*Y. To ensure K is decentralized (block-diagonal), L must also
    % be restricted to a block-diagonal structure.
    L = sdpvar(mtot, ntot, 'full');
    L_constraints = [];
    row_p = 0; col_p = 0;
    
    for i = 1:N
        for j = 1:N
            % If it is not a diagonal block, set the elements to zero
            if i ~= j
                L_constraints = [L_constraints, L(row_p+1:row_p+m_sub(i), col_p+1:col_p+n_sub(j)) == 0];
            end
            col_p = col_p + n_sub(j);
        end
        row_p = row_p + m_sub(i); 
        col_p = 0;
    end

    %% 3. Discrete-Time LMI (Schur Complement Form)
    % Stability condition: F_cl * Y * F_cl' - Y < 0
    % Linearized via Schur Complement: [Y, (F*Y + G*L)'; (F*Y + G*L), Y] > 0
    epsilon = 1e-6;
    M = [Y, (Ftot*Y + Gtot*L)'; (Ftot*Y + Gtot*L), Y];
    
    constraints = [M >= epsilon*eye(2*ntot), Y >= epsilon*eye(ntot), L_constraints];
    
    %% 4. Optimization
    options = sdpsettings('verbose', 0, 'solver', 'sedumi');
    sol = optimize(constraints, [], options);
    
    feas = sol.problem; % 0 indicates Success
    
    %% 5. Result Processing and Gain Recovery
    if feas == 0
        disp('Success: Discrete-Time Decentralized LMI satisfied.');
        
        % Recover Decentralized Gain: K = L * Y^-1
        K = double(L) / double(Y);
        
        % Calculate Closed-Loop Spectral Radius
        cl_eigenvalues = eig(Ftot + Gtot*K);
        rho = max(abs(cl_eigenvalues));
        
        fprintf('Closed-Loop Spectral Radius (rho): %.4f\n', rho);
        
        if rho < 1
            disp('Congratulations! The system is stabilized in discrete-time under a decentralized structure.');
        end
    else
        disp('ERROR: Discrete-Time Decentralized LMI Optimization Failed.');
        K = []; rho = [];
    end
end