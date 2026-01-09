function [K, rho, feas] = solve_discrete_lmi_distributed(Ftot, Gdec, Hdec, N, ContStruc)
    % solve_discrete_lmi_distributed: Computes a state-feedback gain K 
    % in discrete-time subject to a specific communication topology.
    %
    % Ftot: Global discrete-time system matrix (F)
    % Gdec: Cell array of discrete-time input matrices {G1, G2, G3}
    % Hdec: Cell array of discrete-time output/state matrices {H1, H2, H3}
    % ContStruc: NxN binary matrix (1 = communication allowed, 0 = no link)

    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    
    disp('--- Solving Discrete-Time Distributed LMI ---');
    fprintf('Information Structure: %dx%d ContStruc matrix applied.\n', N, N);
    
    yalmip clear

    %% 1. Define Block-Diagonal Lyapunov Matrix (Y)
    % To maintain localized stability conditions, Y is kept block-diagonal.
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    %% 2. Define Structured Auxiliary Matrix (L)
    % L = K*Y. The sparsity of L matches the communication structure (ContStruc).
    L = sdpvar(mtot, ntot, 'full');
    L_constraints = [];
    row_p = 0;
    
    for i = 1:N
        col_p = 0;
        for j = 1:N
            % If ContStruc(i,j) is 0, subsystem i cannot use states from subsystem j.
            if ContStruc(i,j) == 0
                L_constraints = [L_constraints, L(row_p+1:row_p+m_sub(i), col_p+1:col_p+n_sub(j)) == 0];
            end
            col_p = col_p + n_sub(j);
        end
        row_p = row_p + m_sub(i);
    end

    %% 3. Discrete-Time LMI (Schur Complement Form)
    % Linearized Stability Condition: [Y, (F*Y + G*L)'; (F*Y + G*L), Y] > 0
    epsilon = 1e-6;
    M = [Y, (Ftot*Y + Gtot*L)'; (Ftot*Y + Gtot*L), Y];
    
    constraints = [M >= epsilon*eye(2*ntot), Y >= epsilon*eye(ntot), L_constraints];
    
    %% 4. Optimization
    options = sdpsettings('verbose', 0, 'solver', 'sedumi');
    sol = optimize(constraints, [], options);
    
    feas = sol.problem; % 0 indicates Success
    
    %% 5. Result Processing and Gain Recovery
    if feas == 0
        disp('Success: Discrete-Time Distributed LMI satisfied.');
        
        % Recover Structured Gain: K = L * Y^-1
        K = double(L) / double(Y);
        
        % Calculate Closed-Loop Spectral Radius
        cl_eigenvalues = eig(Ftot + Gtot*K);
        rho = max(abs(cl_eigenvalues));
        
        fprintf('Closed-Loop Spectral Radius (rho): %.4f\n', rho);
    else
        disp('ERROR: Discrete-Time Distributed LMI Optimization Failed.');
        K = []; rho = [];
    end
end