function [K, rho, feas] = solve_discrete_lmi_RHP_only(Ftot, Gdec, Hdec, N, ContStruc, r)
    % solve_discrete_lmi_RHP_only: Forces poles to be inside radius 'r' 
    % AND strictly in the Right-Half Plane (Re(z) >= 0) to avoid ringing.
    
    [Gtot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Gdec, Hdec, N);
    
    fprintf('\n--- Solving DT LMI: Radius < %.2f & Re(z) >= 0 (No Ringing) ---\n', r);
    yalmip clear

    %% 1. Define Variables
    Y = [];
    for i=1:N
        Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
    end
    
    L = sdpvar(mtot, ntot, 'full');
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

    %% 2. LMI Constraints
    epsilon = 1e-7;
    FY_GL = Ftot*Y + Gtot*L;

    % Constraint 1: Circular Disk (Radius < r)
    LMI_circ = [ -r*Y,         FY_GL; 
                  FY_GL',      -r*Y ] <= -epsilon*eye(2*ntot);

    % Constraint 2: Right-Half Plane Constraint (Re(z) >= 0)
    % Formula: (Fcl*Y + Y*Fcl') >= 0
    LMI_RHP = (FY_GL + FY_GL') >= epsilon*eye(ntot);

    constraints = [LMI_circ, LMI_RHP, Y >= epsilon*eye(ntot), L_constraints];

    %% 3. Optimization
    options = sdpsettings('solver', 'sedumi', 'verbose', 0);
    sol = optimize(constraints, [], options);
    
    feas = sol.problem;

    %% 4. Results
    if feas == 0
        K = double(L) / double(Y);
        cl_eig = eig(Ftot + Gtot*K);
        rho = max(abs(cl_eig));
        fprintf('Success! Max Radius: %.4f | Min Real Part: %.4f\n', rho, min(real(cl_eig)));
    else
        disp('ERROR: Infeasible. Re(z) >= 0 constraint might be too tight for this structure.');
        K = []; rho = [];
    end
end