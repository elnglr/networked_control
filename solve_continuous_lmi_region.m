function [K, rho, feas] = solve_continuous_lmi_region(Atot, Bdec, Cdec, N, ContStruc, alpha, zeta)
% solve_lmi_region_distributed_CT: Places poles in a region while respecting 
% a specific information structure (Centralized, Decentralized, or Distributed).
%
% Atot: Global system matrix (A)
% Bdec, Cdec: Subsystem matrix cells
% ContStruc: NxN binary matrix (1 = link, 0 = no link)
% alpha: Decay rate (Re(s) < -alpha)
% zeta: Damping ratio (zeta >= cos(theta))

% Prepare dimensions from sub-cells
[Btot, ntot, mtot, m_sub, n_sub] = prepare_dimensions(Bdec, Cdec, N);
theta = acos(zeta); 

fprintf('\n--- Solving Structured Continuous-Time LMI Region ---\n');
fprintf('Structure: Distributed/Decentralized (ContStruc applied)\n');
fprintf('Target: Re(s) < -%.2f and Damping > %.2f\n', alpha, zeta);

yalmip clear

%% 1. Define Block-Diagonal Lyapunov Matrix (Y)
% To maintain independent Lyapunov certificates, Y must be block-diagonal.
Y = [];
for i=1:N
    Y = blkdiag(Y, sdpvar(n_sub(i), n_sub(i), 'symmetric'));
end

%% 2. Define Structured Auxiliary Matrix (L)
% L = K*Y. We force L(i,j) = 0 where ContStruc(i,j) = 0.
L = sdpvar(mtot, ntot, 'full');
L_constraints = [];
row_ptr = 0;
for i = 1:N
    col_ptr = 0;
    for j = 1:N
        if ContStruc(i,j) == 0
            % If no communication from subsystem j to i, block must be zero
            L_constraints = [L_constraints, L(row_ptr+1:row_ptr+m_sub(i), col_ptr+1:col_ptr+n_sub(j)) == 0];
        end
        col_ptr = col_ptr + n_sub(j);
    end
    row_ptr = row_ptr + m_sub(i);
end

%% 3. Define Regional LMI Constraints
epsilon = 1e-7;
AY_BL = Atot*Y + Btot*L;

% Constraint 1: Speed (Alpha-stability) -> Re(s) < -alpha
LMI_speed = (AY_BL + AY_BL') + 2*alpha*Y <= -epsilon*eye(ntot);

% Constraint 2: Damping (Conic Sector) -> zeta
LMI_damping = [ sin(theta)*(AY_BL + AY_BL'),  cos(theta)*(AY_BL - AY_BL');
               -cos(theta)*(AY_BL - AY_BL'),  sin(theta)*(AY_BL + AY_BL') ] <= -epsilon*eye(2*ntot);

% Positive Definiteness
LMI_pos = Y >= epsilon*eye(ntot);
if zeta ==0
constraints = [LMI_speed, LMI_pos, L_constraints];

else

constraints = [LMI_speed, LMI_damping, LMI_pos, L_constraints];

end

%% 4. Optimization
options = sdpsettings('solver', 'sedumi', 'verbose', 0, 'sedumi.eps', 1e-10);
sol = optimize(constraints, [], options);

feas = sol.problem;

%% 5. Gain Recovery and Verification
if feas == 0
    K = double(L) / double(Y);
    rho = max(real(eig(Atot + Btot*K)));
    fprintf('Success! Structured LMI Region satisfied.\n');
    fprintf('Achieved max(Re(lambda)): %.4f\n', rho);

elseif feas == 4
    disp('Warning: The problem is Numerically Ill-Conditioned.');
    K = double(L) / double(Y);
    rho = max(real(eig(Atot + Btot*K)));
    fprintf('Check rho: max(Re(lambda)): %.4f\n',rho);
else
    disp('ERROR: The problem is INFEASIBLE.');
    disp('Reason: The communication structure is too restrictive for the desired performance.');
    K = []; rho = [];



end
end