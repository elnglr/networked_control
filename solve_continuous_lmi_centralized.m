function K_cont = solve_continuous_lmi_centralized(A, B)
% solve_continuous_lmi_centralized: Computes a stabilizing state-feedback 
% gain K using Linear Matrix Inequalities (LMI).
% Requires YALMIP and an SDP solver (e.g., SDPT3, SeDuMi, or MOSEK).
% 
% A: n x n system matrix (12x12 for your plant)
% B: n x m input matrix (12x6 for your plant)

n = size(A, 1); % State dimension
m = size(B, 2); % Input dimension

disp('--- Solving Continuous-Time Centralized Stability LMI ---');

%% 1. Define YALMIP Variables
% Q is the inverse of the Lyapunov matrix (Q = P^-1). It must be symmetric and Q > 0.
Q = sdpvar(n, n, 'symmetric'); 

% L is an auxiliary variable defined as L = K * Q to linearize the BMI.
L = sdpvar(m, n, 'full');

%% 2. Define LMI Constraints
% Stability Constraint (Lyapunov Inequality): 
% In the continuous-time domain, the closed-loop system is stable if:
% (A + BK)Q + Q(A + BK)' < 0  ==>  AQ + Q'A + BL + L'B' < 0
% We add a small epsilon (1e-6) for strict inequality and numerical robustness.
epsilon = 1e-6;

LMI_stability = (A*Q + Q*A' + B*L + L'*B') <= -epsilon*eye(n);

% Positive Definiteness Constraint: Q > 0
LMI_positivity = Q >= epsilon*eye(n);

constraints = [LMI_stability, LMI_positivity];

%% 3. Optimization
% Since we are looking for a feasible solution (Stability), the objective function is empty.
options = sdpsettings('verbose', 0, 'solver', 'sdpt3'); 
sol = optimize(constraints, [], options);

%% 4. Post-Processing and Gain Recovery
if sol.problem == 0
    disp('Success: Continuous-Time Centralized Stability LMI satisfied.');
    
    % Retrieve numerical values
    Q_val = double(Q);
    L_val = double(L);
    
    % Recover the Control Gain: K = L * Q^-1
    K_cont = L_val / Q_val; 
    
    % Verify Closed-Loop Stability
    A_cl = A + B * K_cont;
    cl_eigenvalues = eig(A_cl);
    spectral_abscissa = max(real(cl_eigenvalues));
    
    fprintf('Closed-Loop Spectral Abscissa (max(Re(lambda))): %.4f\n', spectral_abscissa);
    
    if spectral_abscissa < 0
        disp('System successfully stabilized with Centralized K.');
    else
        disp('WARNING: Closed-loop not stable. Check numerical precision or system rank.');
    end
else
    disp('ERROR: Continuous-Time Centralized LMI Optimization Failed.');
    disp(yalmiperror(sol.problem));
    K_cont = [];
end
end