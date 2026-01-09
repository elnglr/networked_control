function K_disc = solve_discrete_lmi_centralized(F, G)
% solve_discrete_lmi_centralized: Solves the standard Discrete-Time 
% stabilization LMI (Poles |lambda| < 1).
%
% F: Discrete-time system matrix
% G: Discrete-time input matrix

n = size(F, 1); 
m = size(G, 2); 

disp('--- Solving 1. Discrete-Time Stabilization LMI (Poles < 1) ---');

%% 1. Define YALMIP Variables
% Q is the inverse of the Lyapunov matrix (Q = P^-1)
% L is the auxiliary variable (L = K * Q) used for linearization
Q = sdpvar(n, n, 'symmetric'); 
L = sdpvar(m, n, 'full');

%% 2. LMI Constraints (Schur Complement Form)
% The stability condition is F_cl' * P * F_cl - P < 0.
% Using the Schur Complement and the substitution L = KQ, 
% we obtain the linearized LMI form:
% [Q , (F*Q + G*L)'; (F*Q + G*L), Q] > 0

epsilon = 1e-6; 

% Linearized Schur Complement representation:
LMI_stability = [Q, (F*Q + G*L)'; 
                 (F*Q + G*L), Q] >= epsilon*eye(2*n);
         
constraints = [LMI_stability, Q >= epsilon*eye(n)];

%% 3. Optimization
% Solving the feasibility problem using the SDPT3 solver
options = sdpsettings('verbose', 0, 'solver', 'sdpt3'); 
sol = optimize(constraints, [], options);

%% 4. Results and Gain Recovery
if sol.problem == 0
    disp('Success: Discrete-Time Stabilization LMI satisfied.');
    
    Q_val = double(Q);
    L_val = double(L);
    
    % Recover the feedback gain: K = L * Q^-1
    K_disc = L_val / Q_val; 
    
    % Verify Closed-Loop Stability (Spectral Radius)
    F_cl = F + G * K_disc;
    rho = max(abs(eig(F_cl)));
    
    fprintf('Closed-Loop Spectral Radius (rho): %.4f\n', rho);
    
    if rho < 1
        disp('System successfully stabilized in Discrete-Time.');
    end
else
    disp('ERROR: Discrete-Time Stabilization LMI Failed.');
    disp(yalmiperror(sol.problem));
    K_disc = [];
end
end