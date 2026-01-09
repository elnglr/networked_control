
%%

clear;clc;
%%

OpenLoopAnaylsis;

n = size(A,1);  m = size(B,2);
%%

ContStruc_c      = ones(N,N);          % centralized
ContStruc_dec    = diag(ones(N,1));    % decentralized
ContStruc_string = [1 0 0;             % distributed string
                    1 1 0;
                    0 1 1];
ContStruc_3ten1e   = [1 0 0;             % distributed star (hub = 2)
                      1 1 0;
                      1 1 1];

idx1 = 1:4;  % Subsys1, First 4 States
idx2 = 5:8;  % Subsys2, Second 4 States
idx3 = 9:12; % Subsys3, Third 4 States
C1 = eye(4);
C2 = C1;
C3 = C1;


%%
% Tablo 2 Verileri
% [Ff1; Q1; Ff2; Q2; Fr; Q3]
u_min = [0; 0; 0; 0; 0; 0];
u_ss  = [8.33; 10; 0.5; 10; 66.2; 10];
u_max = [100000; 500000; 100000; 500000; 750000; 500000];

% Başlangıç sapması normu (Sistemi ne kadar sarsacağınız)
x0_norm = 1;
Qy1 = diag([1, 0, 0, 0.1]); 
Qy2 = diag([1, 0, 0, 0.1]);
Qy3 = diag([1, 0, 1000, 0]); % XB3'e yüksek ağırlık (10^3)



% Toplam Q_cost matrisini oluştur (C matrisleri üzerinden)
% C1, C2, C3 alt sistemlerin çıktı matrisleridir
Q_cost = blkdiag(C1'*Qy1*C1, C2'*Qy2*C2, C3'*Qy3*C3) + 0.001 * eye(12);

% R_cost matrisi
R_cost = 0.01 * eye(6);
%% 1. PERFORMANCE TARGETS
% Continuous-Time Targets
alpha_target = 2;   % Decay rate: Re(s) < -1
zeta_target  = 0.707; % Damping: zeta > 0.707 (45-degree sector)

% Discrete-Time Targets
r_target     = 0.88;  % Radius: |z| < 0.90
angle_target = 20;    % Angle: Max 20 degrees from real axis
structures = {ContStruc_c, ContStruc_dec, ContStruc_string, ContStruc_3ten1e};
struct_names = {'Centralized', 'Decentralized', 'String', '31'};
%% 2. CONTINUOUS-TIME STRUCTURED CONTROL
fprintf('\n==================================================\n');
fprintf('     CONTINUOUS-TIME REGIONAL CONTROL ANALYSIS\n');
fprintf('==================================================\n');


K_ct_results = cell(1,4);

for i = 1:4
    fprintf('Structure: %s...\n ', struct_names{i});
    % [K, rho, feas] = solve_continuous_lmi_comprehensive(A, Bdec, Cdec, N, structures{i}, alpha_target, zeta_target, u_min, u_ss, u_max, x0_norm,0.1,10);
    [K, trace_S, feas] = solve_continuous_lmi_H2(A, Bdec, Cdec, N, structures{i}, Q_cost, R_cost, alpha_target, zeta_target, u_min, u_ss, u_max, x0_norm);
    % if feas ~= 0 && feas ~= 4
    %     [K, rho, feas] = solve_continuous_lmi_distributed(A, Bdec, Cdec, N, structures{i});
    % end
    K_ct_results{i} = K;
    if ~isempty(K)
        plot_eigenvalues(K, 'CT', struct_names{i}, A, B, alpha_target,zeta_target);
        plot_trajectory_ct(A+B*K,B,K, 0.1, x0_norm*x0, struct_names{i})
    end
end

%% 3. DISCRETE-TIME STRUCTURED CONTROL
fprintf('\n==================================================\n');
fprintf('      DISCRETE-TIME REGIONAL CONTROL ANALYSIS\n');
fprintf('==================================================\n');

K_dt_results = cell(1,4);

for i = 1:4
    fprintf('\nStructure: %s', struct_names{i});
    % Using the Circle + Sector Region for DT
    [K, rho, feas] = solve_discrete_lmi_comprehensive(F, Gdec, Hdec, N, structures{i}, 0, r_target,u_min,u_ss,u_max,x0_norm);
    % if feas == 0 || feas == 4
    %     K_dt_results{i} = K;
    % else
        % fprintf('  -> Region too strict. Attempting basic stabilization...\n');
        % [K, rho, feas] = solve_discrete_lmi_distributed(F, Gdec, Hdec, N, structures{i});
        K_dt_results{i} = K;
    % end
        if ~isempty(K)
        plot_eigenvalues(K, 'DT', struct_names{i}, F, G, r_target,angle_target);
        plot_trajectory_dt(F+G*K,G,K, 0.1, 1*x0, struct_names{i})
    end
end

%% 4. SUMMARY TABLE GENERATION
fprintf('\n\n================== FINAL PERFORMANCE SUMMARY ==================\n');
fprintf('%-15s | %-12s | %-12s\n', 'STRUCTURE', 'CT RHO (ABS)', 'DT RHO (RAD)');
fprintf('---------------------------------------------------------------\n');
for i = 1:4
    if ~isempty(K_ct_results{i})
        ct_rho = max(real(eig(A + B*K_ct_results{i})));
    else
        ct_rho = NaN;
    end
    
    if ~isempty(K_dt_results{i})
        dt_rho = max(abs(eig(F + G*K_dt_results{i})));
    else
        dt_rho = NaN;
    end
    
    fprintf('%-15s | %-12.4f | %-12.4f\n', struct_names{i}, ct_rho, dt_rho);
end
fprintf('===============================================================\n');

%% 5. DATA PREPARATION FOR PLOTTING
% We will select the "Star" distributed structure for the simulation
K_final_ct = K_ct_results{4};
K_final_dt = K_dt_results{4};

disp('Control gains computed. Ready for Closed-Loop Simulation.');

%%
%% ============================================================
% 6. CLOSED-LOOP ANALYSIS & TRAJECTORY SIMULATION
% ============================================================

fprintf('\n==================================================\n');
fprintf('      CLOSED-LOOP ANALYSIS & TRAJECTORIES\n');
fprintf('==================================================\n');

rng(1);                 % Tekrarlanabilirlik için
x0 = 0.1*randn(n,1);    % Ortak rastgele başlangıç koşulu

fprintf('Initial condition norm: %.4f\n', norm(x0));

%% -------- CONTINUOUS-TIME CLOSED LOOP --------
if ~isempty(K_final_ct)
    fprintf('\n--- Continuous-Time Closed Loop (STAR structure) ---\n');

    A_cl = A + B*K_final_ct;

    eig_ct_cl = eig(A_cl);
    rho_ct_cl = max(real(eig_ct_cl));

    fprintf('Closed-loop CT eigenvalues:\n');
    disp(eig_ct_cl.');
    fprintf('Spectral abscissa (max Re): %.4f\n', rho_ct_cl);

    if rho_ct_cl < 0
        disp('=> CT closed-loop is ASYMPTOTICALLY STABLE.');
    else
        disp('=> CT closed-loop is NOT asymptotically stable.');
    end

    % Trajectory plot
    plot_trajectory_ct(A_cl,B,K, h, x0, ...
        'CT Closed-Loop Trajectories (Star Structure)');
else
    warning('No CT controller available for Star structure.');
end


%% -------- DISCRETE-TIME CLOSED LOOP --------
if ~isempty(K_final_dt)
    fprintf('\n--- Discrete-Time Closed Loop (STAR structure) ---\n');

    F_cl = F + G*K_final_dt;

    eig_dt_cl = eig(F_cl);
    rho_dt_cl = max(abs(eig_dt_cl));

    fprintf('Closed-loop DT eigenvalues:\n');
    disp(eig_dt_cl.');
    fprintf('Spectral radius (max |z|): %.4f\n', rho_dt_cl);

    if rho_dt_cl < 1
        disp('=> DT closed-loop is ASYMPTOTICALLY STABLE.');
    else
        disp('=> DT closed-loop is NOT asymptotically stable.');
    end

    % Trajectory plot
    plot_trajectory_dt(F_cl,G,K, h, x0, ...
        'DT Closed-Loop Trajectories (Star Structure)');
else
    warning('No DT controller available for Star structure.');
end
