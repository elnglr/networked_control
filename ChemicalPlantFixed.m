%% CHEMICAL PLANT PROJECT – STEP 1 (Decomposition)
clear all; clc; close all;

% Load plant model
MAT08ChemicalPlant;   % contains A,B,C  (from your uploaded file)

n = size(A,1);
m = size(B,2);

%% Define subsystem indices
idx1 = 1:4;
idx2 = 5:8;
idx3 = 9:12;

% Store decomposition
N = 3;  % three reactors

Bdec = cell(N,1);
Cdec = cell(N,1);

% Subsystem 1
Bdec{1} = B(:,1:2);
Cdec{1} = C(idx1,:);

% Subsystem 2
Bdec{2} = B(:,3:4);
Cdec{2} = C(idx2,:);

% Subsystem 3
Bdec{3} = B(:,5:6);
Cdec{3} = C(idx3,:);

disp('Decomposition completed.');

%% STEP 2 — CT and DT Models

% Continuous-time eigenvalues
eigA = eig(A);
rhoA = max(real(eigA));

disp('Continuous-time eigenvalues:');
disp(eigA);
disp(['Spectral abscissa = ', num2str(rhoA)]);

if rhoA < 0
    disp('System is CT asymptotically stable.');
else
    disp('System is NOT CT asymptotically stable.');
end

%% Discretization
h = 0.1;
F = expm(A*h);
G = A \ (F - eye(n)) * B;   % same formula used in Lab 1
H = C;


% Store decomposition
N = 3;  % three reactors

Gdec = cell(N,1);
Hdec = cell(N,1);

% Subsystem 1
Gdec{1} = G(:,1:2);
Hdec{1} = H(idx1,:);

% Subsystem 2
Gdec{2} = G(:,3:4);
Hdec{2} = H(idx2,:);

% Subsystem 3
Gdec{3} = G(:,5:6);
Hdec{3} = H(idx3,:);

% Gdec = cell(N,1);
% Hdec = cell(N,1);
% 
% Gdec{1} = G(:,1:2);    Hdec{1} = Cdec{1};
% Gdec{2} = G(:,3:4);    Hdec{2} = Cdec{2};
% Gdec{3} = G(:,5:6);    Hdec{3} = Cdec{3};

eigF = eig(F);
rhoF = max(abs(eigF));

disp('Discrete-time eigenvalues:');
disp(eigF);
disp(['Spectral radius = ', num2str(rhoF)]);

if rhoF < 1
    disp('System is DT asymptotically stable.');
else
    disp('System is NOT DT asymptotically stable.');
end
%% STEP 3 — FIXED MODES (CT & DT) FOR CHEMICAL PLANT
disp('============================');
disp('      FIXED MODES TEST      ');
disp('============================');

rounding_n = 3;   % rounding precision, same as professor

%% CENTRALIZED STRUCTURE
ContStruc_c = ones(N,N);

FM_c_CT  = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_c,  rounding_n);
FM_c_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_c,  rounding_n);

disp('Centralized CT fixed modes:');
disp(FM_c_CT);
disp('Centralized DT fixed modes:');
disp(FM_c_DT);

%% DECENTRALIZED STRUCTURE
ContStruc_d = diag(ones(N,1));

FM_d_CT  = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_d,  rounding_n);
FM_d_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_d,  rounding_n);

disp('Decentralized CT fixed modes:');
disp(FM_d_CT);
disp('Decentralized DT fixed modes:');
disp(FM_d_DT);

%% DISTRIBUTED STRING STRUCTURE  (1 <-> 2 <-> 3)
ContStruc_string = [1 1 0;
                    1 1 1;
                    0 1 1];

FM_s_CT  = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_string,  rounding_n);
FM_s_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_string,  rounding_n);

disp('Distributed (string) CT fixed modes:');
disp(FM_s_CT);
disp('Distributed (string) DT fixed modes:');
disp(FM_s_DT);

%% DISTRIBUTED STAR STRUCTURE  (center = subsystem 2)
ContStruc_star = [1 1 0;
                  1 1 1;
                  0 1 1]; 

FM_star_CT = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_star, rounding_n);
FM_star_DT = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_star, rounding_n);

disp('Distributed (star) CT fixed modes:');
disp(FM_star_CT);
disp('Distributed (star) DT fixed modes:');
disp(FM_star_DT);
%% STEP 4 – LMI CONTROLLERS FOR CHEMICAL PLANT
% (run after MAT08ChemicalPlant + decomposition + discretization)

clear all; clc; close all;

%% Load model
MAT08ChemicalPlant;    % A,B,C (C = eye(12))
n = size(A,1);  m = size(B,2);

%% Subsystem indices
idx1 = 1:4;
idx2 = 5:8;
idx3 = 9:12;
N = 3;                 % number of subsystems

%% Decomposition (same as before)
Bdec = cell(N,1);
Cdec = cell(N,1);

Bdec{1} = B(:,1:2);    Cdec{1} = C(idx1,:);
Bdec{2} = B(:,3:4);    Cdec{2} = C(idx2,:);
Bdec{3} = B(:,5:6);    Cdec{3} = C(idx3,:);

%% Discretization
h = 0.1;
F = expm(A*h);
G = A \ (F - eye(n)) * B;
H = C;

Gdec = cell(N,1);
Hdec = cell(N,1);

Gdec{1} = G(:,1:2);    Hdec{1} = Cdec{1};
Gdec{2} = G(:,3:4);    Hdec{2} = Cdec{2};
Gdec{3} = G(:,5:6);    Hdec{3} = Cdec{3};

%% Build "total" B and G (needed for simulations)
Btot = [];
Gtot = [];
for i = 1:N
    Btot = [Btot, Bdec{i}];
    Gtot = [Gtot, Gdec{i}];
end

%% ====== CONTROL STRUCTURES ======
ContStruc_c      = ones(N,N);          % centralized
ContStruc_dec    = diag(ones(N,1));    % decentralized
ContStruc_string = [1 1 0;             % distributed string
                    1 1 1;
                    0 1 1];
ContStruc_star   = [1 1 0;             % distributed star (hub = 2)
                    1 1 1;
                    0 1 1];

%% ====== CONTINUOUS-TIME LMI CONTROLLERS ======
disp('===== CONTINUOUS-TIME CONTROLLERS =====');

[K_c_CT,   rho_c_CT,   feas_c_CT]   = LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_c);
[K_dec_CT, rho_dec_CT, feas_dec_CT] = LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_dec);
[K_str_CT, rho_str_CT, feas_str_CT] = LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_string);
[K_star_CT,rho_star_CT,feas_star_CT]= LMI_CT_DeDicont(A,Bdec,Cdec,N,ContStruc_star);

disp(['Centralized:    feas = ',num2str(feas_c_CT),   ',  rho = ',num2str(rho_c_CT)]);
disp(['Decentralized:  feas = ',num2str(feas_dec_CT), ',  rho = ',num2str(rho_dec_CT)]);
disp(['String distr.:  feas = ',num2str(feas_str_CT), ',  rho = ',num2str(rho_str_CT)]);
disp(['Star distr.:    feas = ',num2str(feas_star_CT),',  rho = ',num2str(rho_star_CT)]);

%% ====== DISCRETE-TIME LMI CONTROLLERS ======
disp('===== DISCRETE-TIME CONTROLLERS =====');

[K_c_DT,   rho_c_DT,   feas_c_DT]   = LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_c);
[K_dec_DT, rho_dec_DT, feas_dec_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_dec);
[K_str_DT, rho_str_DT, feas_str_DT] = LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_string);
[K_star_DT,rho_star_DT,feas_star_DT]= LMI_DT_DeDicont(F,Gdec,Hdec,N,ContStruc_star);

disp(['Centralized:    feas = ',num2str(feas_c_DT),   ',  rho = ',num2str(rho_c_DT)]);
disp(['Decentralized:  feas = ',num2str(feas_dec_DT), ',  rho = ',num2str(rho_dec_DT)]);
disp(['String distr.:  feas = ',num2str(feas_str_DT), ',  rho = ',num2str(rho_str_DT)]);
disp(['Star distr.:    feas = ',num2str(feas_star_DT),',  rho = ',num2str(rho_star_DT)]);


%% STEP 5 – CLOSED-LOOP SIMULATIONS

Tfinal = 10;          % seconds
T      = 0:0.01:Tfinal;
steps  = 1:round(Tfinal/h);

% random initial condition
rng(1);               % to make it repeatable
x0 = randn(n,1);

% Preallocate
x_c_CT    = zeros(n,length(T));
x_dec_CT  = x_c_CT;
x_str_CT  = x_c_CT;
x_star_CT = x_c_CT;

x_c_DT    = zeros(n,length(steps));
x_dec_DT  = x_c_DT;
x_str_DT  = x_c_DT;
x_star_DT = x_c_DT;

% Closed-loop matrices
Acl_c_CT    = A + Btot*K_c_CT;
Acl_dec_CT  = A + Btot*K_dec_CT;
Acl_str_CT  = A + Btot*K_str_CT;
Acl_star_CT = A + Btot*K_star_CT;

Fcl_c_DT    = F + Gtot*K_c_DT;
Fcl_dec_DT  = F + Gtot*K_dec_DT;
Fcl_str_DT  = F + Gtot*K_str_DT;
Fcl_star_DT = F + Gtot*K_star_DT;

% Continuous-time trajectories
for k = 1:length(T)
    t = T(k);
    x_c_CT(:,k)    = expm(Acl_c_CT   * t)*x0;
    x_dec_CT(:,k)  = expm(Acl_dec_CT * t)*x0;
    x_str_CT(:,k)  = expm(Acl_str_CT * t)*x0;
    x_star_CT(:,k) = expm(Acl_star_CT* t)*x0;
end

% Discrete-time trajectories
for k = 1:length(steps)
    x_c_DT(:,k)    = (Fcl_c_DT   ^k)*x0;
    x_dec_DT(:,k)  = (Fcl_dec_DT ^k)*x0;
    x_str_DT(:,k)  = (Fcl_str_DT ^k)*x0;
    x_star_DT(:,k) = (Fcl_star_DT^k)*x0;
end

%% ====== PLOTS: H1, H2, H3 ======
figure;
H_idx = [1 5 9];   % indices of H1, H2, H3

for i = 1:3
    % Continuous-time
    subplot(3,2,2*(i-1)+1);
    hold on; grid on;
    title(['H_{',num2str(i),'} – Continuous-time']);
    plot(T, x_c_CT(H_idx(i),:),   'k','LineWidth',1.2);
    plot(T, x_dec_CT(H_idx(i),:), 'm','LineWidth',1.2);
    plot(T, x_str_CT(H_idx(i),:), 'b','LineWidth',1.2);
    plot(T, x_star_CT(H_idx(i),:),'r','LineWidth',1.2);
    xlabel('t [s]');
    ylabel(['H_',num2str(i)]);
    
    % Discrete-time
    subplot(3,2,2*i);
    hold on; grid on;
    title(['H_{',num2str(i),'} – Discrete-time']);
    tk = steps*h;
    plot(tk, x_c_DT(H_idx(i),:),   'k.-');
    plot(tk, x_dec_DT(H_idx(i),:), 'm.-');
    plot(tk, x_str_DT(H_idx(i),:), 'b.-');
    plot(tk, x_star_DT(H_idx(i),:),'r.-');
    xlabel('t [s]');
    ylabel(['H_',num2str(i)]);
end

legend('Centralized','Decentralized','Distributed (string)','Distributed (star)');

