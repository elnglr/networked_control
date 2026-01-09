%% Project_Part1: 
% System Behaviour, Decomposition, 
% OpenLoop Anaylsis, Eigenvalues and Stability Statements,
% Fixed Modes and the conclusion of the 1st Part
%


% STATE VECTOR x
% x = [
%      H1 xa1 xb1 T1;  ---> Reactor 1 --> Subsys1
%      H2 xa2 xb2 T2;  ---> Reactor 2 --> Subsys2
%      H3 xa3 xb3 T3;  ---> Seperator --> Subsys3
%                    ]'


% H_i   : Liquid level in unit i               [m]
% xa_i  : Mass fraction of component A         [-]
% xb_i  : Mass fraction of component B         [-]
% T_i   : Temperature                          [K]




% Reactor 1 : states 1–4
% Reactor 2 : states 5–8
% Separator : states 9–12


% INPUT VECTOR u
% u = [
%      Ff1 Q1; --> Subsys1
%      Ff2 Q2; --> Subsys2
%      Fr Q3;  --> Subsys3
%              ]'


% Ff1 : Fresh feed flow to Reactor 1           [kg/s]
% Q1  : Heat input to Reactor 1                [kJ/s]
% Ff2 : Fresh feed flow to Reactor 2           [kg/s]
% Q2  : Heat input to Reactor 2                [kJ/s]
% Fr  : Recycle flow from Separator to R1      [kg/s]
% Q3  : Heat input to Separator                [kJ/s]



%%
clc;clear; close all;

% Let's recall the system matrices
MAT08ChemicalPlant;

n = size(A,1); % Num of States
m = size(B,2); % Num of Inputs
h = 0.1; % Sampling Time [s]


% Let's discretize the system
F = expm(A*h); 
G = A \ (expm(A*h) - expm(A*0)) * B; 
H = C;

%% DECOMPOSITION into 3 subsystems acc. to given 08_ChemicalPlant.pdf
% Subsys1 : x1 = [H1, xa1, xb1, T1] , u1 = [Ff1, Q1]
% Subsys2 : x2 = [H2, xa2, xb2, T2] , u2 = [Ff2, Q2]
% Subsys3 : x3 = [H3, xa3, xb3, T3] , u3 = [Fr, Q3]

N = 3; % Num of Subsystem


idx1 = 1:4;  % Subsys1, First 4 States
idx2 = 5:8;  % Subsys2, Second 4 States
idx3 = 9:12; % Subsys3, Third 4 States


% Cont-Time
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



% Discrete-Time
Gdec = cell(N,1);
Hdec = cell(N,1);

% Subsystem 1
Gdec{1} = G(:,1:2);
Hdec{1} = H(idx1,:);

% Subsystem 2
Gdec{2} = G(:,3:5);
Hdec{2} = H(idx2,:);

% Subsystem 3
Gdec{3} = G(:,6);
Hdec{3} = H(idx3,:);



disp('Decomposition completed.');


%% CONTINUOUS-TIME EIGENVALUES

eigA = eig(A);

figure;
plot(real(eigA), imag(eigA), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; hold on;

xline(0,'k-','LineWidth',1);
yline(0,'k-','LineWidth',1);

xlabel('Real Part');
ylabel('Imaginary Part');
title('Open-loop Eigenvalues of A (Continuous Time)');



%% DISCRETE-TIME EIGENVALUES

eigF = eig(F);

figure;
plot(real(eigF), imag(eigF), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on; hold on;

theta = linspace(0,2*pi,500);
plot(cos(theta), sin(theta),'k-','LineWidth',1);

xlabel('Real Part');
ylabel('Imaginary Part');
title('Open-loop Eigenvalues of F (Discrete Time)');
axis equal



%% STABILITY STATEMENTS
rhoA = max(real(eigA));

disp('Continuous-time eigenvalues:');
disp(eigA);
disp(['Spectral abscissa = ', num2str(rhoA)]);

if rhoA < 0
    disp('System is CT asymptotically stable.');
else
    disp('System is NOT CT asymptotically stable.');
end


rhoF = max(abs(eigF));

disp('Discrete-time eigenvalues:');
disp(eigF);
disp(['Spectral radius = ', num2str(rhoF)]);

if rhoF < 1
    disp('System is DT asymptotically stable.');
else
    disp('System is NOT DT asymptotically stable.');
end




%%

Tsim = 50;              
t = linspace(0,Tsim,1000);

x0 = zeros(12,1);
x0(1) = 0.1;   % Different Init. Cond. to one of its state

sys_ol = ss(A, B, eye(12), zeros(12,6));

% Simulate
[y,t,x] = initial(sys_ol, x0, t);



figure; 
subplot(2,2,1)
plot(t,x(:,1),'LineWidth',1.5); grid on
title('H1 (Reactor 1 Level)')

subplot(2,2,2)
plot(t,x(:,5),'LineWidth',1.5); grid on
title('H2 (Reactor 2 Level)')

subplot(2,2,3)
plot(t,x(:,9),'LineWidth',1.5); grid on
title('H3 (Separator Level)')

subplot(2,2,4)
plot(t,x(:,7),'LineWidth',1.5); grid on
title('x_{B3} (Product concentration)')

sgtitle('Open-loop Impulse Response (u = 0)')


%% FIXED MODES
rounding_n=4;

% CENTRALIZED STRUCTURE
ContStruc_c = ones(N,N); % Centralized Structure
FM_c_CT  = di_fixed_modes(A, Bdec, Cdec, N, ContStruc_c,  rounding_n); % Cont-Time
FM_c_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_c,  rounding_n);  % Discrete-Time

disp('Centralized CT fixed modes:');
disp(FM_c_CT);
disp('Centralized DT fixed modes:');
disp(FM_c_DT);




% DECENTRALIZED STRUCTURE
ContStruc_d = diag(ones(N,1)); % Decentralized Structure

FM_d_CT  = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_d,  rounding_n); % Cont-Time
FM_d_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_d,  rounding_n);  % Discrete-Time

disp('Decentralized CT fixed modes:');
disp(FM_d_CT);
disp('Decentralized DT fixed modes:');
disp(FM_d_DT);






% DISTRIBUTED STRUCTURE
ContStruc_string = [1 0 1;
                    1 1 0;
                    0 1 1];


% Here We need to choose the structure, check the paper, may be there will
% some suggestions about the structure of network !!!


FM_s_CT  = di_fixed_modes(A,  Bdec, Cdec, N, ContStruc_string,  rounding_n);
FM_s_DT  = di_fixed_modes(F, Gdec, Hdec, N, ContStruc_string,  rounding_n);

disp('Distributed (string) CT fixed modes:');
disp(FM_s_CT);
disp('Distributed (string) DT fixed modes:');
disp(FM_s_DT);


%%

%% 1. Veri Hazırlığı
%% 1. Veri Hazırlığı ve Ön Analiz
sys_ct = ss(A, B, eye(n), zeros(n,m));
sys_dt = ss(F, G, eye(n), zeros(n,m), h);

%% 1. Dinamik Veri Hazırlığı ve Ts Hesaplama
indices = [1, 5, 9, 7]; % H1, H2, H3, xB3
titles = {'H1 (Reactor 1 Level)', 'H2 (Reactor 2 Level)', ...
          'H3 (Separator Level)', 'x_{B3} (Product Concentration)'};
y_labels = {'Level [m]', 'Level [m]', 'Level [m]', 'Concentration [wt.%]'};

% Simülasyon için başlangıç koşullarını belirle
x0 = zeros(12,1); 
x0(1) = 0.1; % Örnek: Reactor 1 seviye sapması
x0(5) = 0.1; % Örnek: Reactor 2 seviye sapması
x0(9) = 0.1;  % Örnek: Separator seviye sapması



% Ts hesaplamak için önce uzun bir ön simülasyon yap (Örn: 50 saniye)
t_pre = linspace(0, 50, 5000);
[~, ~, x_pre] = initial(sys_ct, x0, t_pre);

Ts_dynamic = zeros(1,4);
for i = 1:4
    sig = x_pre(:, indices(i));
    final_val = sig(end); % Kararlı hal (u=0 olduğu için genelde 0'dır)
    peak_dev = max(abs(sig - final_val));
    
    if peak_dev < 1e-6 % Eğer sinyal zaten hareket etmiyorsa
        Ts_dynamic(i) = 0;
    else
        % %2 kriteri için sondan başa doğru bak
        idx = find(abs(sig - final_val) > 0.02 * peak_dev, 1, 'last');
        if isempty(idx)
            Ts_dynamic(i) = 0;
        else
            Ts_dynamic(i) = t_pre(idx);
        end
    end
end

% Seviye grafikleri (H1, H2, H3) için ortak Y-ölçeği belirle
level_data = x_pre(:, [1, 5, 9]);
y_min_common = min(level_data(:)) * 1.1;
y_max_common = max(level_data(:)) * 1.1;

%% 2. FIGURE 1: Continuous-Time Response (Dinamik Ts)
figure('Name', 'Continuous-Time Response (Dynamic Ts)');
for i = 1:4
    t_end = Ts_dynamic(i) + 2; % Yerleşme süresi + 2 saniye pay
    if t_end < 2, t_end = 5; end % Çok hızlı sistemler için alt sınır
    
    t_plot = linspace(0, t_end, 1000);
    [~, ~, x_temp] = initial(sys_ct, x0, t_plot);
    
    subplot(2, 2, i);
    plot(t_plot, x_temp(:, indices(i)), 'b-', 'LineWidth', 1.5);
    title(sprintf('%s\nTs (2%%) = %.2fs', titles{i}, Ts_dynamic(i)));
    xlabel('Time (s)'); ylabel(y_labels{i});
    xlim([0 t_end]);
    
    if i < 4, ylim([y_min_common, y_max_common]); end
    grid on;
end
sgtitle('Continuous-Time Open-loop Impulse Response (Dynamic Scaling)');

%% 3. FIGURE 2: Discrete-Time Response (Dinamik Ts)
figure('Name', 'Discrete-Time Response (Dynamic Ts)');
for i = 1:4
    t_end = Ts_dynamic(i) + 2;
    if t_end < 2, t_end = 5; end
    
    t_disc = 0:h:t_end;
    [~, ~, x_temp_dt] = initial(sys_dt, x0, t_disc);
    
    subplot(2, 2, i);
    stairs(t_disc, x_temp_dt(:, indices(i)), 'r', 'LineWidth', 1.5);
    title(sprintf('%s\nTs (2%%) = %.2fs', titles{i}, Ts_dynamic(i)));
    xlabel('Time (s)'); ylabel(y_labels{i});
    xlim([0 t_end]);
    
    if i < 4, ylim([y_min_common, y_max_common]); end
    grid on;
end
sgtitle(sprintf('Discrete-Time Open-loop Impulse Response (h = %.1f s)', h));