clear variables;
close all;

% 1. Load Survival Data from Virtual Clinical Trials (N = 10,000)
% Note: Stupp baseline data is loaded into Tnt for comparison formatting
load('Stupp_Tsurv-2026-9-28-18-1_N=10000.mat');
Tnt = Tsurv';
clearvars -except Tnt;

% 1. Stupp(C)-3C-Stupp(A)
load('Stupp(C)_3C_Stupp(A)_Tsurv-2026-9-28-18-22_N=10000.mat');
TStuppC_3C_StuppA = Tsurv';
clearvars -except Tnt TStuppC_3C_StuppA;

% 2. Stupp(C)-3C-Stupp(A)-3C
load('Stupp(C)_3C_Stupp(A)_3C_Tsurv-2026-9-28-18-10_N=10000.mat');
TStuppC_3C_StuppA_3C = Tsurv';
clearvars -except Tnt TStuppC_3C_StuppA TStuppC_3C_StuppA_3C;

% 3. Stupp-3C
load('Stupp_3C_Tsurv-2026-9-28-17-50_N=10000.mat');
TStupp_3C = Tsurv';
clearvars -except Tnt TStuppC_3C_StuppA TStuppC_3C_StuppA_3C TStupp_3C;

% 4. 3C-Stupp
load('3C_Stupp_Tsurv-2026-9-28-17-42_N=10000.mat');
T3C_Stupp = Tsurv';
clearvars -except Tnt TStuppC_3C_StuppA TStuppC_3C_StuppA_3C TStupp_3C T3C_Stupp;

% 5. 3C-Stupp-3C
load('3C_Stupp_3C_Tsurv-2026-9-28-17-22_N=10000.mat');
T3C_Stupp_3C = Tsurv';
clearvars -except Tnt TStuppC_3C_StuppA TStuppC_3C_StuppA_3C TStupp_3C T3C_Stupp T3C_Stupp_3C;

n_patients = 10000;

% 2. Prepare Data for MatSurv Format
EventVar = cell(n_patients * 2, 1);
EventVar(:) = {'DECEASED'}; % All patients reach the lethal threshold

% Prepare independent treatment group vectors for each comparison
TreatVar_1 = cell(n_patients * 2, 1);
TreatVar_1(1:n_patients) = {'Stupp'};
TreatVar_1(n_patients+1 : 2*n_patients) = {'Stupp(C)-3C-Stupp(A)'};

TreatVar_2 = cell(n_patients * 2, 1);
TreatVar_2(1:n_patients) = {'Stupp'};
TreatVar_2(n_patients+1 : 2*n_patients) = {'Stupp(C)-3C-Stupp(A)-3C'};

TreatVar_3 = cell(n_patients * 2, 1);
TreatVar_3(1:n_patients) = {'Stupp'};
TreatVar_3(n_patients+1 : 2*n_patients) = {'Stupp-3C'};

TreatVar_4 = cell(n_patients * 2, 1);
TreatVar_4(1:n_patients) = {'Stupp'};
TreatVar_4(n_patients+1 : 2*n_patients) = {'3C-Stupp'};

TreatVar_5 = cell(n_patients * 2, 1);
TreatVar_5(1:n_patients) = {'Stupp'};
TreatVar_5(n_patients+1 : 2*n_patients) = {'3C-Stupp-3C'};

% Aggregate survival times
times_1 = [Tnt; TStuppC_3C_StuppA];
times_2 = [Tnt; TStuppC_3C_StuppA_3C];
times_3 = [Tnt; TStupp_3C];
times_4 = [Tnt; T3C_Stupp];
times_5 = [Tnt; T3C_Stupp_3C];

% Define custom x-axis ticks
custom_xticks = [0, 500, 1000, 1500, 2000]; 

% 3. Plotting Kaplan-Meier Curves with MatSurv
% NOTE: Order in 'GroupsToUse' adjusted to control Hazard Ratio presentation

% --- 3.1. Stupp(C)-3C-Stupp(A) vs Stupp ---
[p1, fh1, stats1] = MatSurv(times_1, EventVar, TreatVar_1, ...
    'TimeMax', 2000, 'XLim', [0 2000], 'XTicks', custom_xticks, ...
    'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'Stupp(C)-3C-Stupp(A)', 'Stupp'});
fh1.WindowState = 'maximized';
saveas(fh1, 'KM_Stupp_StuppC_3C_StuppA.eps', 'epsc');

% --- 3.2. Stupp(C)-3C-Stupp(A)-3C vs Stupp ---
[p2, fh2, stats2] = MatSurv(times_2, EventVar, TreatVar_2, ...
    'TimeMax', 2000, 'XLim', [0 2000], 'XTicks', custom_xticks, ...
    'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'Stupp(C)-3C-Stupp(A)-3C', 'Stupp'});
fh2.WindowState = 'maximized';
saveas(fh2, 'KM_Stupp_StuppC_3C_StuppA_3C.eps', 'epsc');

% --- 3.3. Stupp-3C vs Stupp ---
[p3, fh3, stats3] = MatSurv(times_3, EventVar, TreatVar_3, ...
    'TimeMax', 2000, 'XLim', [0 2000], 'XTicks', custom_xticks, ...
    'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'Stupp-3C', 'Stupp'});
fh3.WindowState = 'maximized';
saveas(fh3, 'KM_Stupp_Stupp_3C.eps', 'epsc');

% --- 3.4. 3C-Stupp vs Stupp ---
[p4, fh4, stats4] = MatSurv(times_4, EventVar, TreatVar_4, ...
    'TimeMax', 2000, 'XLim', [0 2000], 'XTicks', custom_xticks, ...
    'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'3C-Stupp', 'Stupp'});
fh4.WindowState = 'maximized';
saveas(fh4, 'KM_Stupp_3C_Stupp.eps', 'epsc');

% --- 3.5. 3C-Stupp-3C vs Stupp ---
[p5, fh5, stats5] = MatSurv(times_5, EventVar, TreatVar_5, ...
    'TimeMax', 2000, 'XLim', [0 2000], 'XTicks', custom_xticks, ...
    'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'3C-Stupp-3C', 'Stupp'});
fh5.WindowState = 'maximized';
saveas(fh5, 'KM_Stupp_3C_Stupp_3C.eps', 'epsc');