clear variables;
close all;
clc;

% 1. Load Survival Data from Virtual Clinical Trials (N = 10,000)
% Load No Treatment (NT) baseline data
load('No_treatment_Tsurv-2026-8-26-11-7_N=10000.mat');
Tnt = Tsurv';
clearvars -except Tnt;

% Load Radiotherapy (RT) Only data
load('Radio_only_Tsurv-2026-8-26-10-25_N=10000.mat')
Tradio = Tsurv';
clearvars -except Tnt Tradio;

% Load Temozolomide (TMZ) Only data
load('TMZ_only_Tsurv-2026-8-26-10-28_N=10000.mat');
Ttmz = Tsurv';
clearvars -except Tnt Tradio Ttmz;

% Load Stupp Protocol (RT + TMZ) data
load('Stupp_Tsurv-2026-8-26-10-35_N=10000.mat');
Tstupp = Tsurv';
clearvars -except Tnt Tradio Ttmz Tstupp;

n_patients = 10000;

% 2. Prepare Data for MatSurv Analysis
EventVar = cell(n_patients * 2, 1);
EventVar(:) = {'DECEASED'}; % All patients reach the lethal threshold

% Define treatment groups for each comparison
TreatVar_radio = cell(n_patients * 2, 1);
TreatVar_radio(1:n_patients) = {'NT'};
TreatVar_radio(n_patients+1 : 2*n_patients) = {'RT'};

TreatVar_tmz = cell(n_patients * 2, 1);
TreatVar_tmz(1:n_patients) = {'NT'};
TreatVar_tmz(n_patients+1 : 2*n_patients) = {'TMZ'};

TreatVar_stupp = cell(n_patients * 2, 1);
TreatVar_stupp(1:n_patients) = {'NT'};
TreatVar_stupp(n_patients+1 : 2*n_patients) = {'STUPP'};

% Aggregate survival times
times_radio = [Tnt; Tradio];
times_tmz = [Tnt; Ttmz];
times_stupp = [Tnt; Tstupp];

% Define custom x-axis ticks
custom_xticks = [0, 500, 1000, 1500, 2000]; 

% 3. Plotting Kaplan-Meier Curves using MatSurv
% ORDER: {'Treatment', 'NT'} ensures HR < 1 (protective effect)

% --- 3.1. RT vs NT Curve ---
[p_radio, fh_radio, stats_radio] = MatSurv(times_radio, EventVar, TreatVar_radio, ...
    'TimeMax', 2000, 'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'RT', 'NT'});

all_axes_radio = findall(fh_radio, 'Type', 'axes'); 
for i = 1:length(all_axes_radio)
    xticks(all_axes_radio(i), custom_xticks);
end
fh_radio.WindowState = 'maximized';
drawnow;
saveas(fh_radio, 'KM_NT_RT.eps', 'epsc');

% --- 3.2. TMZ vs NT Curve ---
[p_tmz, fh_tmz, stats_tmz] = MatSurv(times_tmz, EventVar, TreatVar_tmz, ...
    'TimeMax', 2000, 'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'TMZ', 'NT'});

all_axes_tmz = findall(fh_tmz, 'Type', 'axes');
for i = 1:length(all_axes_tmz)
    xticks(all_axes_tmz(i), custom_xticks);
end
fh_tmz.WindowState = 'maximized';
drawnow;
saveas(fh_tmz, 'KM_NT_TMZ.eps', 'epsc');

% --- 3.3. Stupp vs NT Curve ---
[p_stupp, fh_stupp, stats_stupp] = MatSurv(times_stupp, EventVar, TreatVar_stupp, ...
    'TimeMax', 2000, 'DispP', true, 'DispHR', true, 'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'STUPP', 'NT'});

all_axes_stupp = findall(fh_stupp, 'Type', 'axes');
for i = 1:length(all_axes_stupp)
    xticks(all_axes_stupp(i), custom_xticks);
end
fh_stupp.WindowState = 'maximized';
drawnow;
saveas(fh_stupp, 'KM_NT_Stupp.eps', 'epsc');