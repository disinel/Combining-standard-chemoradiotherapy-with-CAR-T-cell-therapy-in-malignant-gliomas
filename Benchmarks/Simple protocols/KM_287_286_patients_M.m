clear variables;
close all;

% 1. Load Survival Data from Virtual Clinical Trials
% Load Radiotherapy (RT) Only data
load('Radio_only_Tsurv-2026-8-26-10-23_N=286.mat');
Tro = Tsurv'; 
Nr = N;
clearvars -except Tro Nr;

% Load Stupp Protocol (RT + TMZ) data
load('Stupp_Tsurv-2026-8-26-10-32_N=287.mat');
Trtmz = Tsurv'; 
Ntmz = N;
clearvars -except Tro Trtmz Nr Ntmz;

% 2. Prepare Data for MatSurv Analysis
n_patients = Nr + Ntmz;

EventVar = cell(n_patients, 1);
EventVar(:) = {'DECEASED'}; % All patients reach the lethal threshold

TreatVar = cell(n_patients, 1);
TreatVar(1:Nr) = {'RT'};
TreatVar(Nr+1:n_patients) = {'Stupp'};

times = [Tro; Trtmz];

% Define custom x-axis ticks (months converted to days: 0, 6, 12, 18, 24, 28)
custom_xticks = 30 .* [0, 6, 12, 18, 24, 28]; 

% 3. Plot Kaplan-Meier Curves using MatSurv
[p, fh, stats] = MatSurv(times, EventVar, TreatVar, ...
    'TimeMax', 28*30, ...
    'XLim', [0, 28*30], ...     
    'DispP', true, ...
    'DispHR', true, ...
    'Xlabel', 'Time (days)', ...
    'GroupsToUse', {'Stupp', 'RT'});

% Apply custom ticks and limits to all generated axes (plot and risk table)
all_axes = findall(fh, 'Type', 'axes'); 
for i = 1:length(all_axes)
    xticks(all_axes(i), custom_xticks);
    xlim(all_axes(i), [0, 28*30]); % Ensure secondary axes respect the limits
end

% 4. Save Figure
fh.WindowState = 'maximized';
drawnow;
saveas(fh, 'KM_RT_TMZplusRT_286_287.eps', 'epsc');