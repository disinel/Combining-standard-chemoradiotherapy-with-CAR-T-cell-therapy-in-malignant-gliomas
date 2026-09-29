% =========================================================================
% Script para calcular la Selección Óptima de Protocolos y 
% los Perfiles Biológicos de las Subpoblaciones (Tablas LaTeX)
% =========================================================================
clear variables; close all; clc;

% =========================================================================
% 1. DEFINICIÓN DE ARCHIVOS
% =========================================================================
file_params = 'Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat';
file_stupp  = 'Stupp_Tsurv-2026-9-28-18-1_N=10000.mat'; 
file_P1     = '3C_Stupp_Tsurv-2026-9-28-17-42_N=10000.mat';       
file_P2     = 'Stupp_3C_Tsurv-2026-9-28-17-50_N=10000.mat';        
file_P3     = '3C_Stupp_3C_Tsurv-2026-9-28-17-22_N=10000.mat';     
file_P4     = 'Stupp(C)_3C_Stupp(A)_Tsurv-2026-9-28-18-22_N=10000.mat'; 
file_P5     = 'Stupp(C)_3C_Stupp(A)_3C_Tsurv-2026-9-28-18-10_N=10000.mat'; 

latex_names = {'Stupp', '3C--Stupp', 'Stupp--3C', '3C--Stupp--3C', ...
               'Stupp(C)--3C--Stupp(A)', 'Stupp(C)--3C--Stupp(A)--3C'};
M = length(latex_names);

% =========================================================================
% 2. CARGA DE DATOS Y MANEJO DE NaNs
% =========================================================================
fprintf('Cargando parámetros de la población virtual...\n');
load(file_params, 'r1val', 'tauval', 'Kival', 'delta1val', 'delta2val', ...
                  'rho1val', 'rho2val', 'rho4val', 'N');

fprintf('Cargando tiempos de supervivencia...\n');
data_stupp = load(file_stupp, 'Tsurv'); Tsurv_Stupp = data_stupp.Tsurv(:);
data_P1    = load(file_P1, 'Tsurv');    Tsurv_P1    = data_P1.Tsurv(:);
data_P2    = load(file_P2, 'Tsurv');    Tsurv_P2    = data_P2.Tsurv(:);
data_P3    = load(file_P3, 'Tsurv');    Tsurv_P3    = data_P3.Tsurv(:);
data_P4    = load(file_P4, 'Tsurv');    Tsurv_P4    = data_P4.Tsurv(:);
data_P5    = load(file_P5, 'Tsurv');    Tsurv_P5    = data_P5.Tsurv(:);

Tfinal = 150000;
Tsurv_Stupp(isnan(Tsurv_Stupp)) = Tfinal;
Tsurv_P1(isnan(Tsurv_P1)) = Tfinal;
Tsurv_P2(isnan(Tsurv_P2)) = Tfinal;
Tsurv_P3(isnan(Tsurv_P3)) = Tfinal;
Tsurv_P4(isnan(Tsurv_P4)) = Tfinal;
Tsurv_P5(isnan(Tsurv_P5)) = Tfinal;

% Matriz total: 10000 pacientes x 6 protocolos
All_Tsurv = [Tsurv_Stupp, Tsurv_P1, Tsurv_P2, Tsurv_P3, Tsurv_P4, Tsurv_P5];

% =========================================================================
% 3. CÁLCULO DE ÓPTIMOS (T_s^*) Y ESTADÍSTICAS
% =========================================================================
% Encontrar el tiempo máximo y el índice del protocolo ganador para cada VP
[Ts_star, best_idx] = max(All_Tsurv, [], 2);
Ts_star_months = Ts_star / 30;

% Mediana global de Ts* y su 95% CI
overall_med = median(Ts_star_months);
sort_T = sort(Ts_star_months);
idx_low  = max(1, floor(N/2 - 1.96 * sqrt(N)/2));
idx_high = min(N, ceil(N/2 + 1.96 * sqrt(N)/2));
ci_low = sort_T(idx_low);
ci_up  = sort_T(idx_high);

% =========================================================================
% 4. GENERACIÓN DE TABLA 1: SELECCIÓN DE PROTOCOLOS
% =========================================================================
fprintf('\n%% =========================================================================\n');
fprintf('%% TABLA 1: PROTOCOL SELECTION (LaTeX)\n');
fprintf('%% =========================================================================\n\n');

fprintf('\\begin{table*}[t!]\n');
fprintf('    \\caption{Distribution of patient-specific optimal protocol selection across the virtual population ($N=%d$). The \\textit{Pop. Median OS} column reports the median OS obtained when each protocol is applied uniformly to the entire cohort, whereas the \\textit{Subpop. Median OS} column reports the median survival within the subset of VPs for whom that protocol is individually optimal. Furthermore, the observed median and 95\\%% CI for $T_s^\\ast$ is $%.2f$ months ($95$\\%% CI: $%.2f-%.2f$ months).}\n', N, overall_med, ci_low, ci_up);
fprintf('    \\label{tab:protocol_selection}\n');
fprintf('    \\centering\n');
fprintf('    \\small\n');
fprintf('    \\begin{ruledtabular}\n');
fprintf('        \\begin{tabular}{lccc}\n');
fprintf('            \\textbf{Protocol ($P$)} & \\textbf{Optimal for VPs (\\%%)} & \\textbf{Pop. Median OS (months)} & \\textbf{Subpop. Median OS ($T_s^\\ast$,months)} \\\\\n');
fprintf('            \\hline\n');

for k = 1:M
    mask = (best_idx == k);
    pct_opt = (sum(mask) / N) * 100;
    pop_med = median(All_Tsurv(:, k)) / 30;
    
    if sum(mask) > 0
        sub_med = median(All_Tsurv(mask, k)) / 30;
        sub_med_str = sprintf('%.2f', sub_med);
    else
        sub_med_str = '--';
    end
    
    % Aplicar negrita si es el protocolo más elegido (> 30% como en tu ejemplo)
    if pct_opt == max(arrayfun(@(x) sum(best_idx==x)/N*100, 1:M))
        fprintf('            %-30s & \\textbf{%.2f} & \\textbf{%.2f} & %s \\\\\n', latex_names{k}, pct_opt, pop_med, sub_med_str);
    else
        fprintf('            %-30s & %.2f          & %.2f          & %s \\\\\n', latex_names{k}, pct_opt, pop_med, sub_med_str);
    end
end

fprintf('        \\end{tabular}\n');
fprintf('    \\end{ruledtabular}\n');
fprintf('\\end{table*}\n');

% =========================================================================
% 5. GENERACIÓN DE TABLA 2: PERFILES BIOLÓGICOS DE LAS SUBPOBLACIONES
% =========================================================================
fprintf('\n%% =========================================================================\n');
fprintf('%% TABLA 2: BIOLOGICAL PROFILES (LaTeX)\n');
fprintf('%% =========================================================================\n\n');

fprintf('\\begin{table*}[t!]\n');
fprintf('    \\caption{Median parameter values across the global virtual population and within the optimal patient subpopulations for each clinical protocol.}\n');
fprintf('    \\label{tab:biological_profiles_subpopulations}\n');
fprintf('    \\centering\n');
fprintf('    \\small\n');
fprintf('    \\begin{ruledtabular}\n');
fprintf('        \\begin{tabular}{lcccccccc}\n');
fprintf('            \\textbf{Optimal Subpopulation} & \\textbf{$r_1$} & \\textbf{$\\tau$} & \\textbf{Ki-67} & \\textbf{$\\delta_1$} & \\textbf{$\\delta_2$} & \\textbf{$\\rho_1$} & \\textbf{$\\rho_2,\\rho_3$} & \\textbf{$\\rho_4$} \\\\\n');
fprintf('            \\hline\n');

% Fila Global
g_r1 = median(r1val); g_tau = median(tauval); g_ki = median(Kival);
g_d1 = median(delta1val); g_d2 = median(delta2val); 
g_rho1 = median(rho1val); g_rho2 = median(rho2val); g_rho4 = median(rho4val);

fprintf('            \\textbf{Global Population ($N=%d$)} & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
    N, g_r1, g_tau, g_ki, g_d1, g_d2, g_rho1, g_rho2, g_rho4);
fprintf('            \\hline\n');

% Filas por Subpoblación (Saltamos Stupp si el porcentaje es 0, como en tu ejemplo original)
for k = 2:M
    mask = (best_idx == k);
    if sum(mask) > 0
        s_r1 = median(r1val(mask)); s_tau = median(tauval(mask)); s_ki = median(Kival(mask));
        s_d1 = median(delta1val(mask)); s_d2 = median(delta2val(mask)); 
        s_rho1 = median(rho1val(mask)); s_rho2 = median(rho2val(mask)); s_rho4 = median(rho4val(mask));
        
        fprintf('            \\textbf{%-28s} & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
            latex_names{k}, s_r1, s_tau, s_ki, s_d1, s_d2, s_rho1, s_rho2, s_rho4);
    end
end

fprintf('        \\end{tabular}\n');
fprintf('    \\end{ruledtabular}\n');
fprintf('\\end{table*}\n');