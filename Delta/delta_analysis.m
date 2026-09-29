clear variables;
close all;
clc;
% =========================================================================
% 1. FILE DEFINITIONS
% =========================================================================
file_params = 'Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat';
file_stupp  = 'Stupp_Tsurv-2026-9-28-18-1_N=10000.mat'; 
file_P1     = '3C_Stupp_Tsurv-2026-9-28-17-42_N=10000.mat';       % 3C-Stupp
file_P2     = 'Stupp_3C_Tsurv-2026-9-28-17-50_N=10000.mat';        % Stupp-3C
file_P3     = '3C_Stupp_3C_Tsurv-2026-9-28-17-22_N=10000.mat';     % 3C-Stupp-3C
file_P4     = 'Stupp(C)_3C_Stupp(A)_Tsurv-2026-9-28-18-22_N=10000.mat'; % StuppC-3C-StuppA
file_P5     = 'Stupp(C)_3C_Stupp(A)_3C_Tsurv-2026-9-28-18-10_N=10000.mat'; % StuppC-3C-StuppA-3C
P_names = {'3C-Stupp', 'Stupp-3C', '3C-Stupp-3C', 'Stupp(C)-3C-Stupp(A)', 'Stupp(C)-3C-Stupp(A)-3C'};
P_names_tex = {'\textrm{3C--Stupp}', '\textrm{Stupp--3C}', '\textrm{3C--Stupp--3C}', ...
               '\textrm{Stupp(C)--3C--Stupp(A)}', '\textrm{Stupp(C)--3C--Stupp(A)--3C}'};
M = length(P_names);
% =========================================================================
% 2. LOAD DATA & HANDLE SURVIVORS
% =========================================================================
fprintf('Loading virtual population parameters...\n');
load(file_params, 'r1val', 'delta1val', 'rho4val', 'N');
fprintf('Loading survival times...\n');
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
All_Tsurv = [Tsurv_P1, Tsurv_P2, Tsurv_P3, Tsurv_P4, Tsurv_P5];
% =========================================================================
% 3. CALCULATE ADDITIVE BENEFIT (DELTA)
% =========================================================================
Delta_matrix = All_Tsurv - Tsurv_Stupp;
Deltas = num2cell(Delta_matrix, 1);
% =========================================================================
% 4. RESPONDER QUANTIFICATION TABLE
% =========================================================================
fprintf('\n========================================================================\n');
fprintf('   RESPONDER QUANTIFICATION: delta_i^P > 60 days AND delta_i^P > 365 days\n');
fprintf('========================================================================\n');
fprintf('%-28s %-18s %-18s\n', 'Protocol', '> 60 days (N, %)', '> 365 days (N, %)');
fprintf('------------------------------------------------------------------------\n');
for i = 1:M
    n_60   = sum(Deltas{i} > 60);
    pct_60 = (n_60 / N) * 100;
    
    n_365   = sum(Deltas{i} > 365);
    pct_365 = (n_365 / N) * 100;
    
    fprintf('%-28s %5d (%5.2f%%)       %5d (%5.2f%%)\n', ...
        P_names{i}, n_60, pct_60, n_365, pct_365);
end
% =========================================================================
% 5. PLOT ADDITIVE BENEFIT DISTRIBUTIONS (Histograms)
% =========================================================================
fprintf('\nGenerating Additive Benefit Distributions...\n');
% Ajustado a formato horizontal extremo (ancho 1500, alto 350)
fig_dist = figure('Name', 'Additive Benefit Distributions', 'Color', 'w', 'Position', [100, 100, 1500, 750]);
for i = 1:M
    % Distribución horizontal en una sola fila (1 fila, 5 columnas)
    subplot(2, 3, i);
    histogram(Deltas{i}, 50, 'FaceColor', '#0072BD');
    set(gca, 'YScale', 'log'); 
    xlabel('$\Delta_i^P \textrm{ (days)}$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('N (Log Scale)', 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
end
exportgraphics(fig_dist, 'Benefit_Distributions.pdf', 'ContentType', 'vector');
fprintf('Saved Benefit_Distributions.pdf\n');
% =========================================================================
% 6. SCATTER PLOTS: ADDITIVE BENEFIT VS PARAMETERS
% =========================================================================
fprintf('Generating Additive Scatter Plots...\n');
y_limits_delta = [-20, 500];
for i = 1:M
    % Ajustado a formato horizontal rectangular (ancho 1200, alto 350)
    fig_scatter = figure('Name', ['Scatter Additive: ', P_names{i}], 'Color', 'w', 'Position', [100, 100, 1200, 350]);
    
    % Panel 1: r1 vs Delta
    subplot(1, 3, 1);
    scatter(r1val, Deltas{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(60, '--r', 'LineWidth', 1.3); hold off;
    ylim(y_limits_delta);
    xlabel('$r_1$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\Delta_i^P \textrm{ (days)}$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    % Panel 2: delta2 (delta1val) vs Delta
    subplot(1, 3, 2);
    scatter(delta1val, Deltas{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(60, '--r', 'LineWidth', 1.3); hold off;
    ylim(y_limits_delta);
    xlabel('$\delta_2$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\Delta_i^P \textrm{ (days)}$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    % Panel 3: rho4 vs Delta
    subplot(1, 3, 3);
    scatter(rho4val, Deltas{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(60, '--r', 'LineWidth', 1.3); hold off;
    ylim(y_limits_delta);
    xlabel('$\rho_4$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\Delta_i^P \textrm{ (days)}$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    safe_name = strrep(strrep(P_names{i}, '(', '_'), ')', '');
    file_name = sprintf('Scatter_plot_%s.pdf', safe_name);
    exportgraphics(fig_scatter, file_name, 'ContentType', 'vector');
    fprintf('Saved %s\n', file_name);
end
fprintf('\nDelta analysis complete. All figures saved.\n');