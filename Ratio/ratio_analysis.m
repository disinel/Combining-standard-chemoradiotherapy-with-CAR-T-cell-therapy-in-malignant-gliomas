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
% 3. CALCULATE MULTIPLICATIVE BENEFIT (RATIO)
% =========================================================================
Ratio_matrix = All_Tsurv ./ Tsurv_Stupp;
Ratios = num2cell(Ratio_matrix, 1);
% =========================================================================
% 4. PLOT MULTIPLICATIVE BENEFIT DISTRIBUTIONS (Histograms)
% =========================================================================
fprintf('\nGenerating Multiplicative Benefit Distributions...\n');
% Formato horizontal: ancho 1500, alto 350
fig_dist = figure('Name', 'Multiplicative Benefit Distributions', 'Color', 'w', 'Position', [100, 100, 1500, 750]);
for i = 1:M
    % Distribución en una sola fila
    subplot(2, 3, i);
    histogram(Ratio_matrix(:, i), 50, 'FaceColor', '#D95319');
    set(gca, 'YScale', 'log'); 
    xlabel('$R_i^P = T_{s,i}^P / T_{s,i}^{\textrm{Stupp}}$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('N (Log Scale)', 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    ylim([0, 10^4]);
    xlim([0.9, max(Ratio_matrix(:, i)) + 0.1]);
end
exportgraphics(fig_dist, 'Multiplicative_Benefit_Distributions.pdf', 'ContentType', 'vector');
fprintf('Saved Multiplicative_Benefit_Distributions.pdf\n');
% =========================================================================
% 5. CORRELATIONS & SCATTER PLOTS: MULTIPLICATIVE BENEFIT
% =========================================================================
fprintf('\n========================================================================\n');
fprintf('   MULTIPLICATIVE BENEFIT CORRELATIONS: T_P / T_{Stupp}\n');
fprintf('========================================================================\n');
params_matrix = [r1val(:), delta1val(:), rho4val(:)];
param_tex     = {'r_1', '\delta_2', '\rho_4'};
y_limits_ratio = [0.9, 1.3]; 
for i = 1:M
    fprintf('\n--- Protocol: %s ---\n', P_names{i});
    for j = 1:3
        [rho_val, p_val] = corr(params_matrix(:,j), Ratios{i}(:), 'Type', 'Spearman');
        fprintf('Parameter %-10s : rho = %8.4f (p = %.2e)\n', param_tex{j}, rho_val, p_val);
    end
    
    % Formato horizontal rectangular: ancho 1200, alto 350
    fig_ratio = figure('Name', ['Multiplicative Benefit: ', P_names{i}], 'Color', 'w', 'Position', [100, 100, 1200, 350]);
    
    % Panel 1: r1 vs Ratio
    subplot(1, 3, 1);
    scatter(r1val, Ratios{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(1.0, '--k', 'LineWidth', 1.0); hold off;
    ylim(y_limits_ratio);
    xlabel('$r_1$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(['$T_{s,i}^{', P_names_tex{i}, '} / T_{s,i}^{\textrm{Stupp}}$'], 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    % Panel 2: delta2 vs Ratio
    subplot(1, 3, 2);
    scatter(delta1val, Ratios{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(1.0, '--k', 'LineWidth', 1.0); hold off;
    ylim(y_limits_ratio);
    xlabel('$\delta_2$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(['$T_{s,i}^{', P_names_tex{i}, '} / T_{s,i}^{\textrm{Stupp}}$'], 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    % Panel 3: rho4 vs Ratio
    subplot(1, 3, 3);
    scatter(rho4val, Ratios{i}, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);
    hold on; yline(1.0, '--k', 'LineWidth', 1.0); hold off;
    ylim(y_limits_ratio);
    xlabel('$\rho_4$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(['$T_{s,i}^{', P_names_tex{i}, '} / T_{s,i}^{\textrm{Stupp}}$'], 'Interpreter', 'latex', 'FontSize', 12);
    title(['$', P_names_tex{i}, '$'], 'Interpreter', 'latex', 'FontSize', 13);
    grid on;
    
    safe_name = strrep(strrep(P_names{i}, '(', '_'), ')', '');
    file_name = sprintf('Scatter_Multiplicative_%s.pdf', safe_name);
    exportgraphics(fig_ratio, file_name, 'ContentType', 'vector');
    fprintf('Saved %s\n', file_name);
end
fprintf('\nRatio analysis complete. All figures saved.\n');
% =========================================================================
% 6. QUANTIFICATION OF RELATIVE BENEFIT (LaTeX TABLE DATA)
% =========================================================================
fprintf('\n=================================================================================================\n');
fprintf('   RESPONDER QUANTIFICATION (RELATIVE): R_i^P > 1.05, > 1.10, > 1.20, > 1.30\n');
fprintf('=================================================================================================\n');
fprintf('%-30s | > 1.05 (N, %%) | > 1.10 (N, %%) | > 1.20 (N, %%) | > 1.30 (N, %%)\n', 'Protocol');
fprintf('-------------------------------------------------------------------------------------------------\n');
latex_names = {'3C--Stupp', 'Stupp--3C', '3C--Stupp--3C', ...
    'Stupp(C)--3C--Stupp(A)', 'Stupp(C)--3C--Stupp(A)--3C'};
for i = 1:M
    resp_105 = sum(Ratio_matrix(:, i) > 1.05);
    resp_110 = sum(Ratio_matrix(:, i) > 1.10);
    resp_120 = sum(Ratio_matrix(:, i) > 1.20);
    resp_130 = sum(Ratio_matrix(:, i) > 1.30);
    pct_105 = (resp_105 / N) * 100;
    pct_110 = (resp_110 / N) * 100;
    pct_120 = (resp_120 / N) * 100;
    pct_130 = (resp_130 / N) * 100;
    fprintf('%-30s | %4d (%5.2f%%) | %4d (%5.2f%%) | %4d (%5.2f%%) | %4d (%5.2f%%)\n', ...
        P_names{i}, resp_105, pct_105, resp_110, pct_110, resp_120, pct_120, resp_130, pct_130);
end
fprintf('=================================================================================================\n');
fprintf('\n%% --- Copia y pega estas filas en tu tabla de LaTeX ---\n');
for i = 1:M
    resp_105 = sum(Ratio_matrix(:, i) > 1.05);
    resp_110 = sum(Ratio_matrix(:, i) > 1.10);
    resp_120 = sum(Ratio_matrix(:, i) > 1.20);
    resp_130 = sum(Ratio_matrix(:, i) > 1.30);
    pct_105 = (resp_105 / N) * 100;
    pct_110 = (resp_110 / N) * 100;
    pct_120 = (resp_120 / N) * 100;
    pct_130 = (resp_130 / N) * 100;
    fprintf('        %-30s & %4d & %5.2f & %4d & %5.2f & %4d & %5.2f & %4d & %5.2f \\\\\n', ...
        latex_names{i}, resp_105, pct_105, resp_110, pct_110, resp_120, pct_120, resp_130, pct_130);
end
fprintf('%% -----------------------------------------------------\n\n');