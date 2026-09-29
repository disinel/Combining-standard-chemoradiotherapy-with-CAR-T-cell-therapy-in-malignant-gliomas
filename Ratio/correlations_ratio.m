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

P_names_tex = {'Stupp', '\textrm{3C--Stupp}', '\textrm{Stupp--3C}', '\textrm{3C--Stupp--3C}', ...
               '\textrm{Stupp(C)--3C--Stupp(A)}', '\textrm{Stupp(C)--3C--Stupp(A)--3C}'};
M = length(P_names_tex);

% =========================================================================
% 2. LOAD DATA
% =========================================================================
fprintf('Loading virtual population parameters...\n');
load(file_params); 

fprintf('Loading survival times...\n');
data_stupp = load(file_stupp, 'Tsurv'); Tsurv_Stupp = data_stupp.Tsurv(:);
data_P1    = load(file_P1, 'Tsurv');    Tsurv_P1    = data_P1.Tsurv(:);
data_P2    = load(file_P2, 'Tsurv');    Tsurv_P2    = data_P2.Tsurv(:);
data_P3    = load(file_P3, 'Tsurv');    Tsurv_P3    = data_P3.Tsurv(:);
data_P4    = load(file_P4, 'Tsurv');    Tsurv_P4    = data_P4.Tsurv(:);
data_P5    = load(file_P5, 'Tsurv');    Tsurv_P5    = data_P5.Tsurv(:);

N = length(Tsurv_Stupp);

% Master matrix containing all raw survival times
All_Tsurv = [Tsurv_Stupp, Tsurv_P1, Tsurv_P2, Tsurv_P3, Tsurv_P4, Tsurv_P5];

% =========================================================================
% 3. BUILD PARAMETERS MATRIX (19 parameters)
% =========================================================================
params_matrix = [T0val(:), delta2val(:), delta1val(:), Kival(:), r1val(:), r2val(:), ...
                 alpha1val(:), alpha3val(:), rho1val(:), rho2val(:), rho3val(:), rho4val(:), ...
                 epsilon1val(:), beta1val(:), beta2val(:), tauval(:), gammaval(:), Aval(:), Bval(:)];

param_tex = {'T_0', '\delta_1', '\delta_2', '\text{Ki-67}', 'r_1', 'r_2', ...
             '\alpha_1', '\alpha_3', '\rho_1', '\rho_2', '\rho_3', '\rho_4', ...
             '\epsilon_1', '\beta_1', '\beta_2', '\tau', '\gamma', 'A', 'B'};
num_params = length(param_tex);

% =========================================================================
% 4. CALCULATE CORRELATIONS AND PRINT LATEX TABLE
% =========================================================================
fprintf('\n========================================================================\n');
fprintf('   READY TO PASTE IN OVERLEAF (SPEARMAN: 19 PARAMS vs RAW TSURV):\n');
fprintf('========================================================================\n\n');

fprintf('\\begin{table}[H]\n');
fprintf('    \\centering\n');
fprintf('    \\small\n');
fprintf('    \\begin{tabular}{llrr}\n');
fprintf('        \\toprule\n');
fprintf('        \\textbf{Protocol} & \\textbf{Parameter} & \\textbf{$\\rho$} & \\textbf{$p$-value} \\\\\n');
fprintf('        \\midrule\n');

for k = 1:M
    rhos = zeros(num_params, 1);
    pvals = zeros(num_params, 1);
    
    % Calculate Spearman correlation for each of the 19 parameters against raw Tsurv
    for p = 1:num_params
        [rho_val, p_val] = corr(params_matrix(:, p), All_Tsurv(:, k), 'Type', 'Spearman');
        rhos(p) = rho_val;
        pvals(p) = p_val;
    end
    
    % Sort by absolute correlation value (descending order)
    [~, sort_idx] = sort(abs(rhos), 'descend');
    
    fprintf('        \\multirow{%d}{*}{%s}\n', num_params, P_names_tex{k});
    
    for j = 1:num_params
        idx = sort_idx(j);
        r = rhos(idx);
        p = pvals(idx);
        
        % Elegant formatting for p-values in LaTeX
        if p < 0.001
            p_str = '<0.001';
        elseif p >= 0.01
            p_str = sprintf('%.3f', p);
        else
            str_e = sprintf('%.2e', p);
            p_str = strrep(str_e, 'e', ' \\times 10^{');
            p_str = [p_str, '}'];
            p_str = strrep(p_str, '10^{-0', '10^{-'); % Clean leading zeros inside exponent
        end
        
        fprintf('        & $%s$ & $%.4f$ & $%s$ \\\\\n', param_tex{idx}, r, p_str);
    end
    
    if k < M
        fprintf('        \\midrule\n');
    end
end

fprintf('        \\bottomrule\n');
fprintf('    \\end{tabular}\n');
fprintf('    \\caption{Complete Spearman correlation analysis ($\\rho$) between all 19 model parameters and raw survival time ($T_{surv}$) across different clinical protocols. Parameters are ordered by absolute correlation strength.}\n');
fprintf('    \\label{tab:spearman_correlation_raw_tsurv}\n');
fprintf('\\end{table}\n\n');