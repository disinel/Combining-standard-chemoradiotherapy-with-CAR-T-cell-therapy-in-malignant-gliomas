% =========================================================================
% Script para calcular Correlación de Spearman con Tsurv (Supervivencia cruda)
% y generar código LaTeX automáticamente.
% =========================================================================
clear variables; close all; clc;

% 1. Archivo de Parámetros
file_params = 'Params_defenitve_with_alphaT_-2026-9-24-17-21_N=10000.mat';

fprintf('Cargando cohorte de pacientes virtuales...\n');
load(file_params);

% 2. Matriz de Parámetros y Nombres LaTeX
% Nota: Aval se mapea a \alpha y Bval a \beta según tu tabla
params_mat = [r1val(:), r2val(:), Kival(:), beta2val(:), Aval(:), Bval(:), ...
              rho1val(:), rho2val(:), rho3val(:), rho4val(:), delta1val(:), delta2val(:), ...
              alpha1val(:), alpha3val(:), epsilon1val(:), gammaval(:), alphaTval(:), beta1val(:), T0val(:)];

param_names = {'$r_1$', '$r_2$', '\text{Ki-67}', '$\beta_2$', '$\alpha$', '$\beta$', ...
               '$\rho_1$', '$\rho_2$', '$\rho_3$', '$\rho_4$', '$\delta_1$', '$\delta_2$', ...
               '$\alpha_1$', '$\alpha_3$', '$\epsilon_1$', '$\gamma$', '$\alpha_T$', '$\beta_1$', '$T_0$'};

% 3. Definición de Protocolos y Archivos
files = {
    'Stupp_Tsurv-2026-9-28-18-1_N=10000.mat', 'Stupp';
    '3C_Stupp_Tsurv-2026-9-28-17-42_N=10000.mat', '3C--Stupp';
    'Stupp_3C_Tsurv-2026-9-28-17-50_N=10000.mat', 'Stupp--3C';
    '3C_Stupp_3C_Tsurv-2026-9-28-17-22_N=10000.mat', '3C--Stupp--3C';
    'Stupp(C)_3C_Stupp(A)_Tsurv-2026-9-28-18-22_N=10000.mat', 'Stupp(C)--3C--Stupp(A)';
    'Stupp(C)_3C_Stupp(A)_3C_Tsurv-2026-9-28-18-10_N=10000.mat', 'Stupp(C)--3C--Stupp(A)--3C'
};

num_protocols = size(files, 1);
Tfinal = 150000; % Para penalizar NaNs si los hubiera

% 4. Generación de Tabla LaTeX
fprintf('\n========================================================================\n');
fprintf('Copia y pega este código en tu editor LaTeX:\n');
fprintf('========================================================================\n\n');

fprintf('\\begin{table}[ht!]\n');
fprintf('    \\caption{Top five parameter groups showing the strongest Spearman rank correlations ($\\rho$) with the survival times $T_{s}$ for each combined protocol. Parameters are ordered by absolute correlation magnitude.}\n');
fprintf('    \\label{tab:spearman_correlation_raw_tsurv}\n');
fprintf('    \\centering\n');
fprintf('    \\small\n');
fprintf('    \\begin{ruledtabular}\n');
fprintf('        \\begin{tabular}{llrr}\n');
fprintf('            \\textbf{Protocol} & \\textbf{Parameter} & \\textbf{$\\rho$} & \\textbf{$p$-value} \\\\\n');
fprintf('            \\hline\n');

for i = 1:num_protocols
    % Cargar Tsurv crudo
    data = load(files{i, 1}, 'Tsurv');
    Tsurv_raw = data.Tsurv(:);
    Tsurv_raw(isnan(Tsurv_raw)) = Tfinal; % Manejo de NaNs
    
    % Calcular correlación de Spearman
    [rho_vals, p_vals] = corr(params_mat, Tsurv_raw, 'Type', 'Spearman', 'Rows', 'complete');
    
    % Ordenar por magnitud absoluta
    [~, sort_idx] = sort(abs(rho_vals), 'descend');
    
    % Agrupar parámetros con correlaciones idénticas (ej. r1 y r2)
    grouped_names = {};
    grouped_rhos  = [];
    grouped_ps    = [];
    skip_idx      = [];
    
    for k = 1:length(sort_idx)
        idx1 = sort_idx(k);
        if ismember(idx1, skip_idx)
            continue;
        end
        
        % Buscar si hay otros parámetros con exactamente la misma correlación
        same_rho_idx = find(abs(rho_vals(sort_idx) - rho_vals(idx1)) < 1e-5);
        actual_indices = sort_idx(same_rho_idx);
        
        % Formatear el nombre agrupado
        if length(actual_indices) > 1
            % Si hay varios, quita los '$' internos para unirlos bien
            clean_names = strrep(param_names(actual_indices), '$', '');
            g_name = ['$', strjoin(clean_names, ','), '$'];
        else
            g_name = param_names{idx1};
        end
        
        grouped_names{end+1} = g_name;
        grouped_rhos(end+1)  = rho_vals(idx1);
        grouped_ps(end+1)    = p_vals(idx1);
        skip_idx = [skip_idx; actual_indices];
    end
    
    % Imprimir fila del protocolo y sus top 5
    prot_name = files{i, 2};
    fprintf('            \\multirow{5}{*}{%s}\n', prot_name);
    
    for top = 1:5
        r_val = grouped_rhos(top);
        p_val = grouped_ps(top);
        
        if p_val < 0.001
            p_str = '$<0.001$';
        else
            p_str = sprintf('$%.3f$', p_val);
        end
        
        % Si es positivo, añadir un espacio invisible o justificar para que alinee con los negativos
        if r_val >= 0
            r_str = sprintf(' $%.4f$', r_val);
        else
            r_str = sprintf('$-%.4f$', abs(r_val));
        end
        
        % La primera línea lleva el '&', las siguientes solo sangría
        if top == 1
            fprintf('            & %-15s & %s & %s \\\\\n', grouped_names{top}, r_str, p_str);
        else
            fprintf('            & %-15s & %s & %s \\\\\n', grouped_names{top}, r_str, p_str);
        end
    end
    
    if i < num_protocols
        fprintf('            \\hline\n');
    end
end

fprintf('        \\end{tabular}\n');
fprintf('    \\end{ruledtabular}\n');
fprintf('\\end{table}\n');