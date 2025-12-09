% Gaussian Coefficients Parser from .txt file
clearvars; clc; close all;

txt_file = '../data/igrf14coeffs.txt';
mat_file = '../data/igrf14coeffs.mat';

opts = detectImportOptions(txt_file);
opts.VariableTypes{1} = 'char'; % Force first col (g/h) to be char

T = readtable(txt_file, opts);

col_names = T.Properties.VariableNames;
val_col   = col_names{end-1}; 
sv_col    = col_names{end};  

% Size = 14x14 to support N=13 with (n+1, m+1) indexing
sz = 14; 
g_nm = zeros(sz, sz); h_nm = zeros(sz, sz);
g_sv = zeros(sz, sz); h_sv = zeros(sz, sz);

% Parse Data
for i = 1:height(T)
    
    type = T{i, 1}{1}; % 'g' or 'h'
    n    = T{i, 2};
    m    = T{i, 3};
    val  = T{i, val_col};
    sv   = T{i, sv_col};
    
    if n > 13, continue; end
    
    if type == 'g'
        g_nm(n+1, m+1) = val;
        g_sv(n+1, m+1) = sv;
    elseif type == 'h'
        h_nm(n+1, m+1) = val;
        h_sv(n+1, m+1) = sv;
    end
end

% Saving
% Format: [g, h, g_sv, h_sv, year]
C_nm = zeros(sz, sz*4 + 1);
C_nm(:, 1:sz)         = g_nm;
C_nm(:, sz+1:2*sz)    = h_nm;
C_nm(:, 2*sz+1:3*sz)  = g_sv;
C_nm(:, 3*sz+1:4*sz)  = h_sv;
C_nm(1, end)          = 2025.0; % Base Year

save(mat_file, 'C_nm');