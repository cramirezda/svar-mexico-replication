%% ============================================================================
%  4SVAR_SignRestrictions.m
%  Replicación: Carrillo & Elizondo (2015), "How Robust Are SVARs at
%  Measuring Monetary Policy in Small Open Economies?" — México (Figuras 8-9)
%
%  Estrategia unificada (MATLAB):
%    • 3SVAR_Mexico.do  → Diagnósticos estadísticos (ADF, selección rezagos, estabilidad)
%    • Este script       → Cholesky (Fig. 7) + Sign restrictions + Augmented (Figs. 8-9)
%
%  Identificación:
%    (A) Standard sign restrictions  — Algorithm 1 del paper
%    (B) Augmented sign restrictions — con cota de elasticidad β_lower = -0.7
%
%  Variables endógenas (misma posición que el orden Cholesky):
%    1: igae        — Brecha del producto (output gap)
%    2: pi_core_gap — Brecha inflación subyacente (CPI core)
%    3: pi_inpp_gap — Brecha inflación productor (PPI)
%    4: dep_tcr     — Crecimiento TCR (+ = depreciación)
%    5: i_gap       — Brecha tasa de interés (← choque monetario)
%    6: m2_gap      — Brecha crecimiento del dinero
%
%  Variables exógenas (block-exógenas):
%    7: y_us_gap    — Brecha producto EE.UU.
%    8: pi_us       — Inflación EE.UU. (PCE)
%    9: pi_crb      — Inflación materias primas (CRB)
%
%  Muestra: Enero 2001 – Marzo 2014 | Frecuencia: Mensual | Rezagos: p=3
%% ============================================================================

clear; clc; close all;
rng(42);   % semilla para reproducibilidad

%% ── DIRECTORIOS ─────────────────────────────────────────────────────────────
base_dir = 'C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras';
data_dir  = fullfile(base_dir, 'raw data');
fig_dir   = fullfile(base_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% ============================================================================
%  SECCIÓN 1: CARGA DE DATOS
%% ============================================================================
fprintf('\n===== SECCIÓN 1: CARGA DE DATOS =====\n');

% Leer CSV exportado por STATA (database_for_matlab.csv)
opts = detectImportOptions(fullfile(data_dir, 'database_for_matlab.csv'));
opts = setvartype(opts, opts.VariableNames, 'double');
data_raw = readtable(fullfile(data_dir, 'database_for_matlab.csv'), ...
    'ReadVariableNames', true);

% Lista de variables en orden
endo_names = {'igae','pi_core_gap','pi_inpp_gap','dep_tcr','i_gap','m2_gap'};
exo_names  = {'y_us_gap','pi_us','pi_crb'};

% Extraer matrices
Y = table2array(data_raw(:, endo_names));   % T x 6
Z = table2array(data_raw(:, exo_names));    % T x 3

[T, n] = size(Y);   % T = observaciones, n = 6 variables endógenas
nz     = size(Z, 2);% número de exógenas

fprintf('Muestra: %d observaciones, %d variables endógenas, %d exógenas\n', T, n, nz);

%% ── DIAGNÓSTICO DE ESCALA Y TEST DE NORMALIZACIÓN ───────────────────────────
% dep_tcr usa ×1200 (% anual); el paper C&E usa ×100 (% mensual).
% Ratio P[dep_tcr]/P[i_gap] = ~40 en nuestra versión; debería ser ~1-3.
% Si test_tcr_monthly=true: dep_tcr ÷ 12 → % mensual → P[dep_tcr] ≈ 1.7
% Evaluar si esto corrige la forma del IRF de output gap (Figs. 8-9).
%
% NOTA: beta_lower NO se ve afectado (usa solo Y(:,2) y Y(:,5)).
% NOTA: Cambio permanente requiere modificar 0Proc_data.R (dep_tcr ×100).
test_tcr_monthly = true;   % <══ CAMBIAR A true PARA EJECUTAR EL TEST
fig_suffix = '';
if test_tcr_monthly
    Y(:, 4) = Y(:, 4) / 12;
    fig_suffix = '_test_monthly';
    fprintf('  [TEST ESCALA] dep_tcr dividido por 12 (tasa mensual %%/mes)\n');
    fprintf('  dep_tcr nuevo: sd=%.4f  min=%.4f  max=%.4f\n', ...
            std(Y(:,4)), min(Y(:,4)), max(Y(:,4)));
end

fprintf('\n--- Diagnóstico de escala por variable ---\n');
for vi = 1:n
    fprintf('  %-15s  sd=%8.4f\n', endo_names{vi}, std(Y(:,vi)));
end
fprintf('\n');

%% ============================================================================
%  SECCIÓN 2: ESTIMACIÓN DEL VAR(3) — MÁX. VEROSIMILITUD (OLS por ecuación)
%  Modelos: (a) Simple VAR, (b) VAR con exógenas
%% ============================================================================
p = 3;    % rezagos del VAR (paper: 3)

function [B, Sigma, resid, X] = estimate_var(Y, Z, p, use_exo)
    % Estima VAR(p) por OLS, con o sin variables exógenas
    % Y: T x n, Z: T x nz, p: rezagos, use_exo: 0/1
    % Retorna: B (coefs), Sigma (cov residuos), resid (T-p x n), X (regresores)
    [T, n] = size(Y);
    nz = size(Z, 2);

    % Construir matriz de regresores X = [const | Y_(t-1),...,Y_(t-p) | Z_(t)]
    Y_lag = [];
    for lag = 1:p
        Y_lag = [Y_lag, [zeros(lag, n); Y(1:end-lag, :)]];
    end
    if use_exo
        X = [ones(T, 1), Y_lag, Z];
    else
        X = [ones(T, 1), Y_lag];
    end

    % Quitar primeras p filas (sin suficientes rezagos)
    Yp = Y(p+1:end, :);
    Xp = X(p+1:end, :);
    Zp = Z(p+1:end, :);

    if use_exo
        Xp = [ones(size(Yp,1),1), Y_lag(p+1:end,:), Zp];
    else
        Xp = [ones(size(Yp,1),1), Y_lag(p+1:end,:)];
    end

    % OLS: B = (X'X)^{-1} X'Y
    B = (Xp' * Xp) \ (Xp' * Yp);   % k x n
    resid = Yp - Xp * B;             % (T-p) x n
    Teff = size(resid, 1);
    k = size(Xp, 2);
    Sigma = (resid' * resid) / (Teff - k); % n x n, sin ajuste df
    X = Xp;
end

fprintf('\n--- Estimando VAR Simple ---\n');
[B_s, Sigma_s, resid_s, X_s] = estimate_var(Y, Z, p, 0);

fprintf('--- Estimando VAR con Exógenas ---\n');
[B_e, Sigma_e, resid_e, X_e] = estimate_var(Y, Z, p, 1);

%% ============================================================================
%  PARÁMETROS GLOBALES
%  Compartidos por Sección 4B (Cholesky) y Sección 5 (restricciones de signo)
%% ============================================================================
h        = 16;    % horizonte de IRF en meses
n_boot   = 2000;  % réplicas bootstrap (paper: 2000)
K        = 100;   % modelos aceptados por draw (paper: 100)
K_max    = 5000;  % candidatos máximos por draw
ci_lvl   = 68;    % nivel de IC en % (paper: 68%)
j_mp     = 5;     % índice del choque monetario (i_gap, posición 5)

% Restricciones de signo — Tabla 2 del paper
% y(-), pi_core(-), pi_inpp(libre), dep_tcr(-), i(+), m2(-)
% 0 = sin restricción; -1 = negativo; +1 = positivo
signs_mp = [-1; -1; 0; -1; 1; -1];

%% ============================================================================
%  SECCIÓN 3: ESTIMACIÓN DE β_lower (Ecuación 11 del paper)
%  Modelo: Δπ_t = β0 + βi·Δi_t + βii·D_t·Δi_t + ω_t, ω_t ~ ARMA(1,2)
%
%  Método principal: OLS + Newey-West HAC (lag=3) — replicando C&E (2015)
%  Método comparación: MLE ARMA(1,2) conjunta (más eficiente, SE más chico)
%
%  El paper reporta β̂⁻ = −0.035 y β_lower = −0.7, consistente con OLS+NW.
%  MLE da SE ~4× más chico → β_lower ~−0.18 (no replica el paper).
%% ============================================================================
fprintf('\n===== SECCIÓN 3: ESTIMACIÓN β_lower (OLS + Newey-West, C&E 2015) =====\n');

dpi    = diff(Y(:, 2));   % Δπ_core  (T-1 x 1)
di_vec = diff(Y(:, 5));   % Δi_gap   (T-1 x 1)
D_t    = double(dpi .* di_vec < 0);
di_dt  = D_t .* di_vec;

T_reg = length(dpi);
X_reg = [ones(T_reg, 1), di_vec, di_dt];

% ── Paso 1: OLS ──────────────────────────────────────────────────────────────
b_ols  = (X_reg' * X_reg) \ (X_reg' * dpi);
omega0 = dpi - X_reg * b_ols;
beta_minus_ols = b_ols(2) + b_ols(3);

% ── Paso 2 (PRIMARIO): Newey-West HAC con lag=3 ──────────────────────────────
% lag=3 determinado por estructura ARMA(1,2) del residuo (C&E 2015)
lags_nw = 3;
XX = X_reg' * X_reg;
S  = XX;
for lag = 1:lags_nw
    wt     = 1 - lag/(lags_nw + 1);
    Xp_lag = X_reg(1+lag:end, :);
    Xp_cur = X_reg(1:end-lag, :);
    ep_lag = omega0(1+lag:end);
    ep_cur = omega0(1:end-lag);
    Gamma  = (Xp_cur .* ep_cur)' * (Xp_lag .* ep_lag);
    S      = S + wt * (Gamma + Gamma');
end
Var_nw        = (XX \ S) / XX;
se_nw         = sqrt(max(0, Var_nw(2,2) + Var_nw(3,3) + 2*Var_nw(2,3)));
beta_lower_nw = beta_minus_ols - 1.645 * se_nw;

fprintf('  [PRIMARIO] OLS + Newey-West (lag=3):\n');
fprintf('    β̂⁻ puntual (βi+βii) = %8.4f  (paper: -0.035)\n', beta_minus_ols);
fprintf('    SE β̂⁻ (NW)          = %8.4f  (paper implícito: ~0.40)\n', se_nw);
fprintf('    β_lower (IC 95%%)   = %8.4f  (paper: -0.700)\n', beta_lower_nw);

% ── Paso 3 (COMPARACIÓN): MLE conjunta ARMA(1,2) ────────────────────────────
fprintf('  [COMPARACIÓN] MLE ARMA(1,2):\n');

phi1_init = sum(omega0(2:end) .* omega0(1:end-1)) / sum(omega0(1:end-1).^2);
phi1_init = max(-0.9, min(0.9, phi1_init));
init_p = [b_ols; phi1_init; 0.10; 0.05; log(var(omega0))];

nll_fn = @(p) beta_lower_nll(p, dpi, X_reg);
opts_mle = optimoptions('fminunc', ...
    'Display',              'off', ...
    'Algorithm',            'quasi-newton', ...
    'MaxIterations',        5000,  ...
    'OptimalityTolerance',  1e-9,  ...
    'StepTolerance',        1e-10);

[p_hat, ~, exitflag_mle] = fminunc(nll_fn, init_p, opts_mle);
if exitflag_mle <= 0
    [p_hat, ~, exitflag_mle] = fminunc(nll_fn, p_hat, opts_mle);
end

phi1_hat   = p_hat(4);
theta1_hat = p_hat(5);
theta2_hat = p_hat(6);
sigma2_hat = exp(p_hat(7));
beta_minus_mle = p_hat(2) + p_hat(3);

h_fd = 1e-5; n_p = length(p_hat);
H = zeros(n_p, n_p); f0 = nll_fn(p_hat);
for ii = 1:n_p
    for jj = ii:n_p
        ei = zeros(n_p,1); ei(ii) = h_fd;
        ej = zeros(n_p,1); ej(jj) = h_fd;
        if ii == jj
            H(ii,ii) = (nll_fn(p_hat+ei) - 2*f0 + nll_fn(p_hat-ei)) / h_fd^2;
        else
            H(ii,jj) = (nll_fn(p_hat+ei+ej) - nll_fn(p_hat+ei-ej) - ...
                        nll_fn(p_hat-ei+ej) + nll_fn(p_hat-ei-ej)) / (4*h_fd^2);
            H(jj,ii) = H(ii,jj);
        end
    end
end
Var_mle = inv(H);
se_mle  = sqrt(max(0, Var_mle(2,2) + Var_mle(3,3) + 2*Var_mle(2,3)));
beta_lower_mle = beta_minus_mle - 1.645 * se_mle;

fprintf('    β̂⁻ puntual (MLE)    = %8.4f\n', beta_minus_mle);
fprintf('    SE β̂⁻ (Hessiano)    = %8.4f\n', se_mle);
fprintf('    β_lower (IC 95%%)   = %8.4f\n', beta_lower_mle);
fprintf('    φ1=%.4f  θ1=%.4f  θ2=%.4f  σ²=%.6f  exitflag=%d\n', ...
    phi1_hat, theta1_hat, theta2_hat, sigma2_hat, exitflag_mle);

% ── Paso 4: Selección del β_lower final (OLS+NW como primario) ───────────────
beta_minus    = beta_minus_ols;
se_beta_minus = se_nw;
beta_lower    = beta_lower_nw;
metodo_usado  = 'OLS+NW';

if beta_lower > 0 || beta_lower < -5
    fprintf('  ⚠ OLS+NW fuera de rango — usando MLE\n');
    beta_minus    = beta_minus_mle;
    se_beta_minus = se_mle;
    beta_lower    = beta_lower_mle;
    metodo_usado  = 'MLE';
    if beta_lower > 0 || beta_lower < -5
        beta_lower = -0.7;
        metodo_usado = 'paper_fallback';
        fprintf('  ⚠ MLE también fuera de rango → β_lower=-0.7 (paper)\n');
    end
end
fprintf('  → β_lower usado: %.4f  [método: %s]\n\n', beta_lower, metodo_usado);

% ── Exportar para 4Informe_replicacion.R ─────────────────────────────────────
exitflag = exitflag_mle;
bl_tbl = table(beta_minus, se_beta_minus, beta_lower, ...
               beta_minus_mle, se_mle, beta_lower_mle, ...
               phi1_hat, theta1_hat, theta2_hat, sigma2_hat, exitflag, ...
    'VariableNames', {'beta_minus','se_beta_minus','beta_lower', ...
                      'beta_minus_mle','se_mle','beta_lower_mle', ...
                      'phi1','theta1','theta2','sigma2','exitflag'});
writetable(bl_tbl, fullfile(data_dir, 'beta_lower_arma12.csv'));
fprintf('  Exportado: %s\n\n', fullfile(data_dir, 'beta_lower_arma12.csv'));

%% ============================================================================
%  SECCIÓN 4: FUNCIONES AUXILIARES
%% ============================================================================

function Q = haar_random_orthogonal(n)
    % Genera una matriz ortogonal aleatoria n×n de la distribución de Haar
    Z = randn(n, n);
    [Q, R] = qr(Z);
    % Ajustar signos para distribución uniforme sobre O(n)
    d = sign(diag(R));
    d(d == 0) = 1;
    Q = Q * diag(d);
end

function irf_mat = compute_irf(B, n, p, shock_vec, h)
    % Propaga la IRF desde el impacto (h=0) hasta h períodos
    % B: (1 + n*p) x n (coefs del VAR, primera fila = constante)
    % shock_vec: n×1 (respuesta de impacto = columna j de A0)
    % h: horizonte de IRF
    
    % Extraer los bloques de coeficientes [A1, A2, ..., Ap]
    % B tiene la forma: [const | A1 | A2 | ... | Ap]
    % donde cada Aj es n×n
    Acoefs = cell(1, p);
    for j = 1:p
        row_start = 1 + (j-1)*n + 1;   % +1 para saltar la constante
        row_end   = 1 + j*n;
        Acoefs{j} = B(row_start:row_end, :)';   % n×n (transponer: B es k×n)
    end

    % Matriz de IRF: (h+1) × n
    irf_mat = zeros(h+1, n);
    irf_mat(1, :) = shock_vec';   % período de impacto h=0

    % Propagar: y_{t+h} = A1*y_{t+h-1} + A2*y_{t+h-2} + ... + Ap*y_{t+h-p}
    for hh = 2:(h+1)
        irf_h = zeros(n, 1);
        for j = 1:p
            if hh - j >= 1
                irf_h = irf_h + Acoefs{j} * irf_mat(hh-j, :)';
            end
        end
        irf_mat(hh, :) = irf_h';
    end
end

function [irf_med, irf_lo, irf_hi] = sign_restricted_svar(...
        Y, Z, p, use_exo, n_boot, K, K_max, h, ci_level, ...
        type, j_mp, signs_mp, beta_lower_val)
    % Algoritmo principal del paper (Algorithm 1):
    %   - Bootstrap de residuos (Runkle 1987)
    %   - Rotaciones ortogonales aleatorias (Haar measure)
    %   - Verificación de signos de impacto (Tabla 2)
    %   - (Aumentado): verificación de cota de elasticidad
    %   - Propagación completa de IRF (h períodos)
    %   - Bandas de incertidumbre al ci_level%
    %
    % Retorna: irf_med, irf_lo, irf_hi — cada uno (h+1) × n
    
    [T, n] = size(Y);
    nfore  = h + 1;
    
    % Estimar VAR base
    [B0, Sigma0, resid0, X0] = estimate_var(Y, Z, p, use_exo);
    P0 = chol(Sigma0, 'lower');   % descomposición de Cholesky de Sigma
    
    % Almacén de IRFs aceptadas: dim = (n_boot * K) × (n * nfore)
    % Para memoria más eficiente, usamos un cell array por variable
    all_irfs = zeros(n_boot * K, n, nfore);
    total_accepted = 0;
    
    fprintf('    Bootstrap [%s, use_exo=%d]: ', type, use_exo);
    
    for b = 1:n_boot
        % 1. Bootstrap de residuos con remplazo (Runkle 1987)
        T_eff   = size(resid0, 1);
        idx     = randi(T_eff, T_eff, 1);
        res_b   = resid0(idx, :);
        
        % 2. Re-estimar Sigma (covarianza muestral de los residuos bootstrap)
        Sigma_b = (res_b' * res_b) / T_eff;
        Sigma_b = (Sigma_b + Sigma_b') / 2;   % simetrizar
        
        % Cholesky del bootstrap (con jitter para robustez numérica)
        jitter = 1e-10 * trace(Sigma_b) / n;
        [P_b, flag] = chol(Sigma_b + jitter*eye(n), 'lower');
        if flag ~= 0
            P_b = P0;   % fallback a la Cholesky del VAR base
        end
        
        % 3. Buscar matrices ortogonales que cumplan los signos
        k_acc = 0;
        irfs_draw = zeros(K, n, nfore);
        
        for k = 1:K_max
            % Generar Q ~ Haar (uniforme sobre O(n))
            Q = haar_random_orthogonal(n);
            
            % Candidato: A0 = P_b * Q
            A0 = P_b * Q;
            
            % Respuesta de impacto al choque monetario (columna j_mp)
            resp = A0(:, j_mp);   % n × 1
            
            % 4a. Verificar restricciones de signos
            valid = true;
            for vv = 1:n
                if signs_mp(vv) ~= 0
                    if signs_mp(vv) * resp(vv) <= 0
                        valid = false;
                        break;
                    end
                end
            end
            if ~valid, continue; end
            
            % 4b. (Solo aumentadas) Verificar cota de elasticidad
            %     β_lower * Δi ≤ Δπ ≤ 0
            %     resp(2) = Δπ (inflación), resp(5) = Δi (tasa)
            if strcmp(type, 'augmented')
                dpi_imp = resp(2);   % impacto en inflación (debe ser < 0)
                di_imp  = resp(j_mp);  % impacto en tasa de interés (debe ser > 0)
                if beta_lower_val * di_imp > dpi_imp
                    continue;   % fuera de la cota
                end
            end
            
            % 5. Calcular IRF completa (h períodos) para este A0
            irf_full = compute_irf(B0, n, p, resp, h);   % (h+1) × n
            
            k_acc = k_acc + 1;
            irfs_draw(k_acc, :, :) = irf_full';   % n × (h+1)
            
            if k_acc >= K, break; end
        end   % fin loop candidatos
        
        % Acumular modelos aceptados en este draw
        if k_acc > 0
            all_irfs(total_accepted+1:total_accepted+k_acc, :, :) = ...
                irfs_draw(1:k_acc, :, :);
            total_accepted = total_accepted + k_acc;
        end
        
        % Progreso cada 200 bootstraps
        if mod(b, 200) == 0
            fprintf('%d..', b);
        end
    end   % fin loop bootstrap
    
    fprintf(' (%d modelos aceptados de %d draws)\n', total_accepted, n_boot);
    
    % 6. Calcular bandas de incertidumbre (percentiles)
    if total_accepted == 0
        warning('Ningún modelo aceptado. Revisa las restricciones de signos.');
        irf_med = zeros(nfore, n);
        irf_lo  = zeros(nfore, n);
        irf_hi  = zeros(nfore, n);
        return;
    end
    
    irf_data = all_irfs(1:total_accepted, :, :);   % accepted × n × nfore
    pct_lo = (100 - ci_level) / 2;
    pct_hi = 100 - pct_lo;
    
    irf_med = zeros(nfore, n);
    irf_lo  = zeros(nfore, n);
    irf_hi  = zeros(nfore, n);

    for vv = 1:n
        for hh = 1:nfore
            col = squeeze(irf_data(:, vv, hh));
            irf_med(hh, vv) = median(col);
            irf_lo(hh, vv)  = prctile(col, pct_lo);
            irf_hi(hh, vv)  = prctile(col, pct_hi);
        end
    end
end

%% ============================================================================
%  SECCIÓN 4B: RESTRICCIONES DE EXCLUSIÓN (CHOLESKY) — Figura 7
%  Identificación recursiva (Cushman & Zha 1997), bootstrap Runkle (1987)
%  Modelos: (a) VAR simple, (b) VAR con exógenas block-exógenas
%% ============================================================================
fprintf('\n===== SECCIÓN 4B: CHOLESKY (FIGURA 7) =====\n');

pct_lo_c = (100 - ci_lvl) / 2;
pct_hi_c = 100 - pct_lo_c;

% ── Runkle (1987) recursive bootstrap: re-estimates VAR in each draw ─────────
% Step 1: generate pseudo-data Y_boot from VAR recursion with resampled residuals
% Step 2: re-estimate VAR on Y_boot to get B_boot, Sigma_boot
% Step 3: Cholesky of Sigma_boot → shock_boot
% Step 4: IRF from (B_boot, shock_boot)
% This captures both coefficient and shock-identification uncertainty.

function Y_boot = runkle_bootstrap_data(B, Y, Z, p, res_boot, use_exo)
    % Generates pseudo-data via VAR recursion (Runkle 1987)
    % B: regressor coefficient matrix (k x n), from estimate_var
    % Y: T x n original data (initial conditions kept for first p rows)
    % Z: T x nz exogenous
    % res_boot: T_eff x n resampled residuals
    [T, n]  = size(Y);
    nz      = size(Z, 2);
    Y_boot  = Y;   % copy — first p rows are initial conditions
    for t = p+1:T
        y_lags = [];
        for lag = 1:p
            y_lags = [y_lags, Y_boot(t-lag, :)];
        end
        if use_exo
            x_t = [1, y_lags, Z(t,:)];
        else
            x_t = [1, y_lags];
        end
        Y_boot(t, :) = x_t * B + res_boot(t-p, :);
    end
end

% ── (a) VAR Simple (sin exógenas) ────────────────────────────────────────────
fprintf('   Bootstrap Cholesky — VAR Simple (%d réplicas, Runkle 1987)...\n', n_boot);
[B_chol_s, Sigma_chol_s, resid_chol_s, ~] = estimate_var(Y, Z, p, 0);
P_chol_s = chol(Sigma_chol_s, 'lower');
T_eff_s  = size(resid_chol_s, 1);

irfs_chol_smp_boot = zeros(n_boot, n, h+1);
for b = 1:n_boot
    idx      = randi(T_eff_s, T_eff_s, 1);
    res_b    = resid_chol_s(idx, :);
    Y_b      = runkle_bootstrap_data(B_chol_s, Y, Z, p, res_b, 0);
    [B_b, Sig_b, ~, ~] = estimate_var(Y_b, Z, p, 0);
    Sig_b    = (Sig_b + Sig_b') / 2;
    [P_b, flag] = chol(Sig_b + 1e-10*eye(n), 'lower');
    if flag ~= 0, P_b = P_chol_s; end
    shock_b  = P_b(:, j_mp);
    irf_b    = compute_irf(B_b, n, p, shock_b, h);
    irfs_chol_smp_boot(b, :, :) = irf_b';
end

irf_chol_smp_med = zeros(h+1, n);
irf_chol_smp_lo  = zeros(h+1, n);
irf_chol_smp_hi  = zeros(h+1, n);
for vv = 1:n
    for hh = 1:h+1
        col = squeeze(irfs_chol_smp_boot(:, vv, hh));
        irf_chol_smp_med(hh, vv) = median(col);
        irf_chol_smp_lo(hh, vv)  = prctile(col, pct_lo_c);
        irf_chol_smp_hi(hh, vv)  = prctile(col, pct_hi_c);
    end
end
fprintf('   → VAR Simple completado.\n');

% ── (b) VAR con exógenas ──────────────────────────────────────────────────────
fprintf('   Bootstrap Cholesky — VAR con exógenas (%d réplicas, Runkle 1987)...\n', n_boot);
[B_chol_e, Sigma_chol_e, resid_chol_e, ~] = estimate_var(Y, Z, p, 1);
P_chol_e = chol(Sigma_chol_e, 'lower');
T_eff_e  = size(resid_chol_e, 1);

irfs_chol_exo_boot = zeros(n_boot, n, h+1);
for b = 1:n_boot
    idx      = randi(T_eff_e, T_eff_e, 1);
    res_b    = resid_chol_e(idx, :);
    Y_b      = runkle_bootstrap_data(B_chol_e, Y, Z, p, res_b, 1);
    [B_b, Sig_b, ~, ~] = estimate_var(Y_b, Z, p, 1);
    Sig_b    = (Sig_b + Sig_b') / 2;
    [P_b, flag] = chol(Sig_b + 1e-10*eye(n), 'lower');
    if flag ~= 0, P_b = P_chol_e; end
    shock_b  = P_b(:, j_mp);
    irf_b    = compute_irf(B_b, n, p, shock_b, h);
    irfs_chol_exo_boot(b, :, :) = irf_b';
end

irf_chol_exo_med = zeros(h+1, n);
irf_chol_exo_lo  = zeros(h+1, n);
irf_chol_exo_hi  = zeros(h+1, n);
for vv = 1:n
    for hh = 1:h+1
        col = squeeze(irfs_chol_exo_boot(:, vv, hh));
        irf_chol_exo_med(hh, vv) = median(col);
        irf_chol_exo_lo(hh, vv)  = prctile(col, pct_lo_c);
        irf_chol_exo_hi(hh, vv)  = prctile(col, pct_hi_c);
    end
end
fprintf('   → VAR con exógenas completado.\n');

%% ============================================================================
%  SECCIÓN 5: EJECUCIÓN DE LOS 4 MODELOS
%  (standard/augmented) × (simple/exo)
%% ============================================================================
fprintf('\n===== SECCIÓN 5: CORRIENDO MODELOS =====\n');

% ── Diagnóstico P diagonal (choque de impacto por variable) ──────────────────
[~, Sigma_diag, ~, ~] = estimate_var(Y, Z, p, 1);
P_diag = chol(Sigma_diag, 'lower');
fprintf('--- P diagonal (SD innovaciones con exógenas) ---\n');
for vi = 1:n
    fprintf('  %-15s  P[%d,%d] = %8.4f\n', endo_names{vi}, vi, vi, P_diag(vi,vi));
end
fprintf('  Ratio P[dep_tcr]/P[i_gap] = %.2f  (paper: ≈ 3–5)\n\n', ...
        P_diag(4,4) / P_diag(5,5));

fprintf('\n[A] Standard Sign Restrictions — VAR con exógenas\n');
tic;
[irf_std_exo_med, irf_std_exo_lo, irf_std_exo_hi] = ...
    sign_restricted_svar(Y, Z, p, 1, n_boot, K, K_max, h, ci_lvl, ...
                         'standard', j_mp, signs_mp, beta_lower);
toc;

fprintf('\n[B] Standard Sign Restrictions — VAR simple\n');
tic;
[irf_std_smp_med, irf_std_smp_lo, irf_std_smp_hi] = ...
    sign_restricted_svar(Y, Z, p, 0, n_boot, K, K_max, h, ci_lvl, ...
                         'standard', j_mp, signs_mp, beta_lower);
toc;

fprintf('\n[C] Augmented Sign Restrictions — VAR con exógenas\n');
tic;
[irf_aug_exo_med, irf_aug_exo_lo, irf_aug_exo_hi] = ...
    sign_restricted_svar(Y, Z, p, 1, n_boot, K, K_max, h, ci_lvl, ...
                         'augmented', j_mp, signs_mp, beta_lower);
toc;

fprintf('\n[D] Augmented Sign Restrictions — VAR simple\n');
tic;
[irf_aug_smp_med, irf_aug_smp_lo, irf_aug_smp_hi] = ...
    sign_restricted_svar(Y, Z, p, 0, n_boot, K, K_max, h, ci_lvl, ...
                         'augmented', j_mp, signs_mp, beta_lower);
toc;

%% ============================================================================
%  SECCIÓN 6: GRÁFICAS — FIGURAS 8 y 9
%% ============================================================================
fprintf('\n===== SECCIÓN 6: GRÁFICAS =====\n');

% Eje x del horizonte
x_axis = 0:h;

% Etiquetas de variables para los paneles
var_labels = {'Output', 'Inflation (core)', 'PPI', 'Real ER dep', ...
              'Nom. interest rate', 'Money growth'};

% Colores
col_exo    = [0.12, 0.47, 0.71];  % azul oscuro (VAR w/ exo)
col_simple = [0.50, 0.50, 0.50];  % gris (Simple VAR)
col_aug    = [0.20, 0.63, 0.17];  % verde (augmented)
col_std    = [0.69, 0.19, 0.38];  % rojo/magenta (standard)

alpha_fill = 0.25;   % transparencia de las bandas

function plot_irf_panel(ax, x, med, lo, hi, col, alpha, line_style, label)
    % Dibuja la IRF con banda sombreada en los ejes ax
    fill(ax, [x, fliplr(x)], [lo, fliplr(hi)], col, ...
         'FaceAlpha', alpha, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x, med, 'Color', col, 'LineWidth', 1.5, ...
         'LineStyle', line_style, 'DisplayName', label);
    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, max(x)]);
    ax.XTick = 0:4:max(x);
    ax.FontSize = 8;
    box(ax, 'on');
    grid(ax, 'off');
end

%% ── FIGURA 7: Exclusion Restrictions (Cholesky) ─────────────────────────────
% Layout: 2 filas × 3 columnas (6 variables)
% Azul navy (sólido) = VAR con exógenas | Gris (punteado) = VAR simple

fig7 = figure('Name', 'Figure 7: Exclusion Restrictions', ...
    'Position', [50, 50, 1400, 600], 'Color', 'white');

for vv = 1:6
    ax = subplot(2, 3, vv);

    % VAR con exógenas (azul navy)
    area_x = [x_axis, fliplr(x_axis)];
    area_y = [irf_chol_exo_lo(:, vv)', fliplr(irf_chol_exo_hi(:, vv)')];
    fill(ax, area_x, area_y, col_exo, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x_axis, irf_chol_exo_med(:, vv), 'Color', col_exo, 'LineWidth', 1.5);

    % VAR simple (gris punteado)
    area_y_s = [irf_chol_smp_lo(:, vv)', fliplr(irf_chol_smp_hi(:, vv)')];
    fill(ax, area_x, area_y_s, col_simple, 'FaceAlpha', 0.20, ...
         'EdgeColor', col_simple, 'LineStyle', '--', 'LineWidth', 0.5);
    plot(ax, x_axis, irf_chol_smp_med(:, vv), '--', 'Color', col_simple, 'LineWidth', 1.2);

    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, h]);
    ax.XTick = 0:4:h;
    ax.FontSize = 7.5;
    box(ax, 'on');

    title(ax, var_labels{vv}, 'FontSize', 8.5, 'FontWeight', 'bold');
    if vv == 1
        ylabel(ax, 'Pct pts from trend', 'FontSize', 7.5);
        xlabel(ax, 'Months after shock', 'FontSize', 7.5);
    end
end

annotation(fig7, 'textbox', [0.68, 0.96, 0.32, 0.04], ...
    'String', '■ Blue = VAR w/ exo   ■ Gray dashed = Simple VAR', ...
    'FontSize', 7.5, 'EdgeColor', 'none');

% sgtitle eliminado — título se pone desde el .tex via \caption{}

exportgraphics(fig7, fullfile(fig_dir, 'fig7_exclusion.pdf'), 'Resolution', 300);
exportgraphics(fig7, fullfile(fig_dir, 'fig7_exclusion.png'), 'Resolution', 200);
fprintf('→ Figura 7 guardada: %s\n', fullfile(fig_dir, 'fig7_exclusion.pdf'));

%% ── FIGURA 8: Standard vs Augmented (Simple VAR arriba, w/Exo abajo) ──────
% Layout: 2 paneles verticales × 6 columnas
% Top row:    Standard sign restrictions — Simple VAR vs. VAR w/ exo
% Bottom row: Augmented sign restrictions — Simple VAR vs. VAR w/ exo

fig8 = figure('Name', 'Figure 8: Sign Restrictions', ...
    'Position', [50, 50, 1400, 700], 'Color', 'white');

% ── Panel Superior: STANDARD SIGN RESTRICTIONS ────────────────────────────
for vv = 1:6
    ax = subplot(2, 6, vv);

    % VAR con exógenas (azul)
    area_x = [x_axis, fliplr(x_axis)];
    area_y = [irf_std_exo_lo(:, vv)', fliplr(irf_std_exo_hi(:, vv)')];
    fill(ax, area_x, area_y, col_exo, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x_axis, irf_std_exo_med(:, vv), 'Color', col_exo, 'LineWidth', 1.5);

    % VAR simple (gris punteado)
    area_y_s = [irf_std_smp_lo(:, vv)', fliplr(irf_std_smp_hi(:, vv)')];
    fill(ax, area_x, area_y_s, col_simple, 'FaceAlpha', 0.2, 'EdgeColor', col_simple, 'LineStyle', '--', 'LineWidth', 0.5);
    plot(ax, x_axis, irf_std_smp_med(:, vv), '--', 'Color', col_simple, 'LineWidth', 1.2);

    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, h]);
    ax.XTick = 0:4:h;
    ax.FontSize = 7.5;
    box(ax, 'on');

    title(ax, var_labels{vv}, 'FontSize', 8.5, 'FontWeight', 'bold');
    if vv == 1
        ylabel(ax, 'Pct pts from trend', 'FontSize', 7.5);
        xlabel(ax, 'Months after shock', 'FontSize', 7.5);
    end
end

% ── Panel Inferior: AUGMENTED SIGN RESTRICTIONS ───────────────────────────
for vv = 1:6
    ax = subplot(2, 6, vv + 6);

    % VAR con exógenas (azul)
    area_x = [x_axis, fliplr(x_axis)];
    area_y = [irf_aug_exo_lo(:, vv)', fliplr(irf_aug_exo_hi(:, vv)')];
    fill(ax, area_x, area_y, col_exo, 'FaceAlpha', alpha_fill, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x_axis, irf_aug_exo_med(:, vv), 'Color', col_exo, 'LineWidth', 1.5);

    % VAR simple (gris)
    area_y_s = [irf_aug_smp_lo(:, vv)', fliplr(irf_aug_smp_hi(:, vv)')];
    fill(ax, area_x, area_y_s, col_simple, 'FaceAlpha', 0.2, 'EdgeColor', col_simple, 'LineStyle', '--', 'LineWidth', 0.5);
    plot(ax, x_axis, irf_aug_smp_med(:, vv), '--', 'Color', col_simple, 'LineWidth', 1.2);

    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, h]);
    ax.XTick = 0:4:h;
    ax.FontSize = 7.5;
    box(ax, 'on');

    title(ax, var_labels{vv}, 'FontSize', 8.5, 'FontWeight', 'bold');
    if vv == 1
        ylabel(ax, 'Pct pts from trend', 'FontSize', 7.5);
        xlabel(ax, 'Months after shock', 'FontSize', 7.5);
    end
end

% Anotaciones de fila
annotation(fig8, 'textbox', [0.01, 0.52, 0.15, 0.04], ...
    'String', '(Standard Sign Restrictions)', 'FontSize', 8, ...
    'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', col_exo);
annotation(fig8, 'textbox', [0.01, 0.03, 0.15, 0.04], ...
    'String', '(Augmented Sign Restrictions)', 'FontSize', 8, ...
    'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', col_exo);

% Leyenda manual
annotation(fig8, 'textbox', [0.78, 0.96, 0.22, 0.04], ...
    'String', '■ Blue = VAR w/ exo   ■ Gray = Simple VAR', ...
    'FontSize', 7.5, 'EdgeColor', 'none');

% sgtitle eliminado — título se pone desde el .tex via \caption{}

% Exportar
exportgraphics(fig8, fullfile(fig_dir, ['fig8_sign_restrictions' fig_suffix '.pdf']), 'Resolution', 300);
exportgraphics(fig8, fullfile(fig_dir, ['fig8_sign_restrictions' fig_suffix '.png']), 'Resolution', 200);
fprintf('→ Figura 8 guardada: %s\n', fullfile(fig_dir, ['fig8_sign_restrictions' fig_suffix '.pdf']));

%% ── FIGURA 9: Comparación entre estrategias (VAR w/ exo solamente) ─────────
% Top panel:    Exclusion (Cholesky) vs Augmented sign restrictions
% Bottom panel: Standard signs vs Augmented sign restrictions
% Nota: irf_chol_exo_* proviene de Sección 4B (ya calculado)

% ── Figura 9 ─────────────────────────────────────────────────────────────────
fig9 = figure('Name', 'Figure 9: Comparison of Identifying Assumptions', ...
    'Position', [50, 50, 1400, 700], 'Color', 'white');

% ── Panel Superior: Exclusion vs. Augmented Signs ────────────────────────
for vv = 1:6
    ax = subplot(2, 6, vv);

    % Exclusion / Cholesky (azul)
    area_x = [x_axis, fliplr(x_axis)];
    area_y = [irf_chol_exo_lo(:, vv)', fliplr(irf_chol_exo_hi(:, vv)')];
    fill(ax, area_x, area_y, col_exo, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x_axis, irf_chol_exo_med(:, vv), 'Color', col_exo, 'LineWidth', 1.5);

    % Augmented signs (verde, línea punteada con banda)
    area_y_a = [irf_aug_exo_lo(:, vv)', fliplr(irf_aug_exo_hi(:, vv)')];
    fill(ax, area_x, area_y_a, col_aug, 'FaceAlpha', 0.15, 'EdgeColor', col_aug, 'LineStyle', ':', 'LineWidth', 0.5);
    plot(ax, x_axis, irf_aug_exo_med(:, vv), ':', 'Color', col_aug, 'LineWidth', 1.8);

    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, h]);
    ax.XTick = 0:4:h;
    ax.FontSize = 7.5;
    box(ax, 'on');

    title(ax, var_labels{vv}, 'FontSize', 8.5, 'FontWeight', 'bold');
    if vv == 1
        ylabel(ax, 'Pct pts from trend', 'FontSize', 7.5);
        xlabel(ax, 'Months after shock', 'FontSize', 7.5);
    end
end

% ── Panel Inferior: Standard Signs vs. Augmented Signs ────────────────────
for vv = 1:6
    ax = subplot(2, 6, vv + 6);

    % Standard signs (rojo/magenta, banda)
    area_x = [x_axis, fliplr(x_axis)];
    area_y_s = [irf_std_exo_lo(:, vv)', fliplr(irf_std_exo_hi(:, vv)')];
    fill(ax, area_x, area_y_s, col_std, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    hold(ax, 'on');
    plot(ax, x_axis, irf_std_exo_med(:, vv), 'Color', col_std, 'LineWidth', 1.5);

    % Augmented signs (verde, punteado)
    area_y_a = [irf_aug_exo_lo(:, vv)', fliplr(irf_aug_exo_hi(:, vv)')];
    fill(ax, area_x, area_y_a, col_aug, 'FaceAlpha', 0.15, 'EdgeColor', col_aug, 'LineStyle', ':', 'LineWidth', 0.5);
    plot(ax, x_axis, irf_aug_exo_med(:, vv), ':', 'Color', col_aug, 'LineWidth', 1.8);

    yline(ax, 0, 'k:', 'LineWidth', 0.8);
    xlim(ax, [0, h]);
    ax.XTick = 0:4:h;
    ax.FontSize = 7.5;
    box(ax, 'on');

    title(ax, var_labels{vv}, 'FontSize', 8.5, 'FontWeight', 'bold');
    if vv == 1
        ylabel(ax, 'Pct pts from trend', 'FontSize', 7.5);
        xlabel(ax, 'Months after shock', 'FontSize', 7.5);
    end
end

% Anotaciones
annotation(fig9, 'textbox', [0.01, 0.52, 0.18, 0.04], ...
    'String', '(Exclusion VS Augmented Signs)', 'FontSize', 8, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');
annotation(fig9, 'textbox', [0.01, 0.03, 0.18, 0.04], ...
    'String', '(Standard Signs VS Augmented Signs)', 'FontSize', 8, ...
    'FontWeight', 'bold', 'EdgeColor', 'none');

annotation(fig9, 'textbox', [0.68, 0.96, 0.32, 0.04], ...
    'String', ...
    'Blue=Exclusion/Standard  |  Green dotted=Augmented Signs', ...
    'FontSize', 7.5, 'EdgeColor', 'none');

% sgtitle eliminado — título se pone desde el .tex via \caption{}

% Exportar
exportgraphics(fig9, fullfile(fig_dir, ['fig9_comparison' fig_suffix '.pdf']), 'Resolution', 300);
exportgraphics(fig9, fullfile(fig_dir, ['fig9_comparison' fig_suffix '.png']), 'Resolution', 200);
fprintf('→ Figura 9 guardada: %s\n', fullfile(fig_dir, ['fig9_comparison' fig_suffix '.pdf']));

%% ============================================================================
%  SECCIÓN 7: RESUMEN NUMÉRICO
%% ============================================================================
fprintf('\n===== SECCIÓN 7: RESUMEN DE RESULTADOS =====\n');
fprintf('\nImpacto del choque monetario contractivo (h=0):\n');
fprintf('%-20s  %-12s  %-12s  %-12s\n', 'Variable', 'Exclusion', 'Std Signs', 'Aug Signs');
fprintf('%s\n', repmat('-', 1, 60));

for vv = 1:n
    fprintf('%-20s  %+8.4f      %+8.4f      %+8.4f\n', ...
        var_labels{vv}, ...
        irf_chol_exo_med(1, vv), ...
        irf_std_exo_med(1, vv), ...
        irf_aug_exo_med(1, vv));
end

fprintf('\n(+) = depreciación/aumento, (-) = apreciación/caída\n');
fprintf('\nβ_lower utilizado: %.4f (paper: -0.7)\n', beta_lower);
fprintf('\nTodas las figuras guardadas en: %s\n', fig_dir);
fprintf('Script completado.\n');

%% ============================================================================
%  FUNCIONES LOCALES
%% ============================================================================

function nll = beta_lower_nll(p, dpi, X_reg)
    % NLL condicional para Δπ = X*β + ω, ω ~ ARMA(1,2)
    % p = [β(3×1); φ1; θ1; θ2; log(σ²)]
    % Inicialización: eps(1)=eps(2)=0 (condicional)
    beta   = p(1:3);
    phi1   = p(4);
    theta1 = p(5);
    theta2 = p(6);
    sigma2 = exp(p(7));

    T     = length(dpi);
    omega = dpi - X_reg * beta;   % residuos de la regresión

    % Construir innovaciones recursivamente
    eps = zeros(T, 1);
    for t = 3:T
        eps(t) = omega(t) - phi1*omega(t-1) ...
                           - theta1*eps(t-1) ...
                           - theta2*eps(t-2);
    end

    % Verosimilitud condicional (descarta primeras 2 obs)
    eps_use = eps(3:end);
    T_eff   = length(eps_use);
    nll     = T_eff/2 * log(2*pi*sigma2) + sum(eps_use.^2) / (2*sigma2);
end
