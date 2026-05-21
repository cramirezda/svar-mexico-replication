%% check_beta_lower.m
%  Script de diagnóstico: solo corre Sección 1 (carga) + Sección 3 (β_lower)
%  de SVAR_SignRestrictions.m para verificar que OLS+NW da β_lower ≈ -0.700
%  Sin bootstrap → termina en segundos.

clear; clc;
rng(42);

base_dir = 'C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras';
data_dir  = fullfile(base_dir, 'raw data');

%% ── Carga de datos ───────────────────────────────────────────────────────────
opts = detectImportOptions(fullfile(data_dir, 'database_for_matlab.csv'));
opts = setvartype(opts, opts.VariableNames, 'double');
data_raw = readtable(fullfile(data_dir, 'database_for_matlab.csv'), 'ReadVariableNames', true);

endo_names = {'igae','pi_core_gap','pi_inpp_gap','dep_tcr','i_gap','m2_gap'};
exo_names  = {'y_us_gap','pi_us','pi_crb'};

Y = table2array(data_raw(:, endo_names));
Z = table2array(data_raw(:, exo_names));
[T, n] = size(Y);
fprintf('Muestra: T=%d, n=%d\n', T, n);

%% ── Sección 3: β_lower (OLS + Newey-West) ───────────────────────────────────
fprintf('\n===== CHECK β_lower =====\n');

dpi    = diff(Y(:, 2));
di_vec = diff(Y(:, 5));
D_t    = double(dpi .* di_vec < 0);
di_dt  = D_t .* di_vec;
T_reg  = length(dpi);
X_reg  = [ones(T_reg, 1), di_vec, di_dt];

% OLS
b_ols  = (X_reg' * X_reg) \ (X_reg' * dpi);
omega0 = dpi - X_reg * b_ols;
beta_minus_ols = b_ols(2) + b_ols(3);

% Newey-West HAC lag=3
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
Var_nw   = (XX \ S) / XX;
se_nw    = sqrt(max(0, Var_nw(2,2) + Var_nw(3,3) + 2*Var_nw(2,3)));
bl_nw    = beta_minus_ols - 1.645 * se_nw;

fprintf('  OLS+NW (primario):\n');
fprintf('    beta_minus = %.4f   (paper: -0.035)\n', beta_minus_ols);
fprintf('    SE (NW)    = %.4f   (paper implícito: ~0.404)\n', se_nw);
fprintf('    beta_lower = %.4f   (paper: -0.700)\n', bl_nw);

% MLE ARMA(1,2) comparación
phi1_init = sum(omega0(2:end) .* omega0(1:end-1)) / sum(omega0(1:end-1).^2);
phi1_init = max(-0.9, min(0.9, phi1_init));
init_p = [b_ols; phi1_init; 0.10; 0.05; log(var(omega0))];

nll_fn = @(p) beta_lower_nll(p, dpi, X_reg);
opts_mle = optimoptions('fminunc','Display','off','Algorithm','quasi-newton',...
    'MaxIterations',5000,'OptimalityTolerance',1e-9,'StepTolerance',1e-10);
[p_hat, ~, exitflag] = fminunc(nll_fn, init_p, opts_mle);

beta_minus_mle = p_hat(2) + p_hat(3);
h_fd = 1e-5; np = length(p_hat); H = zeros(np); f0 = nll_fn(p_hat);
for ii=1:np
    for jj=ii:np
        ei=zeros(np,1); ei(ii)=h_fd; ej=zeros(np,1); ej(jj)=h_fd;
        if ii==jj
            H(ii,ii)=(nll_fn(p_hat+ei)-2*f0+nll_fn(p_hat-ei))/h_fd^2;
        else
            H(ii,jj)=(nll_fn(p_hat+ei+ej)-nll_fn(p_hat+ei-ej)-...
                      nll_fn(p_hat-ei+ej)+nll_fn(p_hat-ei-ej))/(4*h_fd^2);
            H(jj,ii)=H(ii,jj);
        end
    end
end
Var_mle = inv(H);
se_mle  = sqrt(max(0, Var_mle(2,2) + Var_mle(3,3) + 2*Var_mle(2,3)));
bl_mle  = beta_minus_mle - 1.645*se_mle;

fprintf('\n  MLE ARMA(1,2) (comparación):\n');
fprintf('    beta_minus = %.4f   exitflag=%d\n', beta_minus_mle, exitflag);
fprintf('    SE (Hess)  = %.4f\n', se_mle);
fprintf('    beta_lower = %.4f\n', bl_mle);

fprintf('\n  → beta_lower USADO (OLS+NW): %.4f\n', bl_nw);
fprintf('  ✓ Match paper: %s\n', string(abs(bl_nw - (-0.7)) < 0.15));

%% ── Función local ────────────────────────────────────────────────────────────
function nll = beta_lower_nll(p, dpi, X_reg)
    beta   = p(1:3);
    phi1   = p(4); theta1 = p(5); theta2 = p(6);
    sigma2 = exp(p(7));
    T = length(dpi);
    omega = dpi - X_reg * beta;
    eps = zeros(T, 1);
    for t = 3:T
        eps(t) = omega(t) - phi1*omega(t-1) - theta1*eps(t-1) - theta2*eps(t-2);
    end
    eps_use = eps(3:end);
    T_eff   = length(eps_use);
    nll = T_eff/2*log(2*pi*sigma2) + sum(eps_use.^2)/(2*sigma2);
end
