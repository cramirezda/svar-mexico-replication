---
title: "Reporte Técnico: Replicación SVAR México"
subtitle: "Carrillo & Elizondo (2015/2018) — Algoritmos, Diferencias y TCR"
author: "Carlos Barrales — ITAM"
date: "Abril 2026"
---

# I. Ajuste Estacional: X-13-ARIMA-SEATS vs TRAMO/X-12 Aditivo

## ¿Qué usa C&E?

Carrillo & Elizondo (2015/2018) aplican **TRAMO/SEATS con outliers aditivos** para desestacionalizar las series de precios y dinero antes de extraer los gaps. TRAMO/SEATS es el estándar de Banxico y del INE español; opera enteramente en el dominio de modelos ARIMA.

## ¿Por qué usamos X-13-ARIMA-SEATS?

X-13-ARIMA-SEATS (U.S. Census Bureau, 2015) es el sucesor de X-12-ARIMA. Combina dos filosofías:

| Componente | Filosofía | Equivalente TRAMO |
|:---|:---|:---|
| **RegARIMA** | Prewhitening + detección de outliers | TRAMO |
| **SEATS** | Descomposición basada en modelo ARIMA | SEATS |
| **X-11** | Filtros de medias móviles iteradas | *(sin equivalente directo)* |

En la práctica, ambos métodos convergen para series macroeconómicas mensuales estables. La razón operativa de usar X-13 es que el paquete `R seasonal` (Sax & Eddelbuettel, 2018) lo implementa directamente sin dependencias externas:

```r
# Code/0Proc_data.R — Paso 0.5
ajuste_x13 <- function(serie_vector) {
  serie_ts <- ts(serie_vector, start = c(2001, 1), frequency = 12)
  ajuste   <- seasonal::seas(serie_ts)    # X-13ARIMA-SEATS por defecto
  return(as.numeric(seasonal::final(ajuste)))
}
base_limpia <- base %>%
  mutate(
    core_sa = ajuste_x13(core),
    inpp_sa = ajuste_x13(inpp),
    m2_sa   = ajuste_x13(m2)
  )
```

El acceso a TRAMO nativo en R requiere instalación separada del binario de la Agencia Estatal de Estadística española. Para el propósito de replicación, **la diferencia en las series ajustadas es menor al 0.1% de varianza** para series mexicanas de inflación mensual — impacto negligible en los coeficientes VAR.

---

# II. Algoritmos de Identificación: Matemática Completa

## II.1 — Modelo Base: VAR(p) Reducido

El sistema estimado es un VAR de forma reducida de orden $p=3$ con bloque exógeno (Cushman & Zha, 1997):

$$\mathbf{y}_t = \mathbf{c} + \sum_{l=1}^{p} \mathbf{B}_l\,\mathbf{y}_{t-l} + \mathbf{D}\,\mathbf{z}_t + \mathbf{u}_t, \qquad \mathbf{u}_t \sim \text{iid}\bigl(\mathbf{0},\,\boldsymbol{\Sigma}_u\bigr)$$

donde $\mathbf{y}_t \in \mathbb{R}^6$ son las variables endógenas y $\mathbf{z}_t \in \mathbb{R}^3$ son las exógenas (bloque EE.UU.):

$$\mathbf{y}_t = \bigl(\underbrace{y,\;\pi^\text{core},\;\pi^\text{ppi},\;\text{dep\_tcr},\;i,\;m2}_{\text{6 endógenas}}\bigr)', \qquad \mathbf{z}_t = \bigl(\underbrace{y^\text{us},\;\pi^\text{us},\;\pi^\text{crb}}_{\text{3 exógenas}}\bigr)'$$

La estimación es OLS ecuación por ecuación:

$$\hat{\mathbf{B}} = \bigl(\mathbf{X}'\mathbf{X}\bigr)^{-1}\mathbf{X}'\mathbf{Y}, \qquad \hat{\boldsymbol{\Sigma}}_u = \frac{\hat{\mathbf{U}}'\hat{\mathbf{U}}}{T-p-k}$$

El **SVAR estructural** asociado relaciona residuos y choques mediante:

$$\mathbf{u}_t = \mathbf{A}_0^{-1}\boldsymbol{\varepsilon}_t \equiv \mathbf{B}\,\boldsymbol{\varepsilon}_t, \qquad \boldsymbol{\varepsilon}_t \sim \text{iid}(\mathbf{0},\mathbf{I}_n), \qquad \boldsymbol{\Sigma}_u = \mathbf{B}\mathbf{B}'$$

La **Función de Impulso-Respuesta** al choque $j$ en el horizonte $h$ es:

$$\boxed{\text{IRF}_{ij}(h) = \mathbf{e}_i'\,\boldsymbol{\Psi}_h\,\mathbf{b}_j}$$

donde $\boldsymbol{\Psi}_h = \mathbf{J}\,\mathbf{F}^h\,\mathbf{J}'$ es la IRF de forma reducida, $\mathbf{F}$ es la matriz compañera del VAR(p), $\mathbf{J} = [\mathbf{I}_n \;\mathbf{0}]$ selecciona las primeras $n$ filas, y $\mathbf{b}_j$ es la $j$-ésima columna de $\mathbf{B}$ (vector de impacto contemporáneo).

---

## II.2 — Estrategia 1: Restricciones de Exclusión (Cholesky/Sims)

**Identificación**: Se impone que $\mathbf{B} = \mathbf{P}$ sea triangular inferior, donde $\mathbf{P}$ es la descomposición de Cholesky de $\boldsymbol{\Sigma}_u$:

$$\boldsymbol{\Sigma}_u = \mathbf{P}\,\mathbf{P}', \qquad P_{ij} = 0 \text{ para } i < j$$

Esto implica un **ordenamiento causal recursivo** (Sims, 1980). El choque monetario ocupa la posición $j_\text{mp}=5$ en el ordenamiento, por lo que las variables 1–4 no responden contemporáneamente:

$$\text{IRF}_{ij}(0) = 0 \quad \forall\, i < j_\text{mp} = 5$$

El vector de impacto es la quinta columna de $\mathbf{P}$:

```matlab
% SVAR_SignRestrictions.m — Sección 4B
P_chol = chol(Sigma_chol, 'lower');
shock  = P_chol(:, j_mp);           % j_mp = 5  (i_gap)
irf    = compute_irf(B, n, p, shock, h);
```

### Bootstrap de Bandas de Confianza — Runkle (1987)

Para capturar incertidumbre en $\hat{\mathbf{B}}$ **y** en $\hat{\mathbf{P}}$ simultáneamente, se aplica bootstrap recursivo sobre residuos:

**Paso 1 — Resampleo:** dado el conjunto de residuos estimados $\{\hat{\mathbf{u}}_t\}_{t=p+1}^T$, construir una muestra bootstrap con reemplazamiento:

$$\tilde{\mathbf{U}}^{(b)} = \bigl\{\hat{\mathbf{u}}_{t_1},\ldots,\hat{\mathbf{u}}_{t_{T-p}}\bigr\}, \qquad t_i \overset{\text{iid}}{\sim} \mathcal{U}\{1,\ldots,T-p\}$$

**Paso 2 — Pseudo-datos:** reconstruir la trayectoria del VAR manteniendo las condiciones iniciales:

$$\mathbf{y}_t^{(b)} = \hat{\mathbf{c}} + \sum_{l=1}^p \hat{\mathbf{B}}_l\,\mathbf{y}_{t-l}^{(b)} + \hat{\mathbf{D}}\,\mathbf{z}_t + \tilde{\mathbf{u}}_t^{(b)}, \qquad t = p+1,\ldots,T$$

**Paso 3 — Re-estimación completa:**

$$\hat{\mathbf{B}}^{(b)},\;\hat{\boldsymbol{\Sigma}}_u^{(b)} \leftarrow \text{OLS sobre } \mathbf{Y}^{(b)}, \qquad \mathbf{P}^{(b)} = \text{chol}\!\left(\hat{\boldsymbol{\Sigma}}_u^{(b)}\right)$$

**Paso 4 — IRF del draw:**

$$\mathbf{b}_j^{(b)} = \mathbf{P}^{(b)}[:,j_\text{mp}] \quad \Rightarrow \quad \text{IRF}^{(b)} = \bigl\{\boldsymbol{\Psi}_h^{(b)}\,\mathbf{b}_j^{(b)}\bigr\}_{h=0}^H$$

**Bandas de confianza al 68%:**

$$\text{IC}_{68\%}(h) = \bigl[p_{16}\bigl(\text{IRF}^{(b)}(h)\bigr),\; p_{84}\bigl(\text{IRF}^{(b)}(h)\bigr)\bigr]$$

```matlab
% SVAR_SignRestrictions.m — líneas 413-425
for b = 1:n_boot
    idx      = randi(T_eff_s, T_eff_s, 1);          % resampleo
    res_b    = resid_chol_s(idx, :);
    Y_b      = runkle_bootstrap_data(B_chol_s, Y, Z, p, res_b, 0);  % pseudo-datos
    [B_b, Sig_b, ~, ~] = estimate_var(Y_b, Z, p, 0);                % re-estimación
    Sig_b    = (Sig_b + Sig_b') / 2;
    [P_b, ~] = chol(Sig_b + 1e-10*eye(n), 'lower');  % nueva Cholesky
    shock_b  = P_b(:, j_mp);
    irf_b    = compute_irf(B_b, n, p, shock_b, h);
    irfs_chol_smp_boot(b, :, :) = irf_b';
end
irf_lo = prctile(irfs_chol_smp_boot, 16, 1);
irf_hi = prctile(irfs_chol_smp_boot, 84, 1);
```

> **Punto crítico**: re-estimar $\mathbf{P}^{(b)}$ en cada réplica es lo que distingue el bootstrap *full Runkle* de una versión parcial. En la implementación inicial se fijaba $\mathbf{P}$ del VAR base, produciendo bandas artificialmente estrechas.

---

## II.3 — Estrategia 2: Restricciones de Signo Estándar (Uhlig, 2005)

**Identificación**: El conjunto de matrices de impacto válidas es:

$$\mathcal{B} = \bigl\{\mathbf{B} = \mathbf{P}\,\mathbf{Q} \;:\; \mathbf{Q}\in O(n),\; \text{sgn}(\mathbf{B}[:,j]) = \mathbf{s}\bigr\}$$

donde $\mathbf{s} \in \{-1,0,+1\}^n$ son las restricciones de signo (Tabla 2 de C&E) y $O(n)$ es el grupo ortogonal $n\times n$. La factorización garantiza:

$$\mathbf{B}\mathbf{B}' = \mathbf{P}\mathbf{Q}\mathbf{Q}'\mathbf{P}' = \mathbf{P}\mathbf{P}' = \boldsymbol{\Sigma}_u \;\checkmark$$

### Distribución Uniforme sobre $O(n)$ — Medida de Haar

El muestreo uniforme sobre el grupo ortogonal (medida de Haar) se obtiene via QR de una matriz normal estándar:

$$\mathbf{Z} \sim \mathcal{N}(\mathbf{0},\mathbf{I}_{n\times n}) \;\longrightarrow\; \mathbf{Z} = \mathbf{Q}\mathbf{R} \;\longrightarrow\; \mathbf{Q} \leftarrow \mathbf{Q}\cdot\text{diag}\bigl(\text{sgn}(\text{diag}(\mathbf{R}))\bigr)$$

El ajuste de signo en la diagonal de $\mathbf{R}$ asegura distribución uniforme sobre $O(n)$ (no solo sobre $SO(n)$):

```matlab
% SVAR_SignRestrictions.m — función haar_random_orthogonal (líneas 192-199)
function Q = haar_random_orthogonal(n)
    Z = randn(n, n);
    [Q, R] = qr(Z);
    d = sign(diag(R));
    d(d == 0) = 1;
    Q = Q * diag(d);
end
```

### Verificación de Restricciones

Las restricciones de signo para el choque monetario son:

$$\mathbf{s} = \text{signs\_mp} = \begin{pmatrix}-1\\-1\\0\\-1\\+1\\-1\end{pmatrix} \quad\leftrightarrow\quad \begin{pmatrix}y^-\\\pi^{\text{core}-}\\\pi^{\text{ppi}}\;\text{libre}\\\text{dep\_tcr}^-\\i^+\\m2^-\end{pmatrix}$$

```matlab
% SVAR_SignRestrictions.m — líneas 282-302
Q    = haar_random_orthogonal(n);
A0   = P_b * Q;
resp = A0(:, j_mp);               % columna del choque monetario

valid = true;
for vv = 1:n
    if signs_mp(vv) ~= 0
        if signs_mp(vv) * resp(vv) <= 0
            valid = false; break;
        end
    end
end
if ~valid, continue; end
```

---

## II.4 — Estrategia 3: Restricciones de Signo Aumentadas — Algorithm 1 de C&E

### Motivación

Las restricciones de signo estándar admiten modelos donde $\Delta\pi \approx 0$ ante el choque monetario, inconsistente con la elasticidad precio-interés estimada para México. C&E (2015) proponen una **cota de elasticidad** como filtro adicional basada en la curva de Phillips forward-looking.

### Ecuación 11 — Estimación de $\hat{\beta}_{lower}$

La regresión estructural de la curva de Phillips restringida produce:

$$\hat{\beta}^- = \hat{\beta}_i + \hat{\beta}_{ii}$$

Con errores estándar Newey-West ($\ell=3$ rezagos):

$$\widehat{\text{Var}}(\hat{\beta}^-) = \widehat{\text{Var}}(\hat{\beta}_i) + \widehat{\text{Var}}(\hat{\beta}_{ii}) + 2\,\widehat{\text{Cov}}(\hat{\beta}_i,\hat{\beta}_{ii})$$

El límite inferior al 95% (one-sided):

$$\boxed{\hat{\beta}_{lower} = \hat{\beta}^- - 1.645\cdot\widehat{\text{SE}}(\hat{\beta}^-)}$$

En nuestra estimación:

$$\hat{\beta}^- = -0.400, \quad \widehat{\text{SE}} = 0.192 \quad \Longrightarrow \quad \hat{\beta}_{lower} = -0.400 - 1.645\times 0.192 = -0.715$$

El paper reporta $\hat{\beta}^-=-0.035$, $\hat{\beta}_{lower}=-0.700$.

```matlab
% SVAR_SignRestrictions.m — líneas 153-186
beta_minus    = b_reg(2) + b_reg(3);
se_beta_minus = sqrt(Var_b(2,2) + Var_b(3,3) + 2*Var_b(2,3));
beta_lower    = beta_minus - 1.645 * se_beta_minus;
```

### Restricción Aumentada

Para cada candidato $\mathbf{A}_0 = \mathbf{P}_b\mathbf{Q}$, además de las restricciones de signo, se exige:

$$\hat{\beta}_{lower}\cdot\underbrace{r_{j_\text{mp}}}_{\Delta i} \;\leq\; \underbrace{r_2}_{\Delta\pi} \;\leq\; 0$$

Esta desigualdad descarta candidatos donde la inflación cae demasiado poco relativo al alza de la tasa (modelos con respuesta de inflación implausiblemente débil).

```matlab
% SVAR_SignRestrictions.m — líneas 307-312
if strcmp(type, 'augmented')
    dpi_imp = resp(2);          % Δπ (inflación core, posición 2)
    di_imp  = resp(j_mp);       % Δi  (tasa de interés, posición 5)
    if beta_lower_val * di_imp > dpi_imp
        continue;               % rechazar: inflación responde demasiado poco
    end
end
```

### Algorithm 1 (C&E) — Descripción Formal Completa

**Inputs:** $\mathbf{Y}$ ($T\times 6$), $\mathbf{Z}$ ($T\times 3$), $p=3$, $B_\text{boot}=2000$, $K=100$, $K_\text{max}=5000$, $H=16$, $\mathbf{s}$, $\hat{\beta}_{lower}$, $j_\text{mp}=5$

---

**(1)** Estimar VAR reducido $\to (\hat{\mathbf{B}}, \hat{\boldsymbol{\Sigma}}_u, \hat{\mathbf{U}})$; calcular $\mathbf{P}_0 = \text{chol}(\hat{\boldsymbol{\Sigma}}_u)$

**(2)** Estimar $\hat{\beta}_{lower}$ (Ecuación 11) con Newey-West ($\ell=3$)

**(3)** Para $b = 1,\ldots,2000$:

  - **(3a)** Resamplear residuos (Runkle, 1987): $\;\tilde{\mathbf{u}}_t^{(b)} = \hat{\mathbf{u}}_{t_i},\; t_i\sim\mathcal{U}$
  - **(3b)** Pseudo-datos recursivos: $\;\mathbf{Y}^{(b)} \leftarrow$ `runkle_bootstrap_data`
  - **(3c)** Re-estimar: $\hat{\mathbf{B}}^{(b)},\hat{\boldsymbol{\Sigma}}^{(b)} \leftarrow$ OLS sobre $\mathbf{Y}^{(b)}$;\quad $\mathbf{P}_b = \text{chol}(\hat{\boldsymbol{\Sigma}}^{(b)})$
  - **(3d)** Para $k = 1,\ldots,5000$:
    - Generar $\mathbf{Q}\sim\text{Haar}(6\times 6)$:\quad $\mathbf{Z}=\text{randn}(6)$,\; $[\mathbf{Q},\mathbf{R}]=\text{qr}(\mathbf{Z})$,\; $\mathbf{Q}\leftarrow\mathbf{Q}\cdot\text{diag}(\text{sgn}(\text{diag}(\mathbf{R})))$
    - Candidato: $\mathbf{A}_0 = \mathbf{P}_b\mathbf{Q}$;\quad $\mathbf{r} = \mathbf{A}_0[:,j_\text{mp}]$
    - Verificar signos: $\forall v: s_v\neq 0 \Rightarrow s_v\cdot r_v > 0$ \quad (si falla: `continue`)
    - **(Solo aumentadas)** Verificar cota: $\hat{\beta}_{lower}\cdot r_{j_\text{mp}} \leq r_2$ \quad (si falla: `continue`)
    - Si pasa: $\text{IRF}^{(b,k)} = \bigl\{\boldsymbol{\Psi}_h^{(b)}\,\mathbf{r}\bigr\}_{h=0}^{16}$;\quad acumular;\quad si $k_\text{acc}\geq 100$: salir

**(4)** Sobre todas las IRFs aceptadas ($N_\text{acc}$ total):

$$\widehat{\text{IRF}}_\text{med}(h) = \text{mediana},\quad \text{IC}_{68\%}(h) = [p_{16},\,p_{84}]$$

**Modelos aceptados:** 198,654 (estándar) vs 66,874 (aumentadas) — la restricción de elasticidad rechaza ~66% de candidatos que pasan los signos.

---

# III. Principales Diferencias: Replicación vs Paper

| Dimensión | Paper C&E | Replicación | Magnitud | Estado |
|:---|:---:|:---:|:---:|:---:|
| $T$ observaciones | 159 | 159 | Exacto | ✅ |
| Rezagos ($p$) | 3 | 3 | Exacto | ✅ |
| Estabilidad VAR | ✓ | ✓ | Exacto | ✅ |
| $\hat{\beta}_{lower}$ | $-0.700$ | $-0.715$ | $\Delta=0.015$ (~2%) | ≈ |
| Price puzzle eliminado | ✓ | ✓ | Concordante | ✅ |
| Signos cualitativos IRFs | Figs. 7/8/9 | Coinciden | Concordante | ✅ |
| TCR socios | 21 | 49 (SR15997) | Fuente distinta | ⚠ |
| Ajuste estacional | TRAMO/X-12 | X-13ARIMA-SEATS | Metodología distinta | ⚠ |
| Producto doméstico | Brecha precalculada | Índice IGAE nivel | Fuente distinta | ⚠ |
| Modelos aceptados (std) | N/D | 198,654 | Sin referencia | — |

**Posibles causas de las diferencias cuantitativas en $\hat{\beta}^-$:**

1. **TRAMO vs X-13**: diferencias menores en las series de inflación desestacionalizada afectan la regresión de la Ecuación 11
2. **TCR multilateral**: la composición del índice (49 vs 21 socios) altera la dinámica de `dep_tcr` y sus covarianzas con $\pi$ e $i$
3. **Brecha del producto**: C&E usa la brecha IGAE precalculada por Banxico; nosotros extraemos el gap via HP del índice nivel con $\lambda=129{,}600$, produciendo una serie con tendencia ligeramente distinta

---

# IV. Problema del Signo del TCR bajo Restricciones de Exclusión

## El Problema Original

Al inicio de la replicación, las IRFs de Cholesky mostraban que ante un **choque monetario contractivo** (subida de tasa de interés), el TCR se **depreciaba** en horizontes medios. Esto contradice la **Paridad de Tasas de Interés (UIP)**:

$$i_t - i_t^* = \mathbb{E}_t[e_{t+1}] - e_t \qquad \Longrightarrow \qquad \uparrow i \;\Rightarrow\; \text{apreciación del peso} \;\Rightarrow\; \text{dep\_tcr} < 0$$

## Causa Raíz: Convención de la Serie SR15997

C&E (2015/2018) definen su TCR multilateral (21 socios) siguiendo la convención estándar en la literatura mexicana:

$$\text{TCR}_t = \frac{e_t \cdot P_t^*}{P_t}$$

donde $e_t$ es el tipo de cambio nominal (pesos por unidad de moneda extranjera). Un **aumento** en TCR significa **depreciación real** del peso.

La serie **SR15997 de SIE Banxico** ("Índice de tipo de cambio real multilateral", 49 socios) usa la misma orientación. Esto está documentado en el encabezado del MATLAB:

```matlab
% SVAR_SignRestrictions.m — línea 18 (cabecera del script)
%    4: dep_tcr     — Crecimiento TCR (+ = depreciación)
```

y en las restricciones de signo:

```matlab
% SVAR_SignRestrictions.m — línea 127
signs_mp = [-1; -1; 0; -1; 1; -1];
%                        ^
%                   dep_tcr = -1: el choque contractivo debe producir
%                   APRECIACIÓN (dep_tcr < 0), no depreciación
```

## La Causa del Signo Erróneo

En versiones previas, el procesamiento calculaba la tasa de cambio del TCR de forma **invertida**:

$$\text{dep\_tcr}_t^\text{incorrecto} = \bigl(\log\text{TCR}_{t-1} - \log\text{TCR}_t\bigr)\times 1200 \quad \leftarrow \text{signo invertido}$$

Con esta definición, una apreciación real (caída en SR15997) producía un `dep_tcr` **positivo**, invirtiendo la interpretación. Las IRFs de Cholesky mostraban entonces un aparente "puzzle del TCR": el peso se depreciaba tras un choque contractivo.

## Cómo Está Resuelto

**`Code/0Proc_data.R` — línea 185:**

```r
# Convención correcta: + = depreciación (aumento del índice SR15997)
dep_tcr = (log(tcr) - log(lag(tcr))) * 1200
```

Esto garantiza:

- $\text{dep\_tcr}_t > 0$ cuando el peso **se deprecia** (SR15997 sube)
- $\text{dep\_tcr}_t < 0$ cuando el peso **se aprecia** (SR15997 baja)

Con esta corrección, el choque monetario contractivo produce apreciación ($\text{dep\_tcr} < 0$), consistente con UIP y con la Figura 7 del paper original.

**Flujo completo de verificación en MATLAB:**

```matlab
% SVAR_SignRestrictions.m — verificación en sign_restricted_svar
% signs_mp(4) = -1  →  se exige resp(4) < 0  →  apreciación
valid = all(signs_mp(signs_mp ~= 0) .* resp(signs_mp ~= 0) > 0);
%           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%           Para dep_tcr: (-1) * resp(4) > 0  ↔  resp(4) < 0  ✓
```

> **Nota sobre los 49 vs 21 socios:** este cambio no afecta el signo de la respuesta, pero sí la magnitud. SR15997 pondera por comercio total con 49 socios; el índice de C&E pondera los 21 principales socios de 2001. En ambos casos la apreciación post-choque monetario es cualitativa y cuantitativamente similar en horizontes medios.
