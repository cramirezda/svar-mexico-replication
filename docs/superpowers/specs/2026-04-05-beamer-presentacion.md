# Diseño: Presentación Beamer — Replicación SVAR México

**Fecha:** 2026-04-05  
**Proyecto:** Replicación Carrillo & Elizondo (2015/2018) — México  
**Estado:** Aprobado por usuario

---

## Objetivo

Generar un archivo LaTeX compilable con presentación Beamer (~15 slides) en español para seminario interno. La audiencia conoce el paper de Carrillo & Elizondo (2015/2018); el énfasis está en las diferencias entre la replicación y el paper original, con código MATLAB ilustrando las estrategias de identificación — especialmente el algoritmo de restricciones de signo aumentadas.

---

## Especificaciones de formato

- **Tema:** Metropolis (`\usetheme{metropolis}`)
- **Idioma:** Español
- **Longitud:** ~15 diapositivas
- **Audiencia:** Seminario interno, conocen C&E
- **Salida:** `presentacion_replicacion.tex` en la raíz del proyecto

---

## Entradas (imágenes a incluir)

Todas las imágenes se incluyen con `\includegraphics`. Las siguientes deben existir en `figures/`:

| Archivo | Descripción |
|---|---|
| `figures/fig_series_originales.png` | Series en niveles con tendencia HP |
| `figures/fig_series_gaps.png` | Series transformadas (gaps/tasas) |
| `figures/varstable.png` | Raíces del VAR (STATA) |
| `figures/fig7_exclusion.png` | IRFs restricciones de exclusión (MATLAB) |
| `figures/fig8_sign_restrictions.png` | IRFs restricciones de signo (MATLAB) |
| `figures/fig9_comparison.png` | Comparación 3 estrategias (MATLAB) |
| `figures/tabla_adf_pp.png` | Tabla ADF/PP |
| `figures/tabla_varsoc.png` | Criterios selección de rezagos |
| `figures/tabla_lm_norm.png` | LM autocorrelación + Jarque-Bera |
| `figures/tabla_descriptivos.png` | Estadísticos descriptivos |

---

## Estructura de diapositivas (15 slides)

### §0 — Portada + Agenda (2 slides)

**Slide 1 — Portada:**
- Título: "Replicación de Carrillo & Elizondo (2015/2018): México"
- Subtítulo: "SVAR con Restricciones de Exclusión y de Signo"
- Autor, institución, fecha
- Período muestral: Enero 2001 – Marzo 2014, T=159

**Slide 2 — Agenda:**
4 bloques con íconos o numeración:
1. Procesamiento de datos
2. Diagnósticos estadísticos
3. Estrategias de identificación
4. Diferencias y conclusiones

---

### §1 — Cambios en procesamiento de datos (3 slides)

**Slide 3 — Tabla de cambios:**

| Dimensión | Paper C&E | Replicación |
|---|---|---|
| Período muestral | 2001:01–2014:03 (T=159) | 2001:01–2014:03 (T=159) |
| TCR bilateral | 21 socios comerciales | 49 socios comerciales (Banxico) |
| Producto doméstico | Brecha IGAE precalculada | Índice IGAE nivel (INEGI) |
| Ajuste estacional | No especificado | X-13ARIMA-SEATS (R `seasonal`) |
| Filtro de tendencia | HP λ=129,600 | HP λ=129,600 |
| β_lower | −0.7 | −0.715 |

**Slide 4 — Series originales:**
- `fig_series_originales.png` (grilla 3×3, series en niveles con tendencia HP)

**Slide 5 — Series transformadas:**
- `fig_series_gaps.png` (grilla 3×3, gaps y tasas anualizadas)

---

### §2 — Diagnósticos estadísticos (2 slides)

**Slide 6 — Pruebas de raíz unitaria + rezagos + estabilidad:**
- Columna izquierda: `tabla_adf_pp.png` (ADF/PP, todos I(0))
- Columna derecha: `tabla_varsoc.png` (AIC/SIC/LR → p=3) + `varstable.png`

**Slide 7 — Autocorrelación + β_lower:**
- `tabla_lm_norm.png` (LM residuos, Jarque-Bera por ecuación)
- Ecuación 11 del paper: $\beta^- = -0.400$, SE$=0.192$, $\beta_{lower} = -0.715$ (paper: $-0.7$)

---

### §3 — Estrategias de identificación (5 slides)

**Slide 8 — Estrategia 1: Restricciones de Exclusión:**
- `fig7_exclusion.png` (IRFs Cholesky con bandas bootstrap)
- Tabla impactos h=0: paper vs. replicación para las 6 variables endógenas

**Slide 9 — Bootstrap Runkle (1987) — Código MATLAB:**
Chunk de código de la función `runkle_bootstrap_data` del archivo `SVAR_SignRestrictions.m`:
```matlab
% Runkle (1987) recursive bootstrap
function [lower, upper] = runkle_bootstrap_data(...)
    % Re-estima B_0 en cada draw
    ...
end
```
Descripción textual: resampling de residuos, re-estimación de B₀ en cada réplica, percentiles 16/84.

**Slide 10 — Estrategia 2: Restricciones de Signo Estándar:**
- `fig8_sign_restrictions.png`
- Modelos aceptados: 198,654 (std) / 66,874 (aug)
- Método: matrices ortogonales aleatorias Q ~ Haar, aceptación por signo en h=0..2

**Slide 11 — Algoritmo Aumentado: Descripción (Algorithm 1):**
Pseudocódigo paso a paso espejando Algorithm 1 del paper:
1. Estimar VAR reducido → obtener Σ_u, β̂
2. Calcular β_lower (Ec. 11)
3. Para cada draw i = 1...S:
   a. Generar Q ~ Haar(n×n)
   b. Calcular B = chol(Σ_u) · Q
   c. Verificar restricciones de signo estándar en columna del shock monetario
   d. Si pasa: calcular elasticidad β_draw
   e. Si β_draw < β_lower: rechazar (restricción aumentada)
   f. Si pasa todo: guardar IRF_i
4. Reportar mediana e intervalo [16,84] sobre IRF aceptadas

**Slide 12 — Algoritmo Aumentado: Código MATLAB:**
Chunk central de `sign_restricted_svar` con las líneas clave:
- Generación de Q con descomposición QR de matriz normal aleatoria
- Evaluación de restricciones de signo
- Evaluación de restricción aumentada (β < β_lower)
- Acumulación de IRFs aceptadas

**Slide 13 — Comparación Final (Figura 9):**
- `fig9_comparison.png`
- Tabla: 3 estrategias × 6 variables, impactos h=0 (paper vs. replicación)

---

### §4 — Diferencias y conclusiones (2 slides)

**Slide 14 — Tabla: Diferencias principales:**

| Dimensión | Paper | Replicación | Causa probable |
|---|---|---|---|
| T observaciones | 159 | 159 | ✓ Resuelto (filtro extendido a dic-2000) |
| TCR | 21 socios | 49 socios | Datos disponibles Banxico |
| β_lower | −0.7 | −0.715 | Diferencia marginal, misma dirección |
| IRFs Cholesky | [referencia visual] | Muy similares | Bootstrap full Runkle implementado |
| Modelos sign rest. | N/D | 198,654 aceptados | Sin referencia en el paper |
| Producto doméstico | Brecha precalculada | Índice IGAE nivel | Fuente INEGI directo |

**Slide 15 — Conclusiones:**
- La replicación recupera los resultados cualitativos centrales del paper
- Price puzzle eliminado con restricciones de signo (las tres estrategias coinciden)
- Diferencia cuantitativa más relevante: TCR (49 vs 21 socios) → sin impacto cualitativo en IRFs
- β_lower prácticamente idéntico: −0.715 vs −0.7
- Restricciones de signo aumentadas operan como filtro efectivo sobre distribución de modelos aceptados

---

## Paquetes LaTeX requeridos

```latex
\usepackage[spanish]{babel}
\usepackage[utf8]{inputenc}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage{listings}   % para chunks de código MATLAB
\usepackage{xcolor}
\usepackage{amsmath}
\usepackage{colortbl}
\usepackage{multirow}
```

---

## Archivos que se crean/modifican

| Archivo | Acción |
|---|---|
| `presentacion_replicacion.tex` | CREAR — presentación Beamer compilable |
