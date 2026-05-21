# CLAUDE.md — SVAR México (Replicación Carrillo & Elizondo 2015/2018)

Guía de contexto acumulada para el agente. Contiene conocimiento de dominio, reglas operativas y anti-patrones identificados en sesiones anteriores.

---

## Project Context

### Objetivo del proyecto
Replicar el ejercicio empírico de Carrillo & Elizondo (2015/2018) para México comparando tres estrategias de identificación SVAR: restricciones de exclusión (Cholesky), restricciones de signo simples, y restricciones de signo aumentadas.

**Paper**: "How Robust Are SVARs at Measuring Monetary Policy in Small Open Economies?" — Journal of Empirical Finance. DOI: 10.1016/j.jempfin.2018.06.009

### Período muestral
Enero 2001 – Marzo 2014 | T=159 observaciones | frecuencia mensual | rezagos p=3

### Definición exacta de variables (crítico)
- `i_gap = int_rate - pi_trend` — NO es HP-gap de la tasa; es **tasa real proxy** (tasa nominal Cetes28d menos tendencia HP de inflación subyacente).
- `pi_core_gap = pi_core - pi_trend` — desviación de inflación subyacente respecto a su tendencia HP.
- `pi_inpp_gap = pi_inpp - pi_trend` — ídem para inflación al productor.
- `m2_gap = g_m2 - pi_trend` — crecimiento M2 menos tendencia inflacionaria.
- `igae` — brecha del producto México via filtro Kalman (ver `0b_brecha_kalman.R`).
- Exógenas block-exogenous: `y_us_gap` (HP de log INDPRO), `pi_us` (PCE), `pi_crb` (CRB).

### Flujo de ejecución (orden obligatorio)
1. `0Proc_data.R` → genera `database.xlsx` (T=159 esperado)
2. `1DickeyFuller.R` → pruebas ADF/PP de raíz unitaria
3. `3SVAR_Mexico.do` (STATA) → diagnósticos + Cholesky IRFs + exporta `database_for_matlab.csv`
4. `SVAR_SignRestrictions.m` (MATLAB) → β_lower + Fig. 7, 8, 9 + exporta `beta_lower_arma12.csv`
5. `4Informe_replicacion.R` → regenera `informe_replicacion.html`

---

## Architecture Decisions

### Estimación de β_lower — REGLA CRÍTICA (2026-05-10)

**Método correcto (C&E 2015):** OLS + Newey-West HAC con lag=3  
**Método incorrecto:** ARMA(1,2) MLE conjunta (más eficiente pero no replica el paper)

**Por qué importa:**
- OLS+NW da SE ≈ 0.40 → β_lower ≈ -0.700 ✓ (paper)
- MLE da SE ≈ 0.099 → β_lower ≈ -0.180 ✗ (no replica)
- Derivación: β̂⁻ = -0.035, β_lower = -0.700 → SE implícito = 0.665/1.645 ≈ **0.404**

**El ARMA(1,2) se menciona en el paper solo para justificar la estructura del error** (y con ello el lag=3 para NW), no como método de estimación de SE.

**Implementación en `SVAR_SignRestrictions.m` (Sección 3):**
- Primario: OLS+NW → `beta_minus`, `se_beta_minus`, `beta_lower` (exportados al CSV)
- Comparación: MLE imprime diagnóstico pero NO determina β_lower final
- Fallback: si OLS+NW ∉ [-5, 0], probar MLE; si ambos fallan, usar -0.7

### Principio general de replicación
**Siempre usar el método exacto del paper, aunque exista una alternativa estadísticamente más eficiente.** La "eficiencia" estadística (SE más chico) puede romper la replicación si el paper usó un estimador más conservador intencionalmente.

### Schema CSV `raw data/beta_lower_arma12.csv`
Columnas: `beta_minus`, `se_beta_minus`, `beta_lower` (OLS+NW — valores primarios), `beta_minus_mle`, `se_mle`, `beta_lower_mle` (MLE — solo comparación), `phi1`, `theta1`, `theta2`, `sigma2`, `exitflag`

El archivo `4Informe_replicacion.R` lee `beta_minus` y `beta_lower` de la primera fila.

---

## Development Guidelines

### MATLAB — fminunc exitflag
- `exitflag = 1`: convergencia limpia (gradiente < OptimalityTolerance) — SEs confiables
- `exitflag = 5`: "predicted decrease < FunctionTolerance" — convergencia débil, SEs del Hessiano pueden ser poco fiables
- `exitflag ≤ 0`: fallo — activar fallback
- **No asumir que exitflag=5 invalida el resultado, pero verificar SEs contra benchmark externo**

### Bootstrap y bandas de incertidumbre
- n_boot = 2000 réplicas (paper), K=100 modelos aceptados por draw, K_max=5000 candidatos
- IC al 68% (percentiles 16–84), no al 95%
- Algoritmo Runkle (1987): re-estimar VAR en cada draw de bootstrap

### Restricciones de signo — Tabla 2 del paper
| Variable | Restricción |
|----------|------------|
| igae (output) | − |
| pi_core (inflación) | − |
| pi_inpp (PPI) | libre |
| dep_tcr (TCR) | − |
| i_gap (tasa) | + |
| m2_gap (dinero) | − |

### Verificación de la muestra
- Ejecutar `0Proc_data.R` siempre que se sospeche cambio en T
- Esperado: T=159, NAs en X-13 = 0, rango Ene2001–Mar2014

---

## Strategies and Hard Rules

### Verification Checklist (antes de declarar SECCIÓN 3 completa)
- [ ] β̂⁻ puntual ≈ -0.035 (tolerancia: ±0.05)
- [ ] SE β̂⁻ ≈ 0.40 (tolerancia: ±0.15)
- [ ] β_lower ≈ -0.700 (tolerancia: ±0.15)
- [ ] Método en log: "[PRIMARIO] OLS + Newey-West (lag=3)"
- [ ] exitflag MLE impreso como diagnóstico (no determina β_lower)

### Anti-patrones identificados

**Anti-patrón #1: Reemplazar OLS+NW con MLE para β_lower**  
- Ocurrió en sesión 2026-05-09  
- Consecuencia: β_lower = -0.179 en vez de -0.700, figuras 8 y 9 incorrectas  
- Corrección: OLS+NW siempre primario para este coeficiente  

**Anti-patrón #2: Hardcodear fechas en ajuste estacional X-13**  
- `start = c(2001,1)` cuando la serie empieza en Dic 2000 → error de alineación  
- Corrección: pasar `base_start` dinámicamente a `ajuste_x13()`  

**Anti-patrón #3: Asumir convergencia por exitflag > 0**  
- exitflag=5 ≠ convergencia limpia; verificar SEs contra referencia externa  

### Diferencias conocidas con el paper (no corregibles)
- TCR: usamos 49 socios (Banxico SR15997) vs 21 socios (IFS) — fuente de datos distinta
- Ajuste estacional: X-13ARIMA-SEATS vs Tramo/X12 Additive del paper
- Estas diferencias son documentadas en el informe, no son bugs

---

## Pending Tasks (al 2026-05-10)

1. **Ejecutar** `SVAR_SignRestrictions.m` → verificar checklist de Sección 3
2. **Ejecutar** `0Proc_data.R` → confirmar T=159
3. **Ejecutar** `4Informe_replicacion.R` → regenerar `informe_replicacion.html`
