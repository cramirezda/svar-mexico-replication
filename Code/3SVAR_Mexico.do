/*==============================================================================
  3SVAR_Mexico.do  — PARTE 1: DIAGNÓSTICOS + CHOLESKY (Figura 7)
  Replicación: Carrillo & Elizondo (2015), Sección 4 — México

  Estrategia híbrida:
    • Este script: Diagnósticos + Identificación por exclusión (Fig. 7)
    • MATLAB script (4SVAR_SignRestrictions.m): Restricciones de signos (Figs. 8-9)

  Orden Cholesky (lento → rápido en ajuste a choques nominales):
      1. igae        : Brecha del producto
      2. pi_core_gap : Brecha inflación subyacente (CPI core)
      3. pi_inpp_gap : Brecha inflación productor (PPI)
      4. dep_tcr     : Crecimiento del TCR (+ = depreciación)
      5. i_gap       : Brecha tasa de interés nominal  ← CHOQUE
      6. m2_gap      : Brecha crecimiento del dinero

  Variables exógenas (block-exógenas — Cushman & Zha 1997):
      y_us_gap  pi_us  pi_crb

  Muestra: Enero 2001 – Marzo 2014 | Frecuencia: Mensual | Rezagos: 3

  NOTA TCR: SR15997 (Banxico) = 49 socios comerciales. El paper usa 21.
  Convención paper: CAÍDA del índice = APRECIACIÓN del peso.
  dep_tcr > 0 → depreciación (coherente con definición estándar).
==============================================================================*/

clear all
set more off
set scheme s1mono        // gráficas tipo paper

* ── Directorios ──────────────────────────────────────────────────────────────
global wdir `"C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras"'
global ddir "$wdir/raw data"
global fdir "$wdir/figures"
global cdir "$wdir/Code"

capture mkdir "$fdir"
cd "$ddir"

/*==============================================================================
  BLOQUE 0: CARGA Y PREPARACIÓN DE DATOS
==============================================================================*/
di as txt _newline "=================================="
di as txt "BLOQUE 0: CARGA DE DATOS"
di as txt "=================================="

* Importar base procesada en R
import excel "database.xlsx", sheet("base") firstrow clear

* ── Conversión de fecha ───────────────────────────────────────────────────
* Intentar primero como string "YYYY-MM-DD"
capture confirm string variable date
if _rc == 0 {
    gen double date2 = date(date, "YMD")
    drop date
    rename date2 date
}
* Si viene como número serial de Excel (días desde 1899-12-30)
else {
    capture assert date > 36000 & date < 50000
    if _rc == 0 gen date = date + td(30dec1899)
}
format date %td

* Crear variable de tiempo mensual para tsset
gen int t = mofd(date)
format t %tm
tsset t, monthly
label var t "Fecha"

* Forzar numérico en todas las variables
foreach v in igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap ///
             y_us_gap pi_us pi_crb {
    capture destring `v', replace force
    label var `v' "`v'"
}

* Restringir a la muestra del paper
keep if tin(2001m1, 2014m3)

di as txt "Observaciones en muestra: " _N
summarize igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap ///
          y_us_gap pi_us pi_crb

* ── Exportar datos para MATLAB (formato CSV) ──────────────────────────────
* MATLAB leerá este CSV directamente (no necesita instalar nada)
export delimited using "$ddir/database_for_matlab.csv", ///
    delimiter(",") replace
di as txt "→ Datos exportados: $ddir/database_for_matlab.csv"

/*==============================================================================
  BLOQUE 1: PRUEBAS DE DIAGNÓSTICO
  Incluye: raíz unitaria, selección de rezagos, estabilidad, autocorrelación
==============================================================================*/
di as txt _newline "=================================="
di as txt "BLOQUE 1: DIAGNÓSTICOS"
di as txt "=================================="

* ── 1.1 Augmented Dickey-Fuller (HO: raíz unitaria) ─────────────────────
di as txt _newline "--- 1.1 ADF Test (HO: raíz unitaria) ---"
foreach v in igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap ///
             y_us_gap pi_us pi_crb {
    di as txt _newline "  Variable: `v'"
    quietly dfuller `v', lags(3) regress
    di as res "  ADF stat = " %8.4f `r(Zt)' "  p-valor aprox.: " ///
              cond(`r(Zt)' < -3.45, "< 0.01", cond(`r(Zt)' < -2.89, "< 0.05", "> 0.05"))
}

* ── 1.2 Phillips-Perron ─────────────────────────────────────────────────
di as txt _newline "--- 1.2 Phillips-Perron Test (HO: raíz unitaria) ---"
foreach v in igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap ///
             y_us_gap pi_us pi_crb {
    di as txt _newline "  Variable: `v'"
    pperron `v', lags(4)
}

* ── 1.3 Selección de rezagos ────────────────────────────────────────────
di as txt _newline "--- 1.3 Criterios de selección de rezagos (maxlag=6) ---"
di as txt "    (Paper: SIC→1, AIC→2, LR test desempata → 3 rezagos)"
varsoc igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap, ///
    exog(y_us_gap pi_us pi_crb) maxlag(6)

* ── 1.4 Estimación del VAR(3) para diagnósticos ────────────────────────
di as txt _newline "--- 1.4 Estimación VAR(3) con exógenas ---"
var igae pi_core_gap pi_inpp_gap dep_tcr i_gap m2_gap, ///
    lags(1/3) exog(y_us_gap pi_us pi_crb)
estimates store var_exo_diag

* Estabilidad: raíces del polinomio característico (deben ser < 1)
di as txt _newline "--- Estabilidad del VAR: raíces (deben ser < 1) ---"
varstable, graph name(varstable_plot, replace)
graph export "$fdir/varstable.png", name(varstable_plot) replace

* Test LM de autocorrelación de residuos (hasta lag 12)
di as txt _newline "--- Test LM Autocorrelación residuos (hasta lag 12) ---"
di as txt "    (HO: sin autocorrelación; no rechazar = residuos ruido blanco)"
varlmar, mlag(12)

* Test de normalidad de residuos (Jarque-Bera)
di as txt _newline "--- Test de Normalidad residuos (Jarque-Bera) ---"
varnorm

* Prueba de Portmanteau (Q de Ljung-Box)
di as txt _newline "--- Test Portmanteau (Ljung-Box Q) ---"
varwle

/*==============================================================================
  BLOQUE 2: RESTRICCIONES DE EXCLUSIÓN — MIGRADO A MATLAB

  La identificación por restricciones de exclusión (Cholesky) y la generación
  de la Figura 7 fueron migradas al script de MATLAB para unificar el flujo
  de estimación de IRFs en un solo entorno junto con las restricciones de signo.

  → Ver: Code/SVAR_SignRestrictions.m — SECCIÓN 4B y FIGURA 7

  El script MATLAB requiere database_for_matlab.csv (exportado en BLOQUE 0)
  y produce fig7_exclusion.pdf / fig7_exclusion.png en la carpeta figures/.
==============================================================================*/
di as txt _newline "=================================="
di as txt "BLOQUE 2: Migrado a MATLAB"
di as txt "(Ver SVAR_SignRestrictions.m — SECCIÓN 4B)"
di as txt "=================================="

/*==============================================================================
  BLOQUE 3: DIAGNÓSTICO DEL TIPO DE CAMBIO REAL (SR15997)
==============================================================================*/
di as txt _newline "=================================="
di as txt "BLOQUE 3: DIAGNÓSTICO TCR"
di as txt "=================================="

di as txt _newline "Serie SR15997 (Banxico): TCR multilateral, 49 socios, base 1990"
di as txt "Convención: UN AUMENTO en SR15997 = DEPRECIACIÓN del peso"
di as txt "dep_tcr = Δlog(SR15997)*1200 → valor positivo = DEPRECIACIÓN ✓"
di as txt ""
di as txt "Hipótesis sobre IRF invertida:"
di as txt "  (a) Diferencia de cobertura: 49 vs 21 socios → posible ruido"
di as txt "  (b) Diferencia de convención: verificar con correlación abajo"
di as txt ""

* Correlación diagnóstica TCR ~ tasa de interés
di as txt "Correlación dep_tcr vs i_gap (debe ser negativa post-choque):"
correlate dep_tcr i_gap pi_core_gap

* Estadísticas descriptivas
sum dep_tcr, detail

di as txt _newline "Si la IRF del TCR sale con signo opuesto al paper:"
di as txt "→ Agregar en Bloque 0: replace dep_tcr = -dep_tcr"
di as txt "  (esto invierte la convención del índice)"
di as txt ""
di as txt "Para Figuras 8 y 9 (restricciones de signos):"
di as txt "→ Correr el script MATLAB: Code/4SVAR_SignRestrictions.m"
di as txt "  El script de MATLAB lee database_for_matlab.csv"
di as txt "  y produce fig8 y fig9 en la carpeta figures/"

di as txt _newline "=================================="
di as txt "Script completado. Figuras en: $fdir"
di as txt "=================================="
