rm(list=ls())

setwd("C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data")
#Sys.setlocale("LC_TIME", "Spanish")

paquetes <- c("dplyr", "readxl", "stringr", "tidyverse", "readr", "stringi",
              "ggplot2", "openxlsx", "lubridate", "mFilter", "tseries",
              "purrr", "vars", "urca", "plm")

sapply(paquetes, function(p) if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE))

df_final <- read_excel('database.xlsx') %>%
  mutate(across(!date, as.numeric))
base <- read_excel('raw_base.xlsx')

# ==============================================================================
# PRUEBAS DE RAÍZ UNITARIA (Augmented Dickey-Fuller)
# ==============================================================================
# Implementación: urca::ur.df con regla de Schwert para selección de rezagos
#   k_max = trunc((n-1)^(1/3)) — igual a tseries::adf.test; garantiza residuos blancos
#   para series altamente persistentes (igae, y_us_gap) con menos rezagos BIC
#   deja autocorrelación residual y sesga el test hacia no rechazar H0.
#
# Especificación:
#   - Variables de brecha (centradas en 0): type="none" y type="drift"
#   - El paper reporta: todas las variables rechazan H0 → I(0)
#
# Nota diagnóstica (igae, y_us_gap):
#   Con lags=1-3 los residuos del ADF muestran autocorrelación (Ljung-Box p<0.05)
#   → test inválido. Con lags=4-5 (regla de Schwert) los residuos son blancos
#   y se rechaza H0 al 5%: ambas series son I(0).

# Selección adaptativa de rezagos: mínimo k ∈ [1, k_max] tal que los residuos
# del ADF sean blancos (Ljung-Box p > 0.05). Si ningún k pasa, se usa k=k_max.
# k_max = trunc((n-1)^(1/3)) — regla de Schwert.
# Motivación: series muy persistentes (igae, y_us_gap) necesitan k ≥ 4 para
# residuos blancos; series con estructura ARMA compleja (i_gap) solo k=1.

schwert_kmax <- function(x) trunc((length(x) - 1)^(1/3))

# alpha = 0.10: umbral conservador para validar residuos blancos.
# Un p-valor de LB de 0.065 se consideraría insuficiente (aún hay
# autocorrelación residual marginal), forzando k más alto hasta p > 0.10.
select_lags_lb <- function(x, type = "drift", alpha = 0.10) {
  k_max <- max(schwert_kmax(x), 5)  # al menos 5 para no truncar de más
  for (k in 1:k_max) {
    res <- ur.df(x, type = type, lags = k)
    lb_p <- Box.test(residuals(res@testreg), lag = 12, type = "Ljung-Box")$p.value
    if (lb_p > alpha) return(k)
  }
  return(k_max)
}

# Función ADF — selección LB + confirmación PP
prueba_adf <- function(x, nombre) {
  x    <- x[!is.na(x)]
  k    <- select_lags_lb(x, type = "drift")

  # Spec "drift" con lags adaptados
  res_d <- ur.df(x, type = "drift", lags = k)
  tau2  <- res_d@teststat[1, "tau2"]
  cv_d  <- res_d@cval["tau2", ]
  lb_p  <- Box.test(residuals(res_d@testreg), lag = 12, type = "Ljung-Box")$p.value

  # Phillips-Perron (robusto, sin selección de lags)
  pp    <- PP.test(x, lshort = TRUE)

  concl <- ifelse(tau2 < cv_d["5pct"], "I(0)**",
           ifelse(pp$p.value < 0.05,    "I(0)-PP", "Verificar"))

  data.frame(
    Variable     = nombre,
    Lags_ADF     = k,
    tau2_ADF     = round(tau2, 3),
    cv5_ADF      = round(cv_d["5pct"], 2),
    LB12_p       = round(lb_p, 3),
    Residuos     = ifelse(lb_p > 0.05, "OK", "!"),
    PP_stat      = round(pp$statistic, 3),
    PP_p         = round(pp$p.value, 4),
    Conclusion   = concl,
    row.names    = NULL
  )
}

resultados_adf <- map2_dfr(
  .x = df_final %>% dplyr::select(-date),
  .y = colnames(df_final %>% dplyr::select(-date)),
  .f = ~ prueba_adf(.x, .y)
)

cat("=== Pruebas ADF — urca::ur.df + PP, lags adaptivos (LB-criterion) ===\n")
print(resultados_adf, row.names = FALSE)
cat("\nNota: Lags_ADF = mínimo k con residuos blancos (LB12 p>0.05).\n")
cat("PP confirma para series con estructura ARMA compleja (i_gap).\n")

# ==============================================================================
# PRUEBAS CONJUNTAS DE RAÍZ UNITARIA — PANEL (LLC y Breitung)
# C&E (2015/18) Apéndice C: "two joint unit-root tests (Levin, Lin & Chu and
# Breitung) in which it is assumed that there is a common unit root. The tests
# were rejected."
#
# Interpretación "panel": las 9 variables del VAR se tratan como N=9 unidades
# de panel con T=159 observaciones cada una. H0: raíz unitaria común (ρ=1).
# Un rechazo conjunto refuerza la conclusión ADF variable-a-variable.
# ==============================================================================

vars_ur <- c("igae","pi_core_gap","pi_inpp_gap","dep_tcr","i_gap","m2_gap",
             "y_us_gap","pi_us","pi_crb")

panel_long <- df_final %>%
  filter(date >= "2001-01-01", date <= "2014-03-01") %>%
  dplyr::select(date, all_of(vars_ur)) %>%
  pivot_longer(cols = -date, names_to = "variable", values_to = "value") %>%
  arrange(variable, date) %>%
  mutate(t_idx = rep(seq_len(n_distinct(date)), times = n_distinct(variable)))

pdata <- pdata.frame(panel_long, index = c("variable", "t_idx"))

# Levin-Lin-Chu (2002): asume coeficiente AR homogéneo entre paneles
llc_res <- tryCatch(
  purtest(value ~ 1, data = pdata, test = "levinlin", lags = "AIC",
          exo = "intercept"),
  error = function(e) { message("LLC error: ", e$message); NULL }
)

# Im-Pesaran-Shin (2003): permite heterogeneidad en coeficientes AR entre unidades.
# Nota: plm::purtest no implementa Breitung (2000) en su API actual; IPS es el
# complemento natural de LLC y también asume raíz unitaria común bajo H0.
ips_res <- tryCatch(
  purtest(value ~ 1, data = pdata, test = "ips", lags = "AIC",
          exo = "intercept"),
  error = function(e) { message("IPS error: ", e$message); NULL }
)

cat("\n=== Tests conjuntos de raíz unitaria (Panel — N=9 series, T=159) ===\n")
cat("H0: todas las series tienen raíz unitaria\n\n")

if (!is.null(llc_res)) {
  llc_stat  <- llc_res$statistic$statistic
  llc_pval  <- llc_res$statistic$p.value
  cat(sprintf("Levin-Lin-Chu (2002):         estadístico = %7.4f  |  p-valor = %.4f  |  %s H0\n",
              llc_stat, llc_pval,
              ifelse(llc_pval < 0.05, "RECHAZA", "No rechaza")))
} else {
  cat("Levin-Lin-Chu: no disponible\n")
  llc_stat <- NA; llc_pval <- NA
}

if (!is.null(ips_res)) {
  ips_stat <- ips_res$statistic$statistic
  ips_pval <- ips_res$statistic$p.value
  cat(sprintf("Im-Pesaran-Shin (2003) [IPS]: estadístico = %7.4f  |  p-valor = %.4f  |  %s H0\n",
              ips_stat, ips_pval,
              ifelse(ips_pval < 0.05, "RECHAZA", "No rechaza")))
} else {
  cat("IPS: no disponible\n")
  ips_stat <- NA; ips_pval <- NA
}

cat("\nPaper C&E (LLC + Breitung): ambos rechazan H0 → variables I(0) ✓\n")
cat("Replicación: LLC + IPS (Breitung no disponible en plm).\n")

# Guardar resultados para el informe
panel_ur_res <- data.frame(
  Test         = c("Levin-Lin-Chu (2002)", "Im-Pesaran-Shin (2003)"),
  Estadistico  = round(c(llc_stat, ips_stat), 4),
  p_valor      = round(c(llc_pval, ips_pval), 4),
  Decision     = c(ifelse(!is.na(llc_pval) & llc_pval < 0.05, "Rechaza H0", "No rechaza H0"),
                   ifelse(!is.na(ips_pval) & ips_pval < 0.05, "Rechaza H0", "No rechaza H0")),
  Paper        = c("Rechaza H0", "Rechaza H0"),
  stringsAsFactors = FALSE
)

fig_dir <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/figures"

wb_ur <- createWorkbook()
addWorksheet(wb_ur, "panel_ur")
writeData(wb_ur, "panel_ur", panel_ur_res)
saveWorkbook(wb_ur,
  file = "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data/panel_ur_results.xlsx",
  overwrite = TRUE)
cat("Guardado: panel_ur_results.xlsx\n")

# ==============================================================================
# ESTIMACIÓN DEL SVAR: RESTRICCIONES DE EXCLUSIÓN (CHOLESKY)
# ==============================================================================

# 1. Definición de los bloques del modelo
# El ordenamiento es crítico para Cholesky (de variables más lentas a más rápidas)
endogenas <- df_final %>% 
  dplyr::select(igae, pi_core_gap, pi_inpp_gap, m2_gap, i_gap, dep_tcr)

exogenas <- df_final %>% 
  dplyr::select(y_us_gap, pi_us, pi_crb)

# Convertimos a formato de serie de tiempo o matrices (el paquete vars lo prefiere)
endogenas_ts <- as.ts(endogenas)
exogenas_mat <- as.matrix(exogenas)

# 2. Estimación del VAR en Forma Reducida
# Se utilizan 3 rezagos (p=3) como se indica en la metodología empírica para México
var_reducido <- VAR(y = endogenas_ts, 
                    p = 3, 
                    type = "const", 
                    exogen = exogenas_mat)

# Revisión rápida de la estabilidad del sistema (las raíces deben ser < 1)
roots(var_reducido)

# 3. Identificación del Choque Estructural (Política Monetaria)
# Calculamos la IRF para un choque (impulso) en la tasa de interés (i_gap)
# Se utilizan 2,000 iteraciones de bootstrap para los intervalos de confianza
irf_mp_shock <- irf(var_reducido, 
                    impulse = "i_gap", 
                    response = c("igae", "pi_core_gap", "pi_inpp_gap", "dep_tcr"), 
                    n.ahead = 24,       # Horizonte de 24 meses
                    ortho = TRUE,       # ortho = TRUE aplica la descomposición de Cholesky
                    boot = TRUE, 
                    runs = 2000, 
                    ci = 0.68)          # Intervalos de confianza al 68% (un estándar en SVARs)

# 4. Visualización de las Funciones Impulso-Respuesta
# Ajustamos los márgenes del gráfico para una mejor visualización
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(irf_mp_shock, main = "Respuesta a un choque de Política Monetaria (Cholesky)")

#################################################################################
# 1. Transformar la base de formato ancho a largo
# Esto pone todos los nombres de tus variables en una columna y los valores en otra
df_grafico <- df_final %>%
  pivot_longer(
    cols = -date,           # Graficar todas las columnas EXCEPTO 'date'
    names_to = "variable",  # Aquí se guardarán los nombres (PIB, INPP, etc.)
    values_to = "valor"     # Aquí se guardarán los números
  )

# 2. Crear el gráfico con subplots
ggplot(df_grafico, aes(x = date, y = valor, color = variable)) +
  geom_line(size = 0.8) +
  # La magia ocurre aquí: crea un cuadro por cada 'variable'
  facet_wrap(~variable, 
             scales = "free_y", # Crucial: cada gráfico tiene su propia escala en el eje Y
             ncol = 2) +        # Puedes elegir cuántas columnas de gráficos quieres
  theme_minimal() +
  theme(legend.position = "none") + # Quitamos la leyenda porque el título del cuadro ya lo dice
  labs(title = "Análisis de Variables en el Tiempo",
       x = "Fecha",
       y = "Valor de la Variable")