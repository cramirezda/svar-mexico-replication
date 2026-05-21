# ==============================================================================
# 0b_brecha_kalman.R
# Estimación de la brecha del producto mensual via Filtro de Kalman (MFK)
# Siguiendo: Elizondo (2012), Carrillo & Elizondo (2015/18)
#
# Modelo espacio-estado (univariado):
#   Ec. de transición: X_t = A*X_{t-1} + w_t,  w_t ~ N(0, Q)
#   Ec. de medición:   Z_t = H*X_t    + v_t,  v_t ~ N(0, R)
#
#   X_t = dif[log(PIB_t)]  (tasa de crecimiento mensual del PIB, no observable)
#   Z_t = dif[log(IGAE_t)] (tasa de crecimiento del IGAE, observable)
#
# Parámetros {A, H, Q, R} estimados conjuntamente por máxima verosimilitud.
# Método MFK (Método de Filtro de Kalman, Elizondo 2012 Sección 4.3.2).
#
# Pipeline:
#   1. Leer series limpias de kalman_input.xlsx
#   2. Calcular tasas de crecimiento del IGAE (Z_t)
#   3. Estimar {A, H, Q, R} via dlmMLE
#   4. Suavizado de Kalman → X_hat_t (crec. mensual PIB estimado)
#   5. Reconstruir log PIB mensual (ancla: Q1-1993 = PIB trimestral observado)
#   6. Filtro HP (lambda=129,600) → tendencia
#   7. Brecha = (log_PIB - tendencia) * 100
#   8. Comparación vs brecha INEGI y serie actual
#   9. Guardar igae_kalman_gap.xlsx
#
# Input:  raw data/brecha_prod/kalman_input.xlsx
# Output: raw data/brecha_prod/igae_kalman_gap.xlsx
# ==============================================================================

rm(list = ls())

paquetes <- c("readxl", "dplyr", "lubridate", "dlm", "mFilter",
              "openxlsx", "ggplot2", "tidyr")
sapply(paquetes, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
})

RAW     <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data/brecha_prod/"
FIGDIR  <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/figures/"

# ==============================================================================
# 1. Leer series limpias
# ==============================================================================
igae_df  <- read_excel(paste0(RAW, "kalman_input.xlsx"), sheet = "igae_mensual")
pib_df   <- read_excel(paste0(RAW, "kalman_input.xlsx"), sheet = "pib_trimestral")
inegi_df <- read_excel(paste0(RAW, "kalman_input.xlsx"), sheet = "inegi_brecha")

igae_df  <- igae_df  %>% arrange(date)
pib_df   <- pib_df   %>% arrange(date)
inegi_df <- inegi_df %>% arrange(date)

cat("IGAE mensual: ", format(min(igae_df$date), "%b-%Y"), "—",
    format(max(igae_df$date), "%b-%Y"), "\n")
cat("PIB trimestral:", format(min(pib_df$date), "%b-%Y"), "—",
    format(max(pib_df$date), "%b-%Y"), "\n\n")

# ==============================================================================
# 2. Tasa de crecimiento del IGAE — observable Z_t
# ==============================================================================
# Z_t = log(IGAE_t) - log(IGAE_{t-1})
Z <- diff(log(igae_df$igae))

# Fechas para Z (empiezan en el segundo mes de los datos IGAE)
dates_Z <- igae_df$date[-1]

cat("Z (dif log IGAE): n =", length(Z), "\n")
cat("Período Z:", format(min(dates_Z), "%b-%Y"), "—",
    format(max(dates_Z), "%b-%Y"), "\n")
cat("Media:", round(mean(Z), 6), " SD:", round(sd(Z), 6), "\n\n")

# ==============================================================================
# 3. Modelo MFK — Estimación por máxima verosimilitud (paquete dlm)
# ==============================================================================
# Parametrización:
#   par[1] = A      (coeficiente AR, libre; limitamos con bounds)
#   par[2] = H      (loading de medición, libre)
#   par[3] = log(Q) (varianza ruido de estado, log para positividad)
#   par[4] = log(R) (varianza ruido de medición, log para positividad)

build_mfk <- function(par) {
  A <- par[1]
  H <- par[2]
  Q <- exp(par[3])
  R <- exp(par[4])

  dlm(
    FF = matrix(H),        # carga observable
    GG = matrix(A),        # transición de estado
    V  = matrix(R),        # varianza ruido de medición
    W  = matrix(Q),        # varianza ruido de estado
    m0 = matrix(0),        # media inicial del estado
    C0 = matrix(1e7)       # covarianza inicial difusa
  )
}

# Valores iniciales razonables
var_Z  <- var(Z)
init_par <- c(
  A    = 0.7,
  H    = 1.0,
  logQ = log(var_Z * 0.3),
  logR = log(var_Z * 0.7)
)

# Nota sobre identificación: H (loading de medición) debe restringirse
# económicamente. Como Z_t = dif[log(IGAE)] y X_t = dif[log(PIB)], ambas son
# tasas de crecimiento mensual comparables → H ∈ [0.5, 2.0].
# Elizondo (2012) estima los cuatro parámetros pero la solución con H libre
# tiende a degenerarse; acotamos H para garantizar identificación.

cat("Estimando modelo MFK via MLE...\n")
fit_mfk <- dlmMLE(
  y     = Z,
  parm  = init_par,
  build = build_mfk,
  lower = c(-0.99, 0.50, -Inf, -Inf),
  upper = c( 0.99, 2.00,  Inf,  Inf),
  control = list(maxit = 2000, factr = 1e7)
)

if (fit_mfk$convergence != 0) {
  warning("MLE no convergió (code ", fit_mfk$convergence, "). Revisar valores iniciales.")
} else {
  cat("Convergencia: OK\n")
}

A_hat <- fit_mfk$par[1]
H_hat <- fit_mfk$par[2]
Q_hat <- exp(fit_mfk$par[3])
R_hat <- exp(fit_mfk$par[4])

# Verificar solución no degenerada
if (Q_hat / R_hat > 1000 | R_hat < 1e-10) {
  message("ADVERTENCIA: solución posiblemente degenerada (Q/R=", round(Q_hat/R_hat,1), ").")
  message("Alternativa: fijando H = 1 (modelo de 3 parámetros).")
  build_mfk_h1 <- function(par) {
    dlm(FF = matrix(1), GG = matrix(par[1]),
        V  = matrix(exp(par[3])), W = matrix(exp(par[2])),
        m0 = matrix(0), C0 = matrix(1e7))
  }
  fit_mfk <- dlmMLE(Z, parm = c(0.7, log(var_Z*0.3), log(var_Z*0.7)),
                    build = build_mfk_h1,
                    lower = c(-0.99, -Inf, -Inf),
                    upper = c( 0.99,  Inf,  Inf),
                    control = list(maxit = 2000))
  build_mfk <- build_mfk_h1
  A_hat <- fit_mfk$par[1]; H_hat <- 1.0
  Q_hat <- exp(fit_mfk$par[2]); R_hat <- exp(fit_mfk$par[3])
  cat("Reestimación con H=1: convergencia =", fit_mfk$convergence, "\n")
}

cat("\n--- Parámetros estimados MFK ---\n")
cat(sprintf("  A  = %8.5f  (persistencia AR del crec. mensual PIB)\n", A_hat))
cat(sprintf("  H  = %8.5f  (loading IGAE sobre PIB)\n",                H_hat))
cat(sprintf("  Q  = %8.6f  (varianza ruido de estado)\n",              Q_hat))
cat(sprintf("  R  = %8.6f  (varianza ruido de medición)\n",            R_hat))
cat(sprintf("  Ratio Q/R = %.4f\n", Q_hat / R_hat))

# ==============================================================================
# 4. Suavizado de Kalman → X_hat (crec. mensual PIB estimado)
# ==============================================================================
mod_fit <- build_mfk(fit_mfk$par)
kf  <- dlmFilter(Z, mod_fit)
ks  <- dlmSmooth(kf)

# ks$s tiene n+1 filas (incluye estado inicial); quitamos la primera
# as.numeric() aplana a vector antes de indexar (compatible con estado 1D)
X_hat <- as.numeric(ks$s)[-1]

cat("\nX_hat (crec. mensual PIB suavizado): n =", length(X_hat), "\n")
cat("Media:", round(mean(X_hat), 6), " SD:", round(sd(X_hat), 6), "\n")

# ==============================================================================
# 5. Reconstrucción del log PIB mensual
# ==============================================================================
# Ancla: log(PIB_Q1_1993) como log PIB en Marzo 1993
# (Elizondo 2012, p.11: "se considera el dato correspondiente al logaritmo
#  del PIB trimestral de marzo de 1993 como valor inicial")

pib_q1_1993 <- pib_df$pib[pib_df$date == as.Date("1993-01-01")]
if (length(pib_q1_1993) == 0) {
  stop("No se encontró el valor del PIB en Q1 1993 (fecha 1993-01-01). Verificar datos.")
}
log_pib_anchor <- log(pib_q1_1993)  # log millones de pesos 2018

# X_hat[t] = dif[log_PIB] para la fecha dates_Z[t]
# X_hat[1] → Feb 1993, X_hat[2] → Mar 1993, etc.
anchor_idx <- which(dates_Z == as.Date("1993-03-01"))

if (length(anchor_idx) == 0) {
  # Si IGAE empieza después de Jan 1993, buscar el trimestre correspondiente
  anchor_idx <- 1
  warning("Ancla en Marzo 1993 no encontrada; usando primera observación disponible.")
}

n_z <- length(X_hat)
log_pib_monthly <- numeric(n_z)

# Anclar en Marzo 1993
log_pib_monthly[anchor_idx] <- log_pib_anchor

# Reconstruir hacia adelante
if (anchor_idx < n_z) {
  for (t in (anchor_idx + 1):n_z) {
    log_pib_monthly[t] <- log_pib_monthly[t - 1] + X_hat[t]
  }
}

# Reconstruir hacia atrás (Ene-Feb 1993)
if (anchor_idx > 1) {
  for (t in (anchor_idx - 1):1) {
    log_pib_monthly[t] <- log_pib_monthly[t + 1] - X_hat[t + 1]
  }
}

pib_monthly_df <- data.frame(
  date      = as.Date(dates_Z),   # POSIXct → Date para evitar comparaciones por TZ
  log_pib   = log_pib_monthly,
  pib_nivel = exp(log_pib_monthly)
)

cat("\nLog PIB mensual reconstruido:\n")
cat("  Período:", format(min(pib_monthly_df$date), "%b-%Y"), "—",
    format(max(pib_monthly_df$date), "%b-%Y"), "\n")
cat("  Ancla (Mar 1993):", round(log_pib_anchor, 4),
    "= log(", round(pib_q1_1993 / 1e6, 3), "billones)\n")

# ==============================================================================
# 6. Filtro Hodrick-Prescott (lambda = 129,600 — estándar mensual)
# ==============================================================================
hp_result <- hpfilter(pib_monthly_df$log_pib, freq = 129600, type = "lambda")

pib_monthly_df$hp_trend <- as.numeric(hp_result$trend)
pib_monthly_df$hp_cycle <- as.numeric(hp_result$cycle)

# Brecha en puntos porcentuales
pib_monthly_df$igae_kalman_gap <- pib_monthly_df$hp_cycle * 100

cat("\nBrecha Kalman (log_pib - HP_trend) * 100:\n")
cat("  Min:", round(min(pib_monthly_df$igae_kalman_gap), 3), "\n")
cat("  Max:", round(max(pib_monthly_df$igae_kalman_gap), 3), "\n")
cat("  SD: ", round(sd(pib_monthly_df$igae_kalman_gap), 3), "\n")

# ==============================================================================
# 7. Comparación: Kalman vs INEGI vs serie original Banxico
# ==============================================================================
# Serie original (Banxico pre-calculada)
banxico_df <- read_excel(
  "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data/Brecha_producto.xlsx"
) %>%
  rename(anio = Años, mes = Meses, igae_banxico = IGAE) %>%
  mutate(
    mes_clean = tolower(stringi::stri_trans_general(trimws(mes), "Latin-ASCII")),
    mes_num = case_when(
      grepl("^ene", mes_clean) ~ 1,  grepl("^feb", mes_clean) ~ 2,
      grepl("^mar", mes_clean) ~ 3,  grepl("^abr", mes_clean) ~ 4,
      grepl("^may", mes_clean) ~ 5,  grepl("^jun", mes_clean) ~ 6,
      grepl("^jul", mes_clean) ~ 7,  grepl("^ago", mes_clean) ~ 8,
      grepl("^sep", mes_clean) ~ 9,  grepl("^oct", mes_clean) ~ 10,
      grepl("^nov", mes_clean) ~ 11, grepl("^dic", mes_clean) ~ 12
    ),
    date = make_date(year = as.numeric(anio), month = mes_num, day = 1)
  ) %>%
  select(date, igae_banxico) %>%
  filter(!is.na(date))

# Unir todo
compare_df <- pib_monthly_df %>%
  select(date, igae_kalman_gap) %>%
  left_join(inegi_df %>% select(date, igae_gap_nopetrol), by = "date") %>%
  left_join(banxico_df, by = "date") %>%
  filter(date >= "2001-01-01", date <= "2014-03-01")

cat("\n--- Estadísticos descriptivos comparados (2001-2014) ---\n")
cat(sprintf("%-25s %8s %8s %8s\n", "Serie", "Media", "SD", "Rango"))
cat(strrep("-", 55), "\n")
for (col in c("igae_kalman_gap", "igae_gap_nopetrol", "igae_banxico")) {
  vals <- compare_df[[col]]
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    rng <- paste0("[", round(min(vals),2), ", ", round(max(vals),2), "]")
    cat(sprintf("%-25s %8.4f %8.4f %15s\n", col,
                mean(vals), sd(vals), rng))
  }
}

# ==============================================================================
# 8. Gráfico comparación
# ==============================================================================
compare_long <- compare_df %>%
  pivot_longer(cols = c(igae_kalman_gap, igae_gap_nopetrol, igae_banxico),
               names_to = "serie", values_to = "valor") %>%
  mutate(serie = recode(serie,
    igae_kalman_gap   = "Kalman MFK (C&E)",
    igae_gap_nopetrol = "INEGI ex-petróleo",
    igae_banxico      = "Banxico (actual)"
  ))

p <- ggplot(compare_long, aes(x = date, y = valor, color = serie, linetype = serie)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_color_manual(values = c("Kalman MFK (C&E)" = "#1f78b4",
                                "INEGI ex-petróleo" = "#e31a1c",
                                "Banxico (actual)"  = "#33a02c")) +
  scale_linetype_manual(values = c("Kalman MFK (C&E)" = "solid",
                                   "INEGI ex-petróleo" = "dashed",
                                   "Banxico (actual)"  = "dotted")) +
  labs(title    = "Brecha del producto México (2001–2014)",
       subtitle = "Comparación: Kalman MFK vs INEGI vs Banxico",
       x = NULL, y = "% desviación del potencial",
       color = NULL, linetype = NULL,
       caption = "Fuentes: INEGI, Banxico. Metodología Kalman: Elizondo (2012) / C&E (2015).") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(paste0(FIGDIR, "brecha_kalman_comparacion.pdf"), p,
       width = 9, height = 5)
ggsave(paste0(FIGDIR, "brecha_kalman_comparacion.png"), p,
       width = 9, height = 5, dpi = 150)
cat("\nGráfico guardado: brecha_kalman_comparacion.pdf / .png\n")

# ==============================================================================
# 9. Guardar serie final (muestra completa + período VAR)
# ==============================================================================
# Serie completa (toda la muestra disponible para usos futuros)
output_full <- pib_monthly_df %>%
  select(date, log_pib, hp_trend, igae_kalman_gap) %>%
  rename(log_pib_kalman = log_pib,
         log_pib_hp_trend = hp_trend,
         igae_kalman = igae_kalman_gap)

# Serie para el VAR (Jan 2001 – Mar 2014, con un período extra de inicio)
output_var <- output_full %>%
  filter(date >= "2000-12-01", date <= "2014-03-01") %>%
  select(date, igae_kalman)

wb <- createWorkbook()
addWorksheet(wb, "serie_completa")
writeData(wb, "serie_completa", output_full)

addWorksheet(wb, "serie_var")
writeData(wb, "serie_var", output_var)

addWorksheet(wb, "parametros_mfk")
writeData(wb, "parametros_mfk", data.frame(
  parametro   = c("A", "H", "Q", "R", "convergencia"),
  valor       = c(A_hat, H_hat, Q_hat, R_hat, fit_mfk$convergence),
  descripcion = c("Coef. AR (transición de estado)",
                  "Loading medición (IGAE sobre PIB)",
                  "Varianza ruido de estado",
                  "Varianza ruido de medición",
                  "0 = convergencia exitosa")
))

saveWorkbook(wb, file = paste0(RAW, "igae_kalman_gap.xlsx"), overwrite = TRUE)

cat("\nGuardado: igae_kalman_gap.xlsx\n")
cat("  - Hoja 'serie_completa': ", nrow(output_full), "obs\n")
cat("  - Hoja 'serie_var':      ", nrow(output_var), "obs (muestra VAR)\n")
cat("  - Hoja 'parametros_mfk': parámetros estimados\n")
cat("\nPróximo paso: actualizar 0Proc_data.R para leer igae_kalman_gap.xlsx\n")
