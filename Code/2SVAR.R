rm(list=ls())

setwd("C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data")
#Sys.setlocale("LC_TIME", "Spanish")

paquetes <- c("dplyr", "readxl", "stringr", "tidyverse", "readr",'stringi',
              "ggplot2",  "openxlsx", "lubridate",'mFilter','tseries','purrr','vars')

sapply(paquetes, function(p) if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE))

df_final <- read_excel('database.xlsx')


endogenas_mat <- df_final %>% 
  dplyr::select(igae, pi_core_gap, pi_inpp_gap, m2_gap, i_gap, dep_tcr) %>% 
  mutate_all(as.numeric) %>% 
  as.matrix()

exogenas_mat <- df_final %>% 
  dplyr::select(y_us_gap, pi_us, pi_crb) %>% 
  mutate_all(as.numeric) %>% 
  as.matrix()

# ==============================================================================
# 1. ESTIMACIÓN DE LOS DOS MODELOS SVAR (Exclusión / Cholesky)
# ==============================================================================

# Matrices (Asegúrate de haber corrido tu código previo para crearlas)
# endogenas_mat y exogenas_mat ya deben estar en el ambiente

# Modelo 1: VAR Simple (Sin variables exógenas)
var_simple <- VAR(y = endogenas_mat, 
                  p = 3, 
                  type = "const")

# Modelo 2: VAR con variables exógenas (El propuesto en el paper)
var_exo <- VAR(y = endogenas_mat, 
               p = 3, 
               type = "const", 
               exogen = exogenas_mat)

# ==============================================================================
# 2. CÁLCULO DE LAS FUNCIONES IMPULSO-RESPUESTA (Horizonte = 16 meses)
# ==============================================================================

# Vector de variables a extraer
vars_to_plot <- c("igae", "pi_core_gap", "pi_inpp_gap", "dep_tcr", "i_gap", "m2_gap")

# IRF para el VAR Simple
irf_simple <- irf(var_simple, 
                  impulse = "i_gap", 
                  response = vars_to_plot, 
                  n.ahead = 16, 
                  ortho = TRUE, 
                  boot = TRUE, runs = 1000, ci = 0.68)

# IRF para el VAR con Exógenas
irf_exo <- irf(var_exo, 
               impulse = "i_gap", 
               response = vars_to_plot, 
               n.ahead = 16, 
               ortho = TRUE, 
               boot = TRUE, runs = 1000, ci = 0.68)

# ==============================================================================
# 3. REPLICACIÓN GRÁFICA DE LA FIGURA 7 (Solo Intervalos de Confianza)
# ==============================================================================

# Títulos de las gráficas
titulos <- c("Output", "Inflation", "PPI", "Real ER dep", "Nominal interest rate", "Money growth")

# Configuración de la cuadrícula 2x3
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

for(i in 1:6){
  v <- vars_to_plot[i]
  
  # 1. Extraer datos del VAR con Exógenas (Solo intervalos)
  low_exo  <- irf_exo$Lower$i_gap[, v]
  upp_exo  <- irf_exo$Upper$i_gap[, v]
  
  # 2. Extraer datos del VAR Simple (Solo intervalos)
  low_simple <- irf_simple$Lower$i_gap[, v]
  upp_simple <- irf_simple$Upper$i_gap[, v]
  
  # 3. Calcular límites dinámicos para el eje Y abarcando AMBOS intervalos
  y_min <- min(c(low_exo, low_simple))
  y_max <- max(c(upp_exo, upp_simple))
  rango <- y_max - y_min
  y_lims <- c(y_min - 0.1 * rango, y_max + 0.1 * rango)
  
  y_label <- ifelse(i %in% c(1, 4), "Percent points from trend", "")
  
  # 4. Crear el lienzo vacío
  plot(0:16, low_exo, type = "n", ylim = y_lims,
       main = titulos[i],
       xlab = "Months after the shock", 
       ylab = y_label,
       frame.plot = TRUE, cex.main = 1.3)
  
  # 5. Polígono del VAR Simple (Gris translúcido)
  polygon(c(0:16, rev(0:16)), c(low_simple, rev(upp_simple)), 
          col = rgb(0.5, 0.5, 0.5, 0.3), border = "darkgray", lty = 2)
  
  # 6. Polígono del VAR con Exógenas (Azul claro translúcido)
  polygon(c(0:16, rev(0:16)), c(low_exo, rev(upp_exo)), 
          col = rgb(0.2, 0.5, 0.9, 0.4), border = NA)
  
  # 7. Línea de referencia en cero
  abline(h = 0, col = "black", lty = 3, lwd = 1.5)
  
  # 8. Agregar leyenda visual
  if(i == 6){
    legend("bottomright", legend = c("CI: VAR w exo", "CI: Simple VAR"),
           fill = c(rgb(0.2, 0.5, 0.9, 0.4), rgb(0.5, 0.5, 0.5, 0.3)),
           border = c(NA, "darkgray"), bty = "n", cex = 1.1)
  }
}

# Título Principal
mtext("SVAR: Comparación de Intervalos de Confianza (Cholesky)", 
      side = 3, outer = TRUE, cex = 1.5, font = 2)