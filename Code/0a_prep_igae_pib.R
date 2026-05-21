# ==============================================================================
# 0a_prep_igae_pib.R
# Limpieza y preparación de las series base para el filtro de Kalman MFK
# Siguiendo: Elizondo (2012), Carrillo & Elizondo (2015/18)
#
# Inputs:
#   raw data/brecha_prod/igae_indice (2).xlsx   — IGAE total desest., base 2018=100
#   raw data/brecha_prod/pibt_cte_valor (2).xlsx — PIB trimestral desest., mmp 2018
#   raw data/brecha_prod/Brecha_producto.csv    — Brecha INEGI (para comparación)
#
# Output:
#   raw data/brecha_prod/kalman_input.xlsx       — Series limpias listas para MFK
# ==============================================================================

rm(list = ls())

paquetes <- c("readxl", "dplyr", "lubridate", "zoo", "openxlsx", "stringi")
sapply(paquetes, function(p) {
  if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
})

RAW <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data/brecha_prod/"

# ==============================================================================
# 1. IGAE total desestacionalizado (mensual, índice base 2018=100)
# ==============================================================================
# Estructura del xlsx (formato ancho, filas = indicadores, columnas = tiempo):
#   Fila col-names (automático): "INDICADOR GLOBAL DE LA ACTIVIDAD ECONÓMICA"
#   df[1,] = "Series desestacionalizadas"
#   df[2,] = "Índice Base 2018 = 100"
#   df[3,] = "Denominación" | NA | "1993" (año, sólo primera col de cada año)
#   df[4,] = NA | NA | "ENE" | "FEB" | ... (meses, repiten)
#   df[5,] = "Indicador Global de la Actividad Económica" | NA | 56.4 | 57.5 | ...

raw_igae <- read_excel(paste0(RAW, "igae_indice (2).xlsx"))

# Extraer vectores de año y mes (columnas 3 en adelante)
year_row  <- as.character(unlist(raw_igae[3, -(1:2)]))
month_row <- as.character(unlist(raw_igae[4, -(1:2)]))

# Rellenar años hacia adelante (fill forward)
year_filled <- na.locf(ifelse(is.na(year_row) | year_row == "NA", NA, year_row))

# Mapeo mes español → número
mes_map <- c(ENE = 1, FEB = 2, MAR = 3, ABR = 4, MAY = 5, JUN = 6,
             JUL = 7, AGO = 8, SEP = 9, OCT = 10, NOV = 11, DIC = 12)
month_num <- mes_map[month_row]

# Construir vector de fechas
dates_igae <- as.Date(paste(year_filled, month_num, "01", sep = "-"), "%Y-%m-%d")

# Extraer valores IGAE total (fila 5 del dataframe = "Indicador Global...")
igae_vals <- suppressWarnings(as.numeric(as.character(unlist(raw_igae[5, -(1:2)]))))

igae_clean <- data.frame(date = dates_igae, igae = igae_vals) %>%
  filter(!is.na(date), !is.na(igae)) %>%
  arrange(date)

cat("=== IGAE total desestacionalizado ===\n")
cat("Período:", format(min(igae_clean$date), "%b %Y"), "—",
    format(max(igae_clean$date), "%b %Y"), "\n")
cat("Observaciones:", nrow(igae_clean), "\n")
cat("Valores iniciales (Jan-Jun 1993):", round(head(igae_clean$igae, 6), 3), "\n")
cat("Rango:", round(range(igae_clean$igae), 3), "\n\n")

# ==============================================================================
# 2. PIB trimestral desestacionalizado (millones de pesos constantes 2018)
# ==============================================================================
# Estructura similar:
#   df[3,] = "Denominación" | NA | "1993" | NA | NA | NA | "1994" | ...
#   df[4,] = NA | NA | "I" | "II" | "III" | "IV" | "I" | ...
#   df[5,] = "Producto interno bruto" | NA | 1380000 | ...

raw_pib <- read_excel(paste0(RAW, "pibt_cte_valor (2).xlsx"))

year_row_p <- as.character(unlist(raw_pib[3, -(1:2)]))
qtr_row_p  <- as.character(unlist(raw_pib[4, -(1:2)]))

year_filled_p <- na.locf(ifelse(is.na(year_row_p) | year_row_p == "NA", NA, year_row_p))

# Trimestre → primer mes del trimestre
qtr_map <- c("I" = 1, "II" = 4, "III" = 7, "IV" = 10)
month_p <- qtr_map[qtr_row_p]

dates_pib <- as.Date(paste(year_filled_p, month_p, "01", sep = "-"), "%Y-%m-%d")

pib_vals <- suppressWarnings(as.numeric(as.character(unlist(raw_pib[5, -(1:2)]))))

pib_clean <- data.frame(date = dates_pib, pib = pib_vals) %>%
  filter(!is.na(date), !is.na(pib)) %>%
  arrange(date)

cat("=== PIB trimestral desestacionalizado ===\n")
cat("Período:", format(min(pib_clean$date), "%b %Y"), "—",
    format(max(pib_clean$date), "%b %Y"), "\n")
cat("Observaciones:", nrow(pib_clean), "trimestres\n")
cat("Valores iniciales (Q1-Q4 1993, mmp):",
    round(head(pib_clean$pib, 4) / 1e6, 4), "\n")
cat("Q1 1993 (ancla Kalman):", pib_clean$pib[1], "mmp de 2018\n\n")

# ==============================================================================
# 3. Brecha INEGI (para comparación — IGAE excluyendo sector petrolero)
# ==============================================================================
# Estructura del CSV (encoding latin1, skip primeras 4 filas):
# Col 1: Año | Col 2: Trimestre | Col 3: Mes
# Col 4: PIB gap total | Col 5: IGAE gap total
# Col 6: PIB gap ex-petrol | Col 7: IGAE gap ex-petrol ← serie en uso

raw_bp <- tryCatch(
  read.csv(paste0(RAW, "Brecha_producto.csv"),
           skip = 4, header = FALSE, fileEncoding = "latin1",
           stringsAsFactors = FALSE),
  error = function(e) {
    read.csv(paste0(RAW, "Brecha_producto.csv"),
             skip = 4, header = FALSE,
             stringsAsFactors = FALSE)
  }
)

# Limpiar nombres
names(raw_bp) <- c("anio", "trimestre", "mes",
                   "pib_gap", "igae_gap_total",
                   "pib_gap_nopetrol", "igae_gap_nopetrol")

# Eliminar filas de notas (al final)
raw_bp <- raw_bp %>%
  filter(!is.na(anio), anio != "",
         !grepl("a\\.", anio, ignore.case = TRUE),
         !grepl("Fuente", anio, ignore.case = TRUE)) %>%
  mutate(
    anio = as.numeric(trimws(anio)),
    mes_clean = tolower(stri_trans_general(trimws(mes), "Latin-ASCII")),
    mes_num = case_when(
      grepl("^ene", mes_clean) ~ 1,  grepl("^feb", mes_clean) ~ 2,
      grepl("^mar", mes_clean) ~ 3,  grepl("^abr", mes_clean) ~ 4,
      grepl("^may", mes_clean) ~ 5,  grepl("^jun", mes_clean) ~ 6,
      grepl("^jul", mes_clean) ~ 7,  grepl("^ago", mes_clean) ~ 8,
      grepl("^sep", mes_clean) ~ 9,  grepl("^oct", mes_clean) ~ 10,
      grepl("^nov", mes_clean) ~ 11, grepl("^dic", mes_clean) ~ 12,
      TRUE ~ NA_real_
    ),
    date = make_date(year = anio, month = mes_num, day = 1),
    igae_gap_total    = suppressWarnings(as.numeric(trimws(igae_gap_total))),
    igae_gap_nopetrol = suppressWarnings(as.numeric(trimws(igae_gap_nopetrol)))
  ) %>%
  filter(!is.na(date), !is.na(igae_gap_nopetrol)) %>%
  select(date, igae_gap_total, igae_gap_nopetrol) %>%
  arrange(date)

cat("=== Brecha INEGI (para comparación) ===\n")
cat("Período:", format(min(raw_bp$date), "%b %Y"), "—",
    format(max(raw_bp$date), "%b %Y"), "\n")
cat("Observaciones:", nrow(raw_bp), "\n")
cat("IGAE gap ex-petróleo (rango):", round(range(raw_bp$igae_gap_nopetrol, na.rm=TRUE), 2), "\n\n")

# ==============================================================================
# 4. Guardar series limpias en kalman_input.xlsx
# ==============================================================================
wb <- createWorkbook()

addWorksheet(wb, "igae_mensual")
writeData(wb, "igae_mensual", igae_clean)

addWorksheet(wb, "pib_trimestral")
writeData(wb, "pib_trimestral", pib_clean)

addWorksheet(wb, "inegi_brecha")
writeData(wb, "inegi_brecha", raw_bp)

saveWorkbook(wb, file = paste0(RAW, "kalman_input.xlsx"), overwrite = TRUE)
cat("Guardado: kalman_input.xlsx\n")
cat("  - Hoja 'igae_mensual':  ", nrow(igae_clean), "obs\n")
cat("  - Hoja 'pib_trimestral':", nrow(pib_clean), "obs\n")
cat("  - Hoja 'inegi_brecha':  ", nrow(raw_bp), "obs\n")
