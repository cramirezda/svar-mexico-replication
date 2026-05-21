rm(list=ls())

setwd("C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data")
#Sys.setlocale("LC_TIME", "Spanish")

paquetes <- c("dplyr", "readxl", "stringr", "tidyverse", "readr",'stringi',
              "ggplot2",  "openxlsx", "lubridate",'mFilter','seasonal')

sapply(paquetes, function(p) if (!require(p, character.only = TRUE)) install.packages(p, dependencies = TRUE))

tcr <- read_excel('TCR.xlsx')%>%
  rename(date = Fecha,
         tcr = SR15997)%>%
  mutate(date = as_date(date, origin = "1899-12-30"))%>%
  filter(date <= '2014-03-01')


pcepilfe <- read_excel('PCEPI.xlsx')%>%
  rename(date= observation_date,
         pce= PCEPI)%>%
  mutate(date= as_date(date))%>%
  filter(date <= '2014-03-01')

m2 <- read_excel('M2.xlsx')%>%
  rename(date=Fecha,
         m2 = SF12722)%>%
  mutate(date = as_date(date, origin = "1899-12-30"))%>%
  filter(date <= '2014-03-01')

indpro <- read_excel('INDPRO.xlsx')%>%
  rename(date = observation_date,
         pibx = INDPRO)%>%
  mutate(date=as_date(date))%>%
  filter(date <= '2014-03-01')

inpc_core <- read_excel('INPC_core.xlsx')%>%
  rename(date = Fecha,
         core = SP74625)%>%
  #mutate(pi_core=((core/lag(core))-1)*100)%>%
  mutate(date = as_date(date, origin = "1899-12-30"))%>%
  filter(date <= '2014-03-01')


crb <- read_excel('CRB.xlsx')%>%
  rename(date = observation_date,
         crb = PALLFNFINDEXM)%>%
  mutate(date= as_date(date))%>%
  filter(date <= '2014-03-01')

cet28 <- read_excel('Cetes_28.xlsx')%>%
  rename(date= Fecha,
         cet28 = SF43936)%>%
  mutate(date = as_date(date))%>%
  mutate(mes_ref = floor_date(date, unit = "month")) %>%
  group_by(mes_ref) %>%
  summarise(int_rate = mean(cet28, na.rm = TRUE)) %>%
  rename(date=mes_ref)%>%
  filter(date <= '2014-03-01')%>%
  ungroup()

# Brecha del producto: Kalman MFK (Elizondo 2012 / Carrillo & Elizondo 2015)
# Estimada en Code/0b_brecha_kalman.R a partir del IGAE desest. (base 2018=100)
# y el PIB trimestral desest. (mmp 2018). Reemplaza la serie precalculada de Banxico.
pib_gap <- read_excel('brecha_prod/igae_kalman_gap.xlsx', sheet = "serie_var") %>%
  rename(igae = igae_kalman) %>%
  mutate(date = as.Date(date)) %>%          # POSIXct → Date para join correcto
  filter(date <= '2014-03-01')
cat(sprintf("pib_gap: %d obs | %s – %s | class: %s\n",
    nrow(pib_gap), format(min(pib_gap$date)), format(max(pib_gap$date)),
    class(pib_gap$date)))

inpp_viejo <- read_excel('INPP.xlsx')%>%
  rename(date = Fecha,
         inpp_viejo = SP6)%>%
  mutate(date = as_date(date, origin = "1899-12-30"))

inpp_nuevo <- read_excel('INPP_w_oil.xlsx')
inpp_nuevo <- as.data.frame(t(inpp_nuevo))
inpp_nuevo <- data.frame(
  date = rownames(inpp_nuevo), 
  inpp_nuevo = inpp_nuevo[, 1])
inpp_nuevo <- inpp_nuevo[-1, ]
inpp_nuevo <- inpp_nuevo%>%
  mutate(date = ym(date))

inpp_datos <- inpp_nuevo %>%
  full_join(inpp_viejo, by = "date")

inpp_datos <- inpp_datos %>%
  mutate(
    inpp_nuevo = as.numeric(inpp_nuevo),
    inpp_viejo = as.numeric(inpp_viejo)) %>%
  arrange(date)

mes_empalme <- ymd("2008-01-01")

valor_nuevo <- inpp_datos %>% filter(date == mes_empalme) %>% pull(inpp_nuevo)
valor_viejo <- inpp_datos %>% filter(date == mes_empalme) %>% pull(inpp_viejo)

factor_enlace <- valor_nuevo / valor_viejo

inpp_datos <- inpp_datos %>%
  mutate(
    inpp= case_when(
      date < mes_empalme  ~ inpp_viejo * factor_enlace,
      date >= mes_empalme ~ inpp_nuevo
    )
  )%>%
  dplyr::select(date,inpp)%>%
  filter(date <= '2014-03-01')
  

lista_dfs <- list(pib_gap, inpc_core, inpp_datos, cet28, m2, tcr, indpro, pcepilfe, crb)

# 2. Unimos iterativamente usando full_join 
base <- lista_dfs %>%
  reduce(full_join, by = "date") %>%
  arrange(date) %>%
  # Filtramos estrictamente para el periodo de la muestra del paper
  filter(date >= "2000-12-01" & date <= "2014-03-01")

cat(sprintf("base tras full_join: %d obs | %s – %s\n",
    nrow(base), format(min(base$date)), format(max(base$date))))
cat(sprintf("  NA en igae (pib_gap col): %d\n", sum(is.na(base$igae))))

wb <- createWorkbook()
addWorksheet(wb, "base")
writeData(wb, "base", base)
saveWorkbook(wb, file =paste0("raw_base.xlsx"), overwrite = TRUE)


base <- base %>% arrange(date)

# ==============================================================================
# PASO 0.5: AJUSTE ESTACIONAL (X-13ARIMA-SEATS)
# ==============================================================================
cat(sprintf("base: %d obs | %s – %s\n",
    nrow(base), format(min(base$date)), format(max(base$date))))

# start_date debe reflejar la fecha real de inicio de la serie para que
# X-13 asigne correctamente los factores estacionales por mes.
ajuste_x13 <- function(serie_vector, start_date) {
  serie_ts <- ts(serie_vector,
                 start     = c(year(start_date), month(start_date)),
                 frequency = 12)
  ajuste <- seasonal::seas(serie_ts)
  return(as.numeric(seasonal::final(ajuste)))
}

base_start <- min(base$date)

base_limpia <- base %>%
  mutate(
    core_sa = ajuste_x13(core, base_start),
    inpp_sa = ajuste_x13(inpp, base_start),
    m2_sa   = ajuste_x13(m2,   base_start)
  )

# Diagnóstico: verificar NAs introducidos por X-13
na_counts <- base_limpia %>%
  summarise(across(c(core_sa, inpp_sa, m2_sa), ~ sum(is.na(.))))
cat("NAs en series X-13:", na_counts$core_sa, "(core)", na_counts$inpp_sa, "(inpp)",
    na_counts$m2_sa, "(m2)\n")

# ==============================================================================
# PASO 1: CRECIMIENTO Y LOGARITMOS (Tasas anualizadas sobre series limpias)
# ==============================================================================
# Ahora utilizamos las series con sufijo '_sa' para calcular las variaciones
df_transform <- base_limpia %>%
  mutate(
    # Precios y Dinero (usando las series desestacionalizadas)
    pi_core = (log(core_sa) - log(lag(core_sa))) * 1200,
    pi_inpp = (log(inpp_sa) - log(lag(inpp_sa))) * 1200,
    g_m2    = (log(m2_sa) - log(lag(m2_sa))) * 1200,
    
    # Variables externas y financieras (se mantienen igual)
    pi_us   = (log(pce) - log(lag(pce))) * 1200,
    pi_crb  = (log(crb) - log(lag(crb))) * 1200,
    dep_tcr = (log(tcr) - log(lag(tcr))) * 1200,
    log_indpro = log(pibx)
  ) %>%
  drop_na()

cat(sprintf("df_transform después de drop_na: %d obs | %s – %s\n",
    nrow(df_transform), format(min(df_transform$date)), format(max(df_transform$date))))

# ==============================================================================
# PASO 2: BRECHA DEL PRODUCTO EXTERNO (EE. UU.)
# ==============================================================================
# Aplicamos el filtro HP a la serie en logaritmos con lambda = 129600

hp_us <- hpfilter(df_transform$log_indpro, freq = 129600, type = "lambda")
df_transform$y_us_gap <- (df_transform$log_indpro - hp_us$trend) * 100

# ==============================================================================
# PASO 3: TENDENCIA NOMINAL Y BRECHAS DOMÉSTICAS
# ==============================================================================
# 3.1 Extraemos la tendencia de la inflación subyacente (el "ancla" nominal)
hp_pi_core <- hpfilter(df_transform$pi_core, freq = 129600, type = "lambda")
df_transform$pi_trend <- hp_pi_core$trend

# 3.2 Construimos los gaps nominales restando esta tendencia común
df_final <- df_transform %>%
  mutate(
    pi_core_gap = pi_core - pi_trend,
    pi_inpp_gap = pi_inpp - pi_trend,
    m2_gap      = g_m2 - pi_trend,
    i_gap       = int_rate - pi_trend 
  ) %>%
  # Seleccionamos el vector de variables finales para el VAR
  dplyr::select(date, 
         igae,         # Brecha de producto México (Ya venía así de Banxico)
         pi_core_gap,  # Brecha de inflación subyacente
         pi_inpp_gap,  # Brecha de inflación al productor
         dep_tcr,      # Depreciación cambiaria
         i_gap,        # Brecha de tasa de interés
         m2_gap,       # Brecha de dinero
         y_us_gap,     # Brecha de producto EE. UU.
         pi_us,        # Inflación EE. UU.
         pi_crb)       # Inflación Materias Primas

# Verificamos los primeros renglones del panel listo para el modelo
cat(sprintf("df_final (database.xlsx): %d obs | %s – %s\n",
    nrow(df_final), format(min(df_final$date)), format(max(df_final$date))))
cat("Objetivo paper: 159 obs (Ene 2001 – Mar 2014)\n")
head(df_final)

wb <- createWorkbook()
addWorksheet(wb, "base")
writeData(wb, "base", df_final)
saveWorkbook(wb, file =paste0("database.xlsx"), overwrite = TRUE)
