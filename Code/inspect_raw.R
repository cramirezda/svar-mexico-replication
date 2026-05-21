library(readxl)

base <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras/raw data/brecha_prod/"

cat("=== IGAE indice ===\n")
igae <- read_excel(paste0(base, "igae_indice (2).xlsx"), sheet = 1)
cat("Dim:", nrow(igae), "x", ncol(igae), "\n")
cat("Columnas:", paste(names(igae), collapse = " | "), "\n")
print(head(igae, 10))
cat("...\n")
print(tail(igae, 5))

cat("\n=== PIB trimestral ===\n")
pib <- read_excel(paste0(base, "pibt_cte_valor (2).xlsx"), sheet = 1)
cat("Dim:", nrow(pib), "x", ncol(pib), "\n")
cat("Columnas:", paste(names(pib), collapse = " | "), "\n")
print(head(pib, 10))
cat("...\n")
print(tail(pib, 5))
