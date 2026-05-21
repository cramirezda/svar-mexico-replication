# ==============================================================================
#  BLOQUE 2: CHOLESKY IRFs — FIGURA 7
#  Replicación: Carrillo & Elizondo (2015) — México
#  Equivalente al Bloque 2 del script 3SVAR_Mexico.do
# ==============================================================================

library(vars)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(readxl)

# ── Directorios ────────────────────────────────────────────────────────────────
wdir  <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras"
ddir  <- file.path(wdir, "raw data")
fdir  <- file.path(wdir, "figures")
dir.create(fdir, showWarnings = FALSE)

# ── Parámetros (igual que el .do) ─────────────────────────────────────────────
p      <- 3      # rezagos
h      <- 16     # horizonte IRF
nboot  <- 2000   # réplicas bootstrap
ci     <- 68     # nivel IC en %
z_val  <- 1.0    # ±1 SD ≈ 68%
set.seed(42)

# ==============================================================================
#  SECCIÓN 1: CARGA DE DATOS
# ==============================================================================
db <- read_excel(file.path(ddir, "database.xlsx"), sheet = "base")

# Convertir fecha
db$date <- as.Date(db$date)
db <- db[db$date >= as.Date("2001-01-01") & db$date <= as.Date("2014-03-31"), ]

endo_vars <- c("igae", "pi_core_gap", "pi_inpp_gap", "dep_tcr", "i_gap", "m2_gap")
exo_vars  <- c("y_us_gap", "pi_us", "pi_crb")

Y <- as.matrix(db[, endo_vars])
Z <- as.matrix(db[, exo_vars])

cat("Observaciones en muestra:", nrow(Y), "\n")

# ==============================================================================
#  SECCIÓN 2: FUNCIÓN DE BOOTSTRAP CHOLESKY
#  Replica exactamente la lógica de irf create ... bs reps() level()
# ==============================================================================

bootstrap_cholesky_irf <- function(Y, Z = NULL, p, impulse_idx, h, nboot, z_val,
                                   use_exo = TRUE, verbose = TRUE) {
  # ── Construir regresores ────────────────────────────────────────────────────
  T_total <- nrow(Y)
  n       <- ncol(Y)
  
  build_X <- function(Y, Z, p, use_exo) {
    T  <- nrow(Y)
    Yp <- Y[(p + 1):T, , drop = FALSE]        # variable dependiente
    X  <- matrix(1, nrow = T - p, ncol = 1)   # constante
    for (lag in 1:p) {
      X <- cbind(X, Y[(p + 1 - lag):(T - lag), , drop = FALSE])
    }
    if (use_exo && !is.null(Z)) {
      X <- cbind(X, Z[(p + 1):T, , drop = FALSE])
    }
    list(Yp = Yp, X = X)
  }
  
  dat   <- build_X(Y, Z, p, use_exo)
  Yp    <- dat$Yp
  X     <- dat$X
  T_eff <- nrow(Yp)
  
  # ── OLS ────────────────────────────────────────────────────────────────────
  B     <- solve(crossprod(X), crossprod(X, Yp))   # k x n
  resid <- Yp - X %*% B                             # T_eff x n
  Sigma <- crossprod(resid) / (T_eff - ncol(X))
  
  # ── Función: propagar IRF desde choque de impacto ──────────────────────────
  propagate_irf <- function(B, shock_vec, n, p, h) {
    # Extraer bloques A1,...,Ap (B tiene forma [const | A1 | A2 | ... | Ap | (exo)])
    # primera fila = constante (índice 1 en R)
    irf_mat <- matrix(0, nrow = h + 1, ncol = n)
    irf_mat[1, ] <- shock_vec
    for (hh in 2:(h + 1)) {
      acc <- numeric(n)
      for (j in 1:p) {
        if (hh - j >= 1) {
          row_start <- 1 + (j - 1) * n + 1   # +1 por constante, base-1 en R
          row_end   <- 1 + j * n
          Aj        <- t(B[row_start:row_end, ])   # n x n
          acc       <- acc + Aj %*% irf_mat[hh - j, ]
        }
      }
      irf_mat[hh, ] <- acc
    }
    irf_mat   # (h+1) x n
  }
  
  # ── Cholesky del VAR base ───────────────────────────────────────────────────
  P0    <- t(chol(Sigma))                    # triangular inferior
  shock <- P0[, impulse_idx]                 # columna j del choque
  
  # IRF puntual (sin bootstrap)
  irf_point <- propagate_irf(B, shock, n, p, h)
  
  # ── Bootstrap de residuos (Runkle 1987) ────────────────────────────────────
  irf_boot <- array(0, dim = c(nboot, h + 1, n))
  
  if (verbose) cat("  Bootstrap progress: ")
  for (b in 1:nboot) {
    if (verbose && b %% 500 == 0) cat(b, "")
    
    # Remuestrear residuos con reemplazo
    idx     <- sample(T_eff, T_eff, replace = TRUE)
    res_b   <- resid[idx, ]
    Sigma_b <- crossprod(res_b) / T_eff
    Sigma_b <- (Sigma_b + t(Sigma_b)) / 2   # simetrizar
    
    # Cholesky con jitter para estabilidad numérica
    jitter  <- 1e-10 * sum(diag(Sigma_b)) / n
    P_b     <- tryCatch(
      t(chol(Sigma_b + jitter * diag(n))),
      error = function(e) P0
    )
    
    shock_b <- P_b[, impulse_idx]
    irf_b   <- propagate_irf(B, shock_b, n, p, h)
    irf_boot[b, , ] <- irf_b
  }
  if (verbose) cat("\n")
  
  # ── Bandas: media ± z * SD (equivale a ±1 SD para 68%) ───────────────────
  irf_sd  <- apply(irf_boot, c(2, 3), sd)
  irf_lo  <- irf_point - z_val * irf_sd
  irf_hi  <- irf_point + z_val * irf_sd
  
  list(
    med  = irf_point,   # (h+1) x n
    lo   = irf_lo,
    hi   = irf_hi,
    boot = irf_boot
  )
}

# ==============================================================================
#  SECCIÓN 3: ESTIMAR AMBOS MODELOS
# ==============================================================================
cat("\nEstimando VAR simple (sin exógenas)...\n")
irf_simple <- bootstrap_cholesky_irf(
  Y, Z, p,
  impulse_idx = 5,   # i_gap es la variable 5
  h, nboot, z_val,
  use_exo = FALSE,
  verbose  = TRUE
)

cat("\nEstimando VAR con exógenas...\n")
irf_exo <- bootstrap_cholesky_irf(
  Y, Z, p,
  impulse_idx = 5,
  h, nboot, z_val,
  use_exo = TRUE,
  verbose  = TRUE
)

# ==============================================================================
#  SECCIÓN 4: CONSTRUIR DATA FRAMES PARA GGPLOT
# ==============================================================================
var_labels <- c(
  igae        = "Output",
  pi_core_gap = "Inflation (core)",
  pi_inpp_gap = "PPI",
  dep_tcr     = "Real ER dep",
  i_gap       = "Nom. interest rate",
  m2_gap      = "Money growth"
)

build_irf_df <- function(irf_obj, model_name, endo_vars) {
  steps <- 0:h
  n     <- length(endo_vars)
  rows  <- list()
  for (v in 1:n) {
    rows[[v]] <- data.frame(
      step     = steps,
      variable = endo_vars[v],
      model    = model_name,
      med      = irf_obj$med[, v],
      lo       = irf_obj$lo[, v],
      hi       = irf_obj$hi[, v]
    )
  }
  bind_rows(rows)
}

df_exo    <- build_irf_df(irf_exo,    "VAR w/ exo",  endo_vars)
df_simple <- build_irf_df(irf_simple, "Simple VAR",  endo_vars)
df_all    <- bind_rows(df_exo, df_simple)

# Etiquetas para los paneles
df_all$var_label <- factor(
  var_labels[df_all$variable],
  levels = var_labels
)

# ==============================================================================
#  SECCIÓN 5: FIGURA 7 — GRILLA 2×3
# ==============================================================================
plot_irf_panel <- function(var_name, df) {
  
  d_exo <- df %>% filter(variable == var_name, model == "VAR w/ exo")
  d_smp <- df %>% filter(variable == var_name, model == "Simple VAR")
  
  ggplot() +
    # Banda VAR w/ exo (navy)
    geom_ribbon(data = d_exo,
                aes(x = step, ymin = lo, ymax = hi),
                fill = "#1f4e79", alpha = 0.25) +
    # Línea VAR w/ exo
    geom_line(data = d_exo,
              aes(x = step, y = med),
              color = "#1f4e79", linewidth = 0.7) +
    # Banda Simple VAR (gris)
    geom_ribbon(data = d_smp,
                aes(x = step, ymin = lo, ymax = hi),
                fill = "#888888", alpha = 0.20) +
    # Línea Simple VAR (punteada)
    geom_line(data = d_smp,
              aes(x = step, y = med),
              color = "#888888", linewidth = 0.6, linetype = "dashed") +
    # Línea cero
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "black", linewidth = 0.4) +
    scale_x_continuous(breaks = seq(0, h, by = 4)) +
    labs(
      title = var_labels[var_name],
      x     = "Months after shock",
      y     = "Pct pts from trend"
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 7),
      axis.text    = element_text(size = 7),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    )
}

# Generar los 6 paneles
panels <- lapply(endo_vars, plot_irf_panel, df = df_all)
names(panels) <- endo_vars

# Combinar con patchwork (equivalente a graph combine)
fig7 <- wrap_plots(panels, nrow = 2, ncol = 3) +
  plot_annotation(
    title   = "Figure 7: Exclusion Restrictions — Mexico",
    caption = "68% bands (±1 SD bootstrap). Navy = VAR w/ exo, Gray dashed = Simple VAR.",
    theme   = theme(
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 7,  hjust = 0)
    )
  )

# ── Exportar ───────────────────────────────────────────────────────────────────
ggsave(file.path(fdir, "fig7_exclusion.pdf"), fig7,
       width = 10, height = 6, units = "in")
ggsave(file.path(fdir, "fig7_exclusion.png"), fig7,
       width = 10, height = 6, units = "in", dpi = 200)

cat("→ Figura 7 guardada en:", fdir, "\n")