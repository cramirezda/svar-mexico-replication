## diag_hump.R — diagnóstico rápido del hump shape en output gap IRF
## Calcula A1[1,:] y P diagonal para entender el mecanismo del hump

library(readxl); library(vars); library(mFilter); library(dplyr)

wdir <- "C:/Users/carlo/OneDrive - INSTITUTO TECNOLOGICO AUTONOMO DE MEXICO/Escritorio/varras"
db   <- read_excel(file.path(wdir, "raw data/database.xlsx"), sheet = "base")

db$date <- as.Date(db$date)
db <- db %>% filter(date >= as.Date("2001-01-01"), date <= as.Date("2014-03-01"))

# Variables endógenas y exógenas
endo_vars <- c("igae","pi_core_gap","pi_inpp_gap","dep_tcr","i_gap","m2_gap")
exo_vars  <- c("y_us_gap","pi_us","pi_crb")

for (v in c(endo_vars, exo_vars)) db[[v]] <- as.numeric(db[[v]])

Y <- as.matrix(db[, endo_vars])
Z <- as.matrix(db[, exo_vars])
T <- nrow(Y); n <- ncol(Y); p <- 3

cat("=== ESTADÍSTICAS BÁSICAS ===\n")
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "igae",        sd(Y[,1]), min(Y[,1]), max(Y[,1])))
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "pi_core_gap", sd(Y[,2]), min(Y[,2]), max(Y[,2])))
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "pi_inpp_gap", sd(Y[,3]), min(Y[,3]), max(Y[,3])))
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "dep_tcr",     sd(Y[,4]), min(Y[,4]), max(Y[,4])))
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "i_gap",       sd(Y[,5]), min(Y[,5]), max(Y[,5])))
cat(sprintf("%-15s  sd=%6.4f  min=%7.4f  max=%7.4f\n",
            "m2_gap",      sd(Y[,6]), min(Y[,6]), max(Y[,6])))

# ── OLS VAR(3) con exógenas (manual, igual al MATLAB) ─────────────────────────
Yp <- Y[(p+1):T, ]
Ylags <- do.call(cbind, lapply(1:p, function(j) Y[(p+1-j):(T-j), ]))
Xp <- cbind(1, Ylags, Z[(p+1):T, ])

B <- solve(t(Xp) %*% Xp) %*% t(Xp) %*% Yp   # k x n
resid <- Yp - Xp %*% B
Teff  <- nrow(Yp); k <- ncol(Xp)
Sigma <- (t(resid) %*% resid) / (Teff - k)
P     <- t(chol(Sigma))   # lower-triangular Cholesky

cat("\n=== CHOLESKY P — DIAGONAL (SDs de innovaciones) ===\n")
cat(sprintf("P[igae,igae]        = %.4f\n", P[1,1]))
cat(sprintf("P[pi_core,pi_core]  = %.4f\n", P[2,2]))
cat(sprintf("P[pi_inpp,pi_inpp]  = %.4f\n", P[3,3]))
cat(sprintf("P[dep_tcr,dep_tcr]  = %.4f\n", P[4,4]))
cat(sprintf("P[i_gap,i_gap]      = %.4f\n", P[5,5]))
cat(sprintf("P[m2_gap,m2_gap]    = %.4f\n", P[6,6]))
cat(sprintf("\nRatio P[dep_tcr]/P[i_gap] = %.2f  (debe ser ~1 para balance)\n",
            P[4,4]/P[5,5]))

# ── Extrae A1 (lag-1 coefficients) ────────────────────────────────────────────
# B es k×n; filas 2:(1+n) son los coeficientes del lag 1
# Acoefs{1} = B(2:7,:)' en MATLAB = t(B[2:7,]) en R
A1 <- t(B[2:(n+1), ])   # n×n

cat("\n=== A1 — PRIMERA FILA (dinámica del output gap) ===\n")
var_names <- endo_vars
cat(sprintf("A1[output, %-15s] = %8.4f\n", var_names, A1[1,]))

cat("\n=== MECANISMO DEL HUMP (h=1) ===\n")
cat("resp(h=0) mediana en sign restrictions (del MATLAB Section 7):\n")
resp0 <- c(-0.1638, -0.2168, -0.8914, -6.2411, 0.1306, -1.8144)
names(resp0) <- var_names

h1_contrib <- A1[1,] * resp0
cat(sprintf("  A1[out,%-12s] * resp0 = %8.4f * %8.4f = %8.4f\n",
            var_names, A1[1,], resp0, h1_contrib))
cat(sprintf("\n  Suma (resp output h=1 predicho): %.4f\n", sum(h1_contrib)))
cat(sprintf("  → Canal dep_tcr (%s): %.4f\n",
            ifelse(h1_contrib[4] < 0, "EMPUJA ABAJO ✓ hump", "empuja arriba"),
            h1_contrib[4]))
cat(sprintf("  → Canal m2_gap  (%s): %.4f\n",
            ifelse(h1_contrib[6] < 0, "EMPUJA ABAJO ✓ hump", "empuja arriba"),
            h1_contrib[6]))

cat("\n=== A2 y A3 — PRIMERA FILA ===\n")
A2 <- t(B[(n+2):(2*n+1), ])
A3 <- t(B[(2*n+2):(3*n+1), ])
cat(sprintf("A2[output, %-15s] = %8.4f\n", var_names, A2[1,]))
cat(sprintf("A3[output, %-15s] = %8.4f\n", var_names, A3[1,]))

cat("\n=== IRF ANALÍTICA COMPLETA (usando resp(h=0) mediana) ===\n")
H <- 24
Acoefs <- list(A1, A2, A3)

irf_analytical <- matrix(0, nrow = H+1, ncol = n)
irf_analytical[1, ] <- resp0

for (hh in 2:(H+1)) {
  irf_h <- numeric(n)
  for (j in 1:p) {
    if (hh - j >= 1) {
      irf_h <- irf_h + Acoefs[[j]] %*% irf_analytical[hh-j, ]
    }
  }
  irf_analytical[hh, ] <- irf_h
}

cat("h   output    pi_core   pi_inpp   dep_tcr   i_gap     m2_gap\n")
for (hh in 1:min(25, H+1)) {
  cat(sprintf("%2d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
              hh-1, irf_analytical[hh,1], irf_analytical[hh,2],
              irf_analytical[hh,3], irf_analytical[hh,4],
              irf_analytical[hh,5], irf_analytical[hh,6]))
}

hump_h <- which.min(irf_analytical[, 1])
cat(sprintf("\nOutput mínimo (más negativo) en h=%d: %.4f\n",
            hump_h - 1, irf_analytical[hump_h, 1]))
cat(sprintf("h=0: %.4f  h=1: %.4f  → ¿hump? %s\n",
            irf_analytical[1,1], irf_analytical[2,1],
            ifelse(abs(irf_analytical[2,1]) > abs(irf_analytical[1,1]),
                   "SÍ (h=1 más negativo)", "NO (h=1 menos negativo)")))
