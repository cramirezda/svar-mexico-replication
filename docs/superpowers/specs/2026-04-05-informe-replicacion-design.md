# Diseño: Informe de Replicación SVAR México

**Fecha:** 2026-04-05  
**Proyecto:** Replicación Carrillo & Elizondo (2015/2018) — México  
**Estado:** Aprobado por usuario

---

## Objetivo

Generar un informe HTML autocontenido que compare, sección por sección, los resultados del paper original de Carrillo & Elizondo (2015/2018) con los de nuestra replicación para México. El informe debe ser portable (un solo archivo `.html`) y visualmente claro para comparabilidad académica.

---

## Arquitectura

### Entradas
- `raw data/raw_base.xlsx` — series intermedias antes de extraer gaps (para series originales + tendencias HP)
- `raw data/database.xlsx` — gaps y tasas finales (para estadísticos descriptivos y diagnósticos)
- `raw data/database_for_matlab.csv` — base de datos procesada (STATA/MATLAB)
- `figures/fig7_exclusion.png` — IRFs Cholesky (MATLAB)
- `figures/fig8_sign_restrictions.png` — Restricciones de signo (MATLAB)
- `figures/fig9_comparison.png` — Comparación estrategias (MATLAB)
- `figures/varstable.png` — Raíces del VAR (STATA)
- Resultados cuantitativos del paper (extraídos del PDF y hard-coded en el script)

### Script generador
**`Code/4Informe_replicacion.R`**

Flujo interno:
1. Leer `raw_base.xlsx` con `readxl` para series originales (pre-transformación)
2. Aplicar filtro HP (λ=129,600) sobre series originales para obtener tendencia y brecha
3. Leer `database.xlsx` para gaps finales (estadísticos descriptivos, diagnósticos)
3. Generar 6 figuras nuevas con `ggplot2` y guardarlas en `figures/`
4. Leer y codificar en base64 todas las imágenes (nuevas + MATLAB + STATA)
5. Ensamblar el HTML completo con `glue` / `cat()`
6. Escribir `informe_replicacion.html` en la raíz del proyecto

### Salida
**`informe_replicacion.html`** — archivo único portable, todas las imágenes embebidas como base64 data URIs.

---

## Figuras nuevas a generar en R

| Nombre archivo | Descripción | Herramienta |
|---|---|---|
| `fig_series_originales.png` | 9 series en niveles con línea de tendencia HP superpuesta. Grilla 3×3. | ggplot2 + mFilter |
| `fig_series_gaps.png` | 9 series transformadas (gaps/tasas). Grilla 3×3. | ggplot2 |
| `tabla_descriptivos.png` | Media, SD, mín, máx, asimetría para 9 variables. | gt o kable |
| `tabla_adf_pp.png` | ADF (3 rezagos) y PP (4 rezagos): estadístico, p-valor, decisión. | gt |
| `tabla_varsoc.png` | Criterios AIC, SIC, LR para p=1..6. Resalta óptimo. | gt |
| `tabla_lm_norm.png` | LM autocorrelación (lags 1–12) + Jarque-Bera por ecuación. | gt |

---

## Estructura del Informe HTML

### Layout general
- **Header fijo:** título del paper, autores, período, frecuencia
- **Scorecard semáforo:** píldoras de color (✓ verde / ~ amarillo / ✗ rojo) para diagnóstico rápido de concordancia
- **Barra de navegación:** anclas a cada sección
- **Cuerpo:** dos columnas permanentes (📄 Paper C&E izq. | 🔬 Replicación der.)
- **Secciones:** separadas por encabezados de color de ancho completo

### Secciones (8)

#### §0 — Scorecard semáforo
Píldoras de concordancia generadas dinámicamente:
- ✓ Estacionariedad, Rezagos p=3, Estabilidad VAR, Signos IRF Cholesky, Price puzzle eliminado
- ~ β_lower (−0.715 vs −0.7 del paper)
- ✗ T observaciones (158 vs 159 del paper)

#### §1 — Descripción de Datos
- Tabla: muestra, variables, fuentes, transformaciones (papel vs. replicación)
- Figuras: series originales con tendencia HP (izq: referencia del paper | der: `fig_series_originales.png`)
- Nota: diferencia en TCR (49 socios comerciales vs. 21 del paper)

#### §2 — Estadísticas Descriptivas
- Tabla comparativa: estadísticos descriptivos del paper (Tabla A1 apéndice) vs. `tabla_descriptivos.png`

#### §3 — Diagnósticos Estadísticos
- ADF / Phillips-Perron: resultados cualitativos del paper vs. `tabla_adf_pp.png`
- Selección de rezagos: descripción del paper vs. `tabla_varsoc.png`
- Estabilidad VAR: descripción vs. `varstable.png` (STATA)
- Autocorrelación y normalidad: descripción vs. `tabla_lm_norm.png`

#### §4 — Estimación de β_lower
- Ecuación 11 del paper, valor puntual y bound del paper vs. estimación propia
- β_minus = −0.400, SE = 0.192, β_lower = −0.715 (paper: β̂ = −0.035, β_lower = −0.7)

#### §5 — Restricciones de Exclusión (Figura 7)
- Figura 7 del paper (imagen del PDF) vs. `fig7_exclusion.png` (MATLAB)
- Tabla resumen: impactos h=0 para las 6 variables (paper vs. replicación)

#### §6 — Restricciones de Signo (Figura 8)
- Figura 8 del paper vs. `fig8_sign_restrictions.png` (MATLAB)
- Modelos aceptados: paper (N/D) vs. replicación (198,654 std / 66,874 aug)

#### §7 — Comparación Final (Figura 9)
- Figura 9 del paper vs. `fig9_comparison.png` (MATLAB)
- Tabla resumen: impactos h=0 de las 3 estrategias para las 6 variables

---

## Paquetes R requeridos

```r
readxl, dplyr, tidyr, ggplot2, patchwork, mFilter,
gt, base64enc, glue, scales, stringr, tseries, FinTS
```

---

## Resultados cuantitativos del paper a extraer (Sección 4 + Apéndice)

Se leerá el PDF `revision literatura/Carrillo&Elizondo 2015_18.pdf` para extraer:
- Tabla 1: descripción de variables y fuentes
- Tabla A1: estadísticos descriptivos
- Tabla A2: resultados de pruebas de raíz unitaria
- Tabla A3 / texto: criterios de selección de rezagos
- Valores de β_lower (Ecuación 11)
- Impactos h=0 de las IRFs reportadas en texto

---

## Convenciones de diseño CSS

- Fondo columna Paper: `#eaf4fb` (azul muy claro) con borde izquierdo `#2980b9`
- Fondo columna Replicación: `#e8f8f5` (verde muy claro) con borde izquierdo `#16a085`
- Encabezados de sección: ancho completo, fondo oscuro, texto blanco
- Semáforo: verde `#27ae60`, amarillo `#f39c12`, rojo `#e74c3c`
- Tipografía: `font-family: 'Georgia', serif` para cuerpo; `monospace` para valores numéricos

---

## Archivos que se crean/modifican

| Archivo | Acción |
|---|---|
| `Code/4Informe_replicacion.R` | CREAR — script generador |
| `informe_replicacion.html` | CREAR — salida del script |
| `figures/fig_series_originales.png` | CREAR — generada por el script |
| `figures/fig_series_gaps.png` | CREAR — generada por el script |
| `figures/tabla_descriptivos.png` | CREAR — generada por el script |
| `figures/tabla_adf_pp.png` | CREAR — generada por el script |
| `figures/tabla_varsoc.png` | CREAR — generada por el script |
| `figures/tabla_lm_norm.png` | CREAR — generada por el script |
