# =====================================================================
# Script 01 - Control de calidad de FASTQ
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
# Autor:      Analia Pastrana
# Fecha:      2026-05-23
#
# Input:      00_data/fastq/*.fastq.gz (16 archivos)
# Output:     01_qc/rqc_report.html (informe agregado)
#             01_qc/qa_object.rds   (objeto R con los datos del QC)
#
# Criterio:   1.1 - Control de calidad de los ficheros (5%)
# =====================================================================

# Cargamos el paquete
#BiocManager::install("Rqc")
library(Rqc)


# --- Localizamos los FASTQ ---
fastq_files <- list.files(
  path       = "00_data/fastq",
  pattern    = "\\.fastq\\.gz$",
  full.names = TRUE
)

# Verificación rápida: deben ser 16
length(fastq_files)
fastq_files

# --- Ejecutamos el análisis de calidad ---
# rqcQA() lee cada FASTQ, muestrea un subconjunto de lecturas
# (por defecto las primeras 1 millón, en nuestro caso será todas
# porque tenemos muchas menos) y calcula las métricas de calidad.
#
# Argumentos:
#   - x:         vector con las rutas a los FASTQ
#   - workers:   nº de núcleos del CPU a usar (1 = secuencial, seguro)
#   - sample:    si muestrear o usar todas las lecturas
#   - pair:      vector que indica qué archivos son pareja
#                (mismo número = misma pareja)

#Creamos el vector pair: los R1 y R2 del mismo personaje
# deben tener el mismo número (1, 1, 2, 2, 3, 3, ...)
pair_vector <- rep(1:8, each = 2)
pair_vector

# Lanzamos el QC
qa <- rqcQA(
  x       = fastq_files,
  workers = 1,
  sample  = FALSE,    # usamos todas las lecturas
  pair    = pair_vector
)


# --- Crear la carpeta de salida ---
dir.create("01_qc", showWarnings = FALSE)

# --- Generar informe HTML interactivo ---
# rqcReport() en versiones recientes solo acepta:
#   - el objeto qa
#   - templateFile (opcional)
# El archivo se genera en una carpeta temporal y devuelve la ruta.
report_path <- rqcReport(qa)

cat("\n✔ Informe generado en:", report_path, "\n")

# --- Lo copiamos a nuestra carpeta 01_qc/ para tenerlo organizado ---
file.copy(
  from      = report_path,
  to        = "01_qc/rqc_report.html",
  overwrite = TRUE
)

# Y también copiamos la carpeta de assets (CSS, JS, imágenes) que necesita
# el HTML para verse bien. Está al lado del archivo HTML, en una carpeta
# con el mismo nombre pero terminada en "_files".
report_dir  <- dirname(report_path)
assets_name <- sub("\\.html$", "_files", basename(report_path))
assets_path <- file.path(report_dir, assets_name)

if (dir.exists(assets_path)) {
  file.copy(
    from      = assets_path,
    to        = "01_qc/",
    recursive = TRUE,
    overwrite = TRUE
  )
  # Renombramos para que coincida con el HTML copiado
  file.rename(
    from = file.path("01_qc", assets_name),
    to   = "01_qc/rqc_report_files"
  )
}

cat("✔ Informe copiado a: 01_qc/rqc_report.html\n")

