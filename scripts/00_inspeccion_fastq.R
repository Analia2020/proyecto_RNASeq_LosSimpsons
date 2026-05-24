# =====================================================================
# Script 01 - Inspección inicial de los archivos FASTQ
# Objetivo: verificar que R encuentra los 16 archivos y echar un
#           primer vistazo a su contenido.
# =====================================================================

# --- 1. Comprobar el directorio de trabajo---

getwd()

# --- 2. Listar los archivos FASTQ ---

fastq_files <- list.files(
  path       = "00_data/fastq",
  pattern    = "\\.fastq\\.gz$",   
  full.names = TRUE
)

# Mostrar la lista
fastq_files

length(fastq_files)  

# --- 4. Separar R1 y R2 para confirmar el emparejamiento ---
r1_files <- list.files("00_data/fastq", pattern = "_R1\\.fastq\\.gz$", full.names = TRUE)
r2_files <- list.files("00_data/fastq", pattern = "_R2\\.fastq\\.gz$", full.names = TRUE)

cat("Archivos R1:", length(r1_files), "\n")  
cat("Archivos R2:", length(r2_files), "\n") 

# --- 5. Extraer el nombre  y conocer su tamaño---
# basename() quita la ruta y deja solo el nombre del archivo
sample_names <- sub("_R1\\.fastq\\.gz$", "", basename(r1_files))
sample_names

fastq_files <- list.files("00_data/fastq", pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Pesos en megabytes
pesos_mb <- file.info(fastq_files)$size / (1024 * 1024)
data.frame(archivo = basename(fastq_files), peso_MB = round(pesos_mb, 1))

# Peso total
cat("Peso total:", round(sum(pesos_mb), 1), "MB\n")