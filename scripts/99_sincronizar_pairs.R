# --- BLOQUE 7: Alineamiento de los FASTQ contra el índice ---

# Carpeta de salida para los BAM
bam_dir <- "02_alineamiento/bam"
dir.create(bam_dir, showWarnings = FALSE)

# Listar los archivos R1 y R2 ordenados
r1_files <- list.files("00_data/fastq", pattern = "_R1\\.fastq\\.gz$", full.names = TRUE)
r2_files <- list.files("00_data/fastq", pattern = "_R2\\.fastq\\.gz$", full.names = TRUE)

# Nombres de muestra 
sample_names <- sub("_R1\\.fastq\\.gz$", "", basename(r1_files))

# Verificación
cat("Muestras a alinear:", length(sample_names), "\n")
print(sample_names)

# Rutas de salida para los BAM 
bam_files <- file.path(bam_dir, paste0(sample_names, ".bam"))

# Lanzar el alineamiento de las 8 muestras 
align(
  index             = indice_prefijo,
  readfile1         = r1_files,
  readfile2         = r2_files,
  output_file       = bam_files,
  type              = "rna",
  nthreads          = 2,
  unique            = FALSE,
  nBestLocations    = 1
)

cat("\n✔ Alineamiento terminado. Archivos BAM en:", bam_dir, "\n")
list.files(bam_dir)