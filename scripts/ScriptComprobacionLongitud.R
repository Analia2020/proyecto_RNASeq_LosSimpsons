# Contar lecturas de R1 y R2 
library(ShortRead)

r1_marge <- "00_data/fastq/MargeSimpson_R1.fastq.gz"
r2_marge <- "00_data/fastq/MargeSimpson_R2.fastq.gz"

cat("Lecturas en R1:", length(readFastq(r1_marge)), "\n")
cat("Lecturas en R2:", length(readFastq(r2_marge)), "\n")

r1_files <- list.files("00_data/fastq", pattern = "_R1\\.fastq\\.gz$", full.names = TRUE)
r2_files <- list.files("00_data/fastq", pattern = "_R2\\.fastq\\.gz$", full.names = TRUE)
sample_names <- sub("_R1\\.fastq\\.gz$", "", basename(r1_files))

cat("Comprobando emparejamiento R1 vs R2...\n\n")

for (i in seq_along(sample_names)) {
  n_r1 <- length(readFastq(r1_files[i]))
  n_r2 <- length(readFastq(r2_files[i]))
  estado <- ifelse(n_r1 == n_r2, "OK", "MISMATCH")
  cat(sprintf("  %-20s R1=%d   R2=%d   %s\n",
              sample_names[i], n_r1, n_r2, estado))
}



# Comparar el orden de los IDs en R1 y R2 de Marge

id_r1 <- sub("/[12]$", "", as.character(id(r1)))
id_r2 <- sub("/[12]$", "", as.character(id(r2)))

cat("Longitud R1:", length(id_r1), "\n")
cat("Longitud R2:", length(id_r2), "\n\n")

# Comparar los primeros 10 de cada uno
cat("=== Primeros 10 IDs de R1 ===\n")
print(head(id_r1, 10))
cat("\n=== Primeros 10 IDs de R2 ===\n")
print(head(id_r2, 10))

# ¿Hasta qué posición coinciden?
min_len <- min(length(id_r1), length(id_r2))
coinciden <- id_r1[1:min_len] == id_r2[1:min_len]
cat("\nCoinciden las primeras", sum(coinciden), "lecturas en orden\n")
cat("Primera discrepancia en posición:", which(!coinciden)[1], "\n")


# Mirar alrededor de la posición 289 en ambos archivos
cat("=== R1 alrededor de la posición 289 ===\n")
print(id_r1[285:295])

cat("\n=== R2 alrededor de la posición 289 ===\n")
print(id_r2[285:295])

# También mirar las últimas lecturas de cada archivo
cat("\n=== Últimas 5 lecturas de R1 ===\n")
print(tail(id_r1, 5))

cat("\n=== Últimas 5 lecturas de R2 ===\n")
print(tail(id_r2, 5))