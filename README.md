# Análisis de Expresión Diferencial de Genes Relacionados con la Obesidad mediante RNA-seq

**Comparativa: Obeso2 vs Normopeso**

Asignatura: Secuenciación y Ómicas de Próxima Generación  
Universidad Internacional de La Rioja (UNIR)  
Fecha: Mayo 2026

---

## Descripción del proyecto

Análisis de expresión diferencial mediante RNA-seq sobre un panel de 37 genes relacionados con la obesidad. Se comparan dos grupos fenotípicos (Obeso2: Marge, Patty, Selma vs Normopeso: Bart, Lisa, Maggie) utilizando datos simulados de 8 personajes del universo de Los Simpson.

## Hipótesis y objetivo

**Hipótesis**: el grupo Obeso2 presenta un perfil de expresión distinguible del Normopeso, con desregulación de genes clave del control del apetito.

**Objetivo**: identificar los genes diferencialmente expresados entre ambos grupos y discutir su relevancia biológica en obesidad.

---

## Estructura del proyecto
proyecto_rnaseq_simpsons/
├── 00_data/
│   ├── fastq/                          # 16 archivos FASTQ paired-end
│   └── referencia/
│       ├── tx2gene.tsv                 # mapeo transcrito → gen (706 transcritos)
│       └── carpetas_genes/             # 37 carpetas con secuencias rna.fna
│
├── 01_qc/                              # Control de calidad
│   ├── rqc_report.html
│   └── qa_object.rds
│
├── 02_alineamiento/
│   ├── referencia_combinada.fa         # FASTA con 706 transcritos
│   ├── indice/                         # índice de Rsubread
│   └── bam/                            # 8 archivos BAM
│
├── 03_resultados/
│   ├── cuentas_genes.csv
│   ├── featurecounts_resultado.rds
│   ├── deseq2_resultados.csv
│   ├── deseq2_DEGs.csv
│   ├── dds_object.rds
│   ├── res_object.rds
│   ├── figuras/
│   │   ├── volcano_plot.png
│   │   └── heatmap_DEGs_obeso2_vs_normopeso.png
│   └── tablas/
│       └── tabla_DEGs_completa.xlsx
│
└── scripts/
├── 00_diagnostico_r1r2.R           # diagnóstico inicial de FASTQ
├── 00_inspeccion_fastq.R           # vistazo inicial a los datos
├── 01_control_de_calidad.R         # QC con Rqc
├── 02_construir_referencia.R       # FASTA combinado + índice Rsubread
├── 03_alineamiento_singleend.R     # alineamiento con Rsubread::align
├── 04_cuantificacion.R             # featureCounts
├── 05_deseq2.R                     # análisis de expresión diferencial
├── 06_visualizacion_volcano.R      # volcano plot
├── 07_heatmap.R                    # heatmap de DEGs
├── 08_tabla_expresion.R            # tabla resumen de DEGs
└── 99_sincronizar_pairs.R          # script auxiliar (no usado en pipeline final)

---

## Pipeline

### 1. Control de calidad
- **Herramienta**: `Rqc` (Bioconductor)
- **Resultado**: calidad media Q≈35, longitud uniforme 151 pb, composición de bases balanceada
- **Decisión**: no se realizó trimming

### 2. Construcción de la referencia
- 37 archivos `rna.fna` (formato NCBI Datasets) combinados en un único FASTA
- 706 transcritos totales unidos en `referencia_combinada.fa`
- Indexación con `Rsubread::buildindex`

### 3. Alineamiento
- **Modo**: single-end (solo R1)
- **Justificación**: se detectaron inconsistencias en el emparejamiento R1/R2 de las 3 muestras del grupo Obeso2 (24 lecturas adicionales por archivo R2, con IDs base no coincidentes a partir de la posición 289). Se procesaron todas las muestras como single-end para mantener homogeneidad metodológica.
- **Herramienta**: `Rsubread::align`
- **Resultado**: mapeo ≥99.7% en las 8 muestras

### 4. Cuantificación
- **Herramienta**: `featureCounts` (paquete Rsubread)
- **Anotación**: SAF construida a partir del TSV transcrito→gen
- **Parámetros clave**: `countMultiMappingReads=TRUE`, `fraction=TRUE` (lecturas multi-mapping fraccionadas entre transcritos)
- **Resultado**: matriz 37 genes × 8 muestras, asignación ~100%

### 5. Análisis de expresión diferencial
- **Herramienta**: `DESeq2`
- **Diseño**: `~ grupo`
- **Contraste**: `c("grupo", "Obeso2", "Normopeso")`
- **Modificación**: estimación gen a gen de la dispersión (`estimateDispersionsGeneEst`) seguida del test de Wald, debido al reducido número de genes (37) que impide el ajuste paramétrico estándar
- **Criterios DEG**: padj < 0.05 (Benjamini-Hochberg) y |log2FoldChange| > 0.58 (cambio ≥1.5×)

### 6. Visualización
- **Volcano plot**: ggplot2 + ggrepel
- **Heatmap**: pheatmap (Z-score sobre VST, clustering jerárquico)
- **Tabla**: openxlsx (formato Excel con colores)

---

## Resultados principales

Se identificaron **8 genes diferencialmente expresados (DEGs)** entre Obeso2 y Normopeso:

### Sobre-expresados en Obeso2

| Gen | log2FC | padj |
|---|---|---|
| MC4R | +1.07 | 1.7×10⁻⁶ |
| FTO | +0.81 | 5.6×10⁻⁵ |
| SIM1 | +0.81 | 5.6×10⁻⁵ |
| UBR2 | +0.81 | 5.6×10⁻⁵ |

### Sub-expresados en Obeso2

| Gen | log2FC | padj |
|---|---|---|
| BDNF | -0.82 | 5.6×10⁻⁵ |
| LEP | -0.81 | 5.6×10⁻⁵ |
| NTRK2 | -0.81 | 5.6×10⁻⁵ |
| CADM2 | -0.81 | 5.6×10⁻⁵ |

El clustering jerárquico no supervisado separó correctamente las muestras por grupo, validando la robustez del perfil identificado.

---

## Software utilizado

- **R**: 4.4.1
- **Bioconductor**: 3.19
- **Rqc**: control de calidad
- **Rsubread**: alineamiento + featureCounts
- **DESeq2**: análisis de expresión diferencial
- **ggplot2 + ggrepel**: visualización
- **pheatmap**: heatmap
- **Biostrings**: manipulación de secuencias FASTA
- **openxlsx**: exportación a Excel

Sistema operativo: Windows 11  
Rtools 4.4 instalado para compilación de paquetes Bioconductor.

---

## Cómo reproducir el análisis

1. Clonar/copiar la carpeta del proyecto.
2. Abrir el archivo `.Rproj` en RStudio.
3. Ejecutar los scripts en orden numérico (00 → 09).
4. Las dependencias se instalan desde `scripts/00_install_packages.R`.

Tiempo total de ejecución estimado: < 5 minutos en hardware estándar (datos simulados ligeros, ~1 MB total).

---

## Notas y limitaciones

- **Datos simulados**: las cuentas y magnitudes de cambio son sintéticas. Los efectos detectados son moderados pero estadísticamente robustos.
- **Tamaño muestral reducido**: 3 réplicas biológicas por grupo. Adecuado para esta práctica académica pero insuficiente para estudios de descubrimiento real.
- **Modo single-end forzado**: la inconsistencia detectada en los archivos R2 del grupo Obeso2 impidió el análisis paired-end. El script `03_sincronizar_pairs.R` se generó como diagnóstico pero no se utilizó en el pipeline final.
- **UBR2**: gen con menor caracterización funcional en obesidad; sus hallazgos requieren validación experimental adicional.

---

## Referencias

Mahmoud, R., Kimonis, V., & Butler, M. G. (2022). Genetics of obesity in humans: A clinical review. *International Journal of Molecular Sciences*, 23(19), 11005.

Abuzzahab, M. J., et al. (2025). Improving the diagnosis of hyperphagia in melanocortin-4 receptor pathway diseases. *Obesity*, 33(7), 1217–1231.

Benak, D., et al. (2024). FTO in health and disease. *Frontiers in Cell and Developmental Biology*, 12, 1500394.

Barde, Y.-A. (2025). The physiopathology of brain-derived neurotrophic factor. *Physiological Reviews*, 105(4), 2073–2140.

Al Zein, M., et al. (2024). Leptin is a potential biomarker of childhood obesity.

Rong, S. S., & Yu, X. (2023). Phenotypic and genetic links between body fat measurements and primary open-angle glaucoma. *IJMS*, 24(4), 3925.

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts. *Bioinformatics*, 30(7), 923–930.

Liao, Y., Smyth, G. K., & Shi, W. (2019). The R package Rsubread. *Nucleic Acids Research*, 47(8), e47.
