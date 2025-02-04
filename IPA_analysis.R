install.packages("magick")
install.packages("pdftools")
# library(magick)

#IPA analysis

canonical <- vroom::vroom(file = "/datos/home/paulinapg/PhD_IPA_analysis/Results/CanonicalPathways-allphd_ROSMAP_DLFPC_DE_full_gene_list_dichoNIAReagan.csv")

dim(canonical)
canonical$pvalue <- 10^(-canonical$`-log(p-value)`)

canonical_sig <- canonical %>% filter(`-log(p-value)` > 1.3) %>% filter(`z-score` != "Err:520")
dim(canonical_sig)

# Ordenar los datos por -log(p-value) para que se visualicen de mayor a menor
canonical_sig <- canonical_sig %>%
  arrange(desc(`-log(p-value)`))

canonic.p <- ggplot(canonical_sig, aes(x = reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`))) +
  # Barras para -log(P-value) con leyenda
  geom_bar(aes(y = `-log(p-value)`, fill = "P-value"), stat = "identity", color = "black") +
  # Línea y puntos para Ratio
  geom_line(aes(y = Ratio * max(`-log(p-value)`) / max(Ratio), color = "Ratio"), group = 1, linetype = "dashed") +
  geom_point(aes(y = Ratio * max(`-log(p-value)`) / max(Ratio), color = "Ratio"), shape = 22, fill = "white", size = 3) +
  coord_flip() +
  # Escala dual de ejes con eje secundario para Ratio
  scale_y_continuous(
    name = "-Log (P-value)",
    sec.axis = sec_axis(~ . * max(canonical_sig$Ratio) / max(canonical_sig$`-log(p-value)`), name = "Ratio")
  ) +
  # Escalas manuales para crear la leyenda
  scale_fill_manual(name = "", values = c("P-value" = "gray")) +
  scale_color_manual(name = "", values = c("Ratio" = "black")) +
  # Etiquetas y tema
  labs(x = NULL, title = "a)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),  # Ajustar tamaño del texto del eje y
    legend.position = "top"                 # Posición de la leyenda
  ) +
  guides(
    fill = guide_legend(override.aes = list(color = "black")),  # Leyenda personalizada para fill
    color = guide_legend(override.aes = list(linetype = c("dashed"), shape = c(22)))  # Leyenda personalizada para color
  )

#Read images

canonical_bubbles <- image_read("/datos/home/paulinapg/PhD_IPA_analysis/Results/CanonicalPathwaysBubleChart-allphd_ROSMAP_DLFPC_DE_full_gene_list_dichoNIAReagan.pdf")
img2 <- image_read("/datos/home/paulinapg/PhD_IPA_analysis/Results/")
img3 <- image_read("/datos/home/paulinapg/PhD_IPA_analysis/Results/")

#Disease and function

dis_fun <- vroom::vroom(file = "/datos/home/paulinapg/PhD_IPA_analysis/Results/D")

