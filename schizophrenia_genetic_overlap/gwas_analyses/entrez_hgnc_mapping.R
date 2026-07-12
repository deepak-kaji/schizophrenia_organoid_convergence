library(org.Hs.eg.db)
library(AnnotationDbi)

magma <- read.table(
  "~/programs/MAGMA/results/SCZ_PGC3_gene_results.genes.out",
  header=TRUE
)

magma$GENE <- as.character(magma$GENE)

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=magma$GENE,
  columns=c("SYMBOL"),
  keytype="ENTREZID"
)

magma_annot <- merge(
  magma,
  gene_map,
  by.x="GENE",
  by.y="ENTREZID",
  all.x=TRUE
)

write.csv(
  magma_annot,
  "~/programs/MAGMA/results/SCZ_PGC3_gene_results_annotated.csv",
  row.names=FALSE
)
