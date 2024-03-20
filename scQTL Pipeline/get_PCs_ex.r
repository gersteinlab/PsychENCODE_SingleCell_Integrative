
num_components = 100


paste("reached L2.3.IT")
dat <- read.table("/Users/dc/Desktop/rsch/pec2/v5_scqtls/cell_type_expr_matrices_syncronized/L2.3.IT.mrged.bed", sep="\t", header=TRUE, row.names=1)
dat_transposed = t(dat)
pca.out_trasp = prcomp(dat_transposed)
matrix_to_print = t(pca.out_trasp$x[,1:num_components])
output_file = "/Users/dc/Desktop/rsch/pec2/v5_scqtls/cell_type_expr_matrices_syncronized/output_PCA_files/L2.3.IT_PCs.dat"
write.table(matrix_to_print, output_file, append = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)


