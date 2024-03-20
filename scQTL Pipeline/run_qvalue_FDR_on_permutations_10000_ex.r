
library(qvalue)

fdr_thresh = 0.05

d = read.table("/Users/dc/Desktop/rsch/pec2/v5_scqtls/output_10000_permuts/permutations_Astro.100_expr_PCs__10000.txt", head=FALSE, stringsAsFactors=FALSE)
d=d[!is.na(d$V20),]
d$qval=qvalue(d$V20)$qvalue
write.table(d[which(d$qval <= fdr_thresh), ], "/Users/dc/Desktop/rsch/pec2/v5_scqtls/output_10000_permuts/permutations_Astro.100_expr_PCs__10000.sig_FDR_0.05.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


