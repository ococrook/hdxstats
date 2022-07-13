library("RColorBrewer")
library("scales")

averaged_protection <- read.csv("vignettes/data/averaged_protection_both.csv")
antibody_selection <- colnames(averaged_protection)[c(2:length(averaged_protection))]

heatmap_data <- list()

for (antibody in antibody_selection) {
  filename_csv <- paste("vignettes/data/hdx_", antibody, "_labelled_residues_pdb_5edv_chainA.csv", sep="")
  df <- data.frame(read.csv(filename_csv))
  
  hdx_values <- c()
  for (i in seq(length(df$aa))){
    if (df$labels[i] == "deprotected") {
      hdx_values <- append(hdx_values, -1*df$hdx_values[i])
    } else{
      hdx_values <- append(hdx_values, df$hdx_values[i])
    }
  }
  heatmap_data[[antibody]] <- hdx_values
}

heatmap_matrix <- data.matrix(heatmap_data)
heatmap_matrix <- matrix(unlist(heatmap_data), ncol = 245, nrow = 18)
rownames(heatmap_matrix) <- antibody_selection
colnames(heatmap_matrix) <- df$aa

colormap <- colorRamp(brewer.pal(8, "RdBu"))
domain <- c(-1,1) # min and max values
color_function <- col_bin(colormap, domain)

#plot_ly(z=heatmap_matrix, color="RdBu")
# heatmap(heatmap_matrix, Colv = NA, Rowv = NA, scale="column", cexRow=1.5, , cexCol=1.5, xlab="HOIP-RBR sequence", ylab="Antibody", col=rainbow(25))
#heatmap(heatmap_matrix, Colv = NA, Rowv = NA, scale="column", xlab="HOIP-RBR sequence", ylab="Antibody", col=rainbow(25), margins=c(5,5))
