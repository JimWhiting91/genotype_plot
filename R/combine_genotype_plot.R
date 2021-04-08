#' Combine Genotype Plot Elements
#' 
#' Function combines all present elements of the genotype_plot() output
#'
#' @param genotype_plot_res List of Genotype Plot elements derived from running genotype_plot function.
#' @param heights Set the relative heights of the snp position labels and genotypes/dendrogram respectively
#' @param widths Set the relative widths of the dendrogram and genotyps/snp position labels respectively

#' @import cowplot
#' @import ggplot2
#' 
#' @return A fully composed genotype plot based on all available elements
#' @export
combine_genotype_plot<-function(genotype_plot_res,heights=c(1,9),widths=c(2,3)){
  
  # Check names match up
  if(class(genotype_plot_res) != "list"){
    stop("Error: Input is not of class list, unrecognised.")
  }
  if(paste(names(genotype_plot_res),collapse = "_") != paste(c("positions","genotypes","dendrogram","dendro_labels","cluster_pca"),collapse = "_")){
    stop("Error: names() of genotype_plot_res input do not match expected. Please provide output from genotype_plot().")
  }
  
  # Check for presence of a dendrogram
  if(is.null(genotype_plot_res$dendrogram)){
    
    # Only plot genotypes and positions
    plot_grid(genotype_plot_res$positions,
              genotype_plot_res$genotypes,
              ncol=1,nrow=2,axis="tblr",align="v",rel_heights = heights)
    
  } else {
    
    # Plot dendrogram, genotypes and positions
    empty_fig <- ggplot() + theme_void()
    plot_grid(plot_grid(empty_fig,genotype_plot_res$positions,rel_widths = widths,axis = "tblr",align="h",ncol=2,nrow=1),
              plot_grid(genotype_plot_res$dendrogram,genotype_plot_res$genotypes,rel_widths = widths,axis = "tblr",align="h",ncol=2,nrow=1),
              ncol=1,nrow=2,axis="tblr",align="v",rel_heights = heights)
    
  }
  
}