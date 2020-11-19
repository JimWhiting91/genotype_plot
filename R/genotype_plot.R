genotype_plot<-function(vcf=NULL,
                        chr=NULL,
                        start=0,
                        end=NULL,
                        popmap=NULL,
                        cluster=FALSE,
                        snp_label_size=500000,
                        colour_scheme=c("#d4b9da","#e7298a","#980043")){
  
  # Usage:
  # vcf = a bgzipped VCF 
  # chr = chr of scaf ID
  # start = start of region
  # end = end of region
  # gff = file with gene annotations ".gff3" format
  # popmap = two column data frame with column 1 for individual IDs as they appear in the VCF and column 2 for pop labels
  # cluster = whether to organise haplotypes by hclust clustering
  # colour_scheme = character vector of colour values
  
  # Note, the gff functionality is currently not used.
  
  # Get what we need
  lib<-c("ggdendro","vcfR","ggplot2","tidyr","data.table")
  lapply(lib,library,character.only=T)
  
  # Subset the VCF on the command-line for the chr of interest
  if(substr(vcf, nchar(vcf)-2+1, nchar(vcf)) == "gz"){
    system(paste0("bcftools view -r ",chr,":",start,"-",end," ",vcf," > gt_plot_tmp.vcf"),wait=TRUE)
  } else {
    stop("VCF needs to be bgzipped")
  }
  
  # Tidy up the popmap and turn into a list
  popmap2 <- lapply(unique(popmap[,2]),function(pop){
    return(as.character(popmap[popmap[,2]==pop,1]))
  })
  names(popmap2)<-unique(popmap[,2])
  
  # Read in the subsetted VCF
  vcf_in<-read.vcfR("gt_plot_tmp.vcf")
  
  # Filter the VCF for individuals...
  vcf_in <- vcf_in[,c("FORMAT",popmap[,1])]
  
  # Remove invariants following filtering
  vcf_in <- vcf_in[is.polymorphic(vcf_in,na.omit=T),]
  
  # Keep tidy
  system("rm -f gt_plot_tmp.vcf")
  
  
  ### Make the chromosome connecting lines
  SNP_pos<-data.frame(chr=chr,
                      pos=as.integer(vcf_in@fix[,2]))
  SNP_pos$y1<--0.5
  SNP_pos$y2<--1.5
  colnames(SNP_pos)<-c("chrom","BP","y1","y2")
  SNP_pos$index3<-seq(min(SNP_pos$BP),max(SNP_pos$BP),
                      by=(max(SNP_pos$BP)-min(SNP_pos$BP))/nrow(SNP_pos))[1:nrow(SNP_pos)]
  
  colnames(SNP_pos)<-c("chrom","BP","y1","y2","GEN_pos")
  
  # We also need a filter for lines if we're plotting whole chr...
  if(max(SNP_pos$BP)-min(SNP_pos$BP) > snp_label_size){
    Mb_vec<-as.integer(seq(start,end,by=snp_label_size))
    
    filtered_SNPs<-data.frame(rbindlist(lapply(Mb_vec,function(x){
      return(SNP_pos[abs(SNP_pos$BP-x) == min(abs(SNP_pos$BP-x)),])
    })))
    
    lines<-ggplot(filtered_SNPs,aes(x=as.factor(BP),y=y))+
      geom_segment(aes(x=BP,xend=GEN_pos,y=y1,yend=y2))+
      theme_classic()+
      theme(axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
      xlim(c(start,end))+
      scale_x_continuous(expand = c(0, 0),position="top",
                         breaks = Mb_vec,
                         labels = Mb_vec/1000000)
    
  } else {
    
    lines<-ggplot(SNP_pos,aes(x=as.factor(BP),y=y))+
      geom_segment(aes(x=BP,xend=GEN_pos,y=y1,yend=y2))+
      theme_classic()+
      theme(axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
      xlim(c(start,end))+
      scale_x_continuous(expand = c(0, 0),position="top")
    
  }
  
  # Now get our genotypes and transform
  vcf2<-t(extract.gt(vcf_in))
  colnames(vcf2)<-as.character(SNP_pos$GEN_pos)
  # Pull genos again
  genos<-reshape::melt(vcf2)
  colnames(genos)<-c('index','snp','GT')
  
  ### split the geno type by / into 2 columns
  genos<-genos %>%
    separate(GT, c("x1", "x2"), "[|/]")
  genos$x1<-as.numeric(as.character(genos$x1))
  genos$x2<-as.numeric(as.character(genos$x2))
  genos$GT <- genos$x1 + genos$x2
  
  ### First, if we are clustering do that
  if(cluster){
    
    # Remake numeric genotype matrix
    geno_matrix <- as.matrix(pivot_wider(genos[,c("index","snp","GT")], names_from = snp, values_from = GT))
    matrix_inds <- geno_matrix[,1]
    geno_matrix <- geno_matrix[,2:ncol(geno_matrix)]
    geno_matrix <- gsub(" ","",geno_matrix)
    colnames(geno_matrix) <- SNP_pos$BP
    rownames(geno_matrix) <- matrix_inds
    class(geno_matrix) <- "numeric"
    
    # Cluster with hclust
    clust_genos <- hclust(dist(geno_matrix))
    
    # Plot dendrogram
    dendro <- ggdendrogram(clust_genos,rotate = T,labels = T)+
      theme_dendro()+
      scale_y_reverse()
    
    dendro_labels <- clust_genos$labels[clust_genos$order]
  } else {
    dendro <- NULL
    dendro_labels <- NULL
  }    
    
  # Reorder
  if(cluster == FALSE){
    name_order<-data.frame(names=unlist(popmap2),
                           index=length(unlist(popmap2)):1)
    name_order<-name_order[order(name_order$names),]
  } else {
    name_order<-data.frame(names=dendro_labels,
                           index=1:length(unlist(popmap2)))
    name_order<-name_order[order(name_order$names),]
  }
  
  # Take our genotypes again but now reorder based on clustering/no clustering index
  vcf2<-vcf2[order(rownames(vcf2)),]
  rownames(vcf2)<-name_order$index
  genos<-reshape::melt(vcf2)
  colnames(genos)<-c('index','snp','GT')

  ### split the geno type by / into 2 columns
  genos<-genos %>%
    separate(GT, c("x1", "x2"), "[|/]")
  genos$x1<-as.numeric(as.character(genos$x1))
  genos$x2<-as.numeric(as.character(genos$x2))

  ### sum the 2 columns
  # 0 = homoz for ref, 1 = heteroz, 2 = homoz for alt
  genos$GT<-genos$x1 + genos$x2

  # Pop cutoffs
  cutoffs<-data.frame(pops=names(popmap2),
                      cutoffs=NA,
                      labs=NA)
  cutoffs[1,2]<-length(popmap2[[1]])
  cutoffs[1,3]<-length(popmap2[[1]])/2
  
  if(nrow(cutoffs) > 1){
    for(i in 2:nrow(cutoffs)){
      cutoffs[i,2]<-cutoffs[i-1,2]+length(popmap2[[i]])
      cutoffs[i,3]<-mean(c(cutoffs[i-1,2],cutoffs[i,2]))
    }
  }
  cutoffs$cutoffs2<-length(unlist(popmap2))-cutoffs$cutoffs+0.5
  
  # Tidy up to get final label pos
  cutoffs$labs2<-NA
  cutoffs$labs2[1]<-mean(c(cutoffs$cutoffs2[1],length(unlist(popmap2))))
  if(nrow(cutoffs) > 1){
    for(i in 2:nrow(cutoffs)){
      cutoffs[i,"labs2"]<-mean(c(cutoffs[i,"cutoffs2"],cutoffs[i-1,"cutoffs2"]))
    }
  }
  
  # Catch bad labelling for individuals
  if(length(unique(popmap[,2])) == nrow(popmap)){
    cutoffs[1,"labs2"] <- nrow(popmap)
  }
  
  # Plot genotypes based on clustering...
  if(cluster == FALSE){
    genotypes<-ggplot(data=genos,aes(x=snp,y=index))+
      #geom_raster(aes(fill=factor(GT)))+
      geom_tile(aes(fill=factor(GT)))+
      scale_fill_manual(values=colour_scheme,name="Genotype",
                        breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))+
      theme_bw()+
      theme(panel.grid = element_blank(), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.ticks  = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=15),
            legend.position = "bottom",
            legend.title=element_text(size=14),
            legend.text = element_text(size=14),
            panel.border = element_blank())+
      geom_hline(yintercept = c(length(unlist(popmap2))+0.5,cutoffs$cutoffs2))+
      scale_y_continuous(breaks = as.integer(cutoffs$labs2),
                         labels = cutoffs$pops)+
      scale_x_continuous(expand = c(0, 0))
    
  } else {
    genotypes<-ggplot(data=genos,aes(x=snp,y=index))+
      #geom_raster(aes(fill=factor(GT)))+
      geom_tile(aes(fill=factor(GT)))+
      scale_fill_manual(values=colour_scheme,name="Genotype",
                        breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))+
      theme_bw()+
      theme(panel.grid = element_blank(), 
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), 
            axis.ticks  = element_blank(),
            axis.text = element_blank(),
            legend.position = "bottom",
            legend.title=element_text(size=14),
            legend.text = element_text(size=14),
            panel.border = element_blank())+
      scale_x_continuous(expand = c(0, 0))
  }
  
  # Return everything for plots...
  output<-list(lines,genotypes,dendro,dendro_labels)
  names(output) <- c("positions","genotypes","dendrogram","dendro_labels")
  return(output)
  
}
