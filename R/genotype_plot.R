

genotype_plot<-function(vcf=NULL,
                        chr=NULL,
                        start=NULL,
                        end=NULL,
                        gff=NULL,
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
  # popmap = list of vectors corresponding to individuals in the VCF. Names of list should be pop names.
  # cluster = whether to organise haplotypes by hclust clustering
  # colour_scheme = character vector of colour values
  
  # Get what we need
  lib<-c("ggdendro","ggrepel","vcfR","ggplot2","reshape2","stringr","gggenes","cowplot","dplyr","grid","parallel","tidyr","data.table","ggpubr")
  lapply(lib,library,character.only=T)
  
  # Subset the VCF on the command-line for the chr of interest
  if(substr(vcf, nchar(vcf)-2+1, nchar(vcf)) == "gz"){
    system(paste0("bcftools view -r ",chr,":",start,"-",end," ",vcf," > gt_plot_tmp.vcf"),wait=TRUE)
  } else {
    stop("VCF needs tp bgzipped pal...")
    return()
  }
  
  # Read in the subsetted VCF
  vcf_in<-read.vcfR("gt_plot_tmp.vcf")
  
  # Filter the VCF for individuals...
  vcf_in <- vcf_in[,c("FORMAT",unlist(popmap))]
  
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
    Mb_vec<-as.integer(seq(0,max(SNP_pos$BP),by=snp_label_size))
    
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
  
  ### First, if we are clustering do that
  if(cluster==TRUE){
    # Pull genos again
    genos2 <- extract.gt(vcf_in,as.numeric = T)
    clust_genos <- hclust(dist(t(genos2)))
    
    # Plot dendrogram
    dendro <- ggdendrogram(clust_genos,rotate = T,labels = T)+
      theme_dendro()+
      scale_y_reverse()
    
    dendro_labels <- clust_genos$labels[clust_genos$order]
  } else {
    dendro <- NULL
    dendro_labels <- NULL
  }    
  
  ### open the region vcf
  vcf2<-t(extract.gt(vcf_in))
  colnames(vcf2)<-as.character(SNP_pos$GEN_pos)
  
  # Reorder
  if(cluster == FALSE){
  name_order<-data.frame(names=unlist(popmap),
                         index=length(unlist(popmap)):1)
  name_order<-name_order[order(name_order$names),]
  } else {
  name_order<-data.frame(names=dendro_labels,
                           index=1:length(unlist(popmap)))
  name_order<-name_order[order(name_order$names),]
  }
  
  vcf2<-vcf2[order(rownames(vcf2)),]
  rownames(vcf2)<-name_order$index
  genos<-reshape::melt(vcf2)
  colnames(genos)<-c('index','snp','GT')
  
  ### split the geno type by / into 2 columns
  genos<-genos %>%
    separate(GT, c("x1", "x2"), "/")
  genos$x1<-as.numeric(as.character(genos$x1))
  genos$x2<-as.numeric(as.character(genos$x2))
  
  ### sum the 2 columns
  # 0 = homoz for ref, 1 = heteroz, 2 = homoz for alt
  genos$GT<-genos$x1 + genos$x2
  
  # Pop cutoffs
  cutoffs<-data.frame(pops=names(popmap),
                      cutoffs=NA,
                      labs=NA)
  cutoffs[1,2]<-length(popmap[[1]])
  cutoffs[1,3]<-length(popmap[[1]])/2
  
  for(i in 2:nrow(cutoffs)){
    cutoffs[i,2]<-cutoffs[i-1,2]+length(popmap[[i]])
    cutoffs[i,3]<-mean(c(cutoffs[i-1,2],cutoffs[i,2]))
  }
  cutoffs$cutoffs2<-length(unlist(popmap))-cutoffs$cutoffs+0.5
  
  # Tidy up to get final label pos
  cutoffs$labs2<-NA
  cutoffs$labs2[1]<-mean(c(cutoffs$cutoffs2[1],length(unlist(popmap))))
  for(i in 2:nrow(cutoffs)){
    cutoffs[i,"labs2"]<-mean(c(cutoffs[i,"cutoffs2"],cutoffs[i-1,"cutoffs2"]))
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
    geom_hline(yintercept = c(length(unlist(popmap))+0.5,cutoffs$cutoffs2))+
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
  
  if(length(gff)==0){
   gene_plot <- NULL
  } else {
    
    #### Make genes plot for top
    # Read in the geneinfo and tidy
    genes<-data.frame(fread(gff))
    colnames(genes)[1]<-"chr"
    genes<-genes[genes$V3 == "gene" & genes$chr == chr,]
    genes<-genes[genes$V4 < end & genes$V5 > start,]
    genes$molecule<-chr
    colnames(genes)[4:5]<-c("startPoint","endPoint")
    genes$label_x<-rowMeans(genes[,c("endPoint","startPoint")])
    
    
    # get names
    genes<-genes %>%
      separate(V9, c("names1", "names2"), "source_gene_common_name=")
    genes<-genes %>%
      separate(names2, c("name", "bin"), ";")
    
    # Filter for uninformative
    genes[genes$name == "None","name"]<-NA
    
    # Plot
    gene_plot<-ggplot(genes, aes(xmin = startPoint, xmax = endPoint,y = molecule,label=name)) +
      geom_gene_arrow(arrowhead_height = unit(0,"mm"),arrowhead_width = unit(0,"mm"),arrow_body_height = unit(5,"mm")) +
      geom_text_repel(data=genes,aes(x=label_x),nudge_y      = 0.05,
                      direction    = "x",
                      angle        = 90,
                      vjust        = -0.5,
                      segment.size = 0.2,size=5) +
      #facet_wrap(~ molecule, scales = "free", ncol = 1) +
      scale_fill_brewer(palette = "Set3")+
      theme_nothing()+
      theme(panel.grid=element_blank(),
            #axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position="top",
            legend.text = element_text(size=18))+
      labs(fill="")+
      xlim(c(start,end))+
      coord_cartesian(ylim = c(12,0))+
      scale_x_continuous(expand = c(0, 0))
  }
    
    # Return everything for plots...
    output<-list(gene_plot,lines,genotypes,dendro,dendro_labels,gene_plot)
    names(output) <- c("genes","positions","genotypes","dendrogram","dendro_labels","gene_plot")
    return(output)

}

genotype_plot_pool<-function(vcf=NULL,chr=NULL,start=NULL,end=NULL,gff=NULL,popmap=NULL){
  
  # Usage:
  # vcf = a bgzipped VCF 
  # chr = chr of scaf ID
  # start = start of region
  # end = end of region
  # gff = file with gene annotations ".gff3" format
  # popmap = list of vectors corresponding to individuals in the VCF. Names of list should be pop names.
  
  # Get what we need
  lib<-c("ggrepel","vcfR","ggplot2","reshape2","stringr","gggenes","cowplot","dplyr","grid","parallel","tidyr","data.table","ggpubr","viridis")
  lapply(lib,library,character.only=T)
  
  # Subset the VCF on the command-line for the chr of interest
  if(substr(vcf, nchar(vcf)-2+1, nchar(vcf)) == "gz"){
    system(paste0("bcftools view -r ",chr,":",start,"-",end," ",vcf," > gt_plot_tmp.vcf"),wait=TRUE)
  } else {
    print("VCF needs to be bgzipped pal...")
    return()
  }
  
  # Read in the subsetted VCF
  vcf_in<-read.vcfR("gt_plot_tmp.vcf")
  
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
  if(max(SNP_pos$BP)-min(SNP_pos$BP) > 10000000){
    Mb_vec<-as.integer(seq(0,max(SNP_pos$BP),by=1000000))
  
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
  
  #### Because this is pooled, we'll plot depth...
  ### open the region vcf
  vcf2<-t(extract.gt(vcf_in,element = "AD"))
  colnames(vcf2)<-as.character(SNP_pos$GEN_pos)
  
  # Reorder
  name_order<-data.frame(names=unlist(popmap),
                         index=length(unlist(popmap)):1)
  name_order<-name_order[order(name_order$names),]
  vcf2<-vcf2[unlist(popmap),]
  vcf2<-vcf2[order(rownames(vcf2)),]
  rownames(vcf2)<-name_order$index
  genos<-melt(vcf2)
  colnames(genos)<-c('index','snp','GT')
  
  ### split the geno type by / into 2 columns
  genos<-genos %>%
    separate(GT, c("x1", "x2"), ",")
  genos$x1<-as.numeric(as.character(genos$x1))
  genos$x2<-as.numeric(as.character(genos$x2))
  
  ### sum the 2 columns
  # 0 = homoz for ref, 1 = heteroz, 2 = homoz for alt
  genos$GT<-genos$x1/(genos$x1 + genos$x2)
  
  # Pop cutoffs
  cutoffs<-data.frame(pops=names(popmap),
                      cutoffs=NA,
                      labs=NA)
  cutoffs[1,2]<-length(popmap[[1]])
  cutoffs[1,3]<-length(popmap[[1]])/2
  
  for(i in 2:nrow(cutoffs)){
    cutoffs[i,2]<-cutoffs[i-1,2]+length(popmap[[i]])
    cutoffs[i,3]<-mean(c(cutoffs[i-1,2],cutoffs[i,2]))
  }
  cutoffs$cutoffs2<-length(unlist(popmap))-cutoffs$cutoffs+0.5
  
  # Tidy up to get final label pos
  cutoffs$labs2<-NA
  cutoffs$labs2[1]<-length(unlist(popmap))
  for(i in 2:nrow(cutoffs)){
    cutoffs[i,"labs2"]<-mean(c(cutoffs[i,"cutoffs2"],cutoffs[i-1,"cutoffs2"]))
  }
  
  # Plot genotypes
  
  genotypes<-ggplot(data=genos,aes(x=snp,y=index))+
    #geom_raster(aes(fill=factor(GT)))+
    geom_tile(aes(fill=GT))+
    #scale_fill_manual(values=c("#d4b9da","#e7298a","#980043"),name="Genotype",
    #                  breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT"))+
    scale_fill_viridis(option = "inferno")+
    theme_bw()+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.ticks  = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=15),
          legend.position = "bottom",
          legend.title=element_text(size=14),
          legend.text = element_text(size=12,angle=45,hjust=1),
          panel.border = element_blank())+
    geom_hline(yintercept = c(length(unlist(popmap))+0.5,cutoffs$cutoffs2))+
    scale_y_continuous(breaks = as.integer(cutoffs$labs2),
                       labels = cutoffs$pops)+
    scale_x_continuous(expand = c(0, 0))+
    labs(fill = "REF/ALT")

  
  if(length(gff)==0){
    gene_plot <- NULL
  } else {
    
    #### Make genes plot for top
    # Read in the geneinfo and tidy
    genes<-data.frame(fread(gff))
    colnames(genes)[1]<-"chr"
    genes<-genes[genes$V3 == "gene" & genes$chr == chr,]
    genes<-genes[genes$V4 < end & genes$V5 > start,]
    genes$molecule<-chr
    colnames(genes)[4:5]<-c("startPoint","endPoint")
    genes$label_x<-rowMeans(genes[,c("endPoint","startPoint")])
    
    
    # get names
    genes<-genes %>%
      separate(V9, c("names1", "names2"), "source_gene_common_name=")
    genes<-genes %>%
      separate(names2, c("name", "bin"), ";")
    
    # Filter for uninformative
    genes[genes$name == "None","name"]<-NA
    
    # Plot
    gene_plot<-ggplot(genes, aes(xmin = startPoint, xmax = endPoint,y = molecule,label=name)) +
      geom_gene_arrow(arrowhead_height = unit(0,"mm"),arrowhead_width = unit(0,"mm"),arrow_body_height = unit(5,"mm")) +
      geom_text_repel(data=genes,aes(x=label_x),nudge_y      = 0.05,
                      direction    = "x",
                      angle        = 90,
                      vjust        = -0.5,
                      segment.size = 0.2,size=5) +
      #facet_wrap(~ molecule, scales = "free", ncol = 1) +
      scale_fill_brewer(palette = "Set3")+
      theme_nothing()+
      theme(panel.grid=element_blank(),
            #axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position="top",
            legend.text = element_text(size=18))+
      labs(fill="")+
      xlim(c(start,end))+
      coord_cartesian(ylim = c(12,0))+
      scale_x_continuous(expand = c(0, 0))
    
    return(list(gene_plot,lines,genotypes))
  }
}

# # Test pool
# vcf.t = "~/Exeter/tmp/poolseq_freebayes.final.vcf.gz" 
# chr.t = "chr12"
# start.t = 1
# end.t = 40000000
# #gff = file with gene annotations ".gff3" format
# popmap.t = list(c("IF10"),c("IF6"),c("IF8"),c("IF9"))
# names(popmap.t)<-c("IF10","IF6","IF8","IF9")
# 
# p1<-genotype_plot_pool(vcf=vcf.t,chr=chr.t,start=start.t,end=end.t,popmap=popmap.t)
# 
# # TEST pool
# pdf("~/Desktop/pool_genotype_test_chr12.pdf",width=16,height=6)
# plot_grid(p1[[1]],
#           p1[[2]],
#            rel_heights = c(1,4),ncol=1,nrow=2,align = "v",axis = "lr")
# dev.off()

# 
# # Set up popmap and name populations
# pops <- c("TACHP","TACLP","GHP","GLP","APHP","APLP","OHP","OLP","MADHP","MADLP")
# pops2 <- c("LT","UT","GH","GL","APHP","APLP","LO","UQ","LMD","UMD")
# 
# popmap <- lapply(pops2,function(x){
#   as.character(read.table(paste0("~/Exeter/VCFs/",x,".popmap"))[,1])
# })
# names(popmap)<-pops
# tester<-genotype_plot(vcf="~/Exeter/VCFs/five_aside_STAR_3033083_final.vcf.gz",
#                       chr="chr20",
#                       start=1,
#                       end=800000,
#                       popmap = popmap,
#                       cluster=TRUE)

################
# # Josie region 1test
# paria<-c("Paria_F12",
#          "Paria_F16",
#          "Paria_F2",
#          "Paria_F2_RNA",
#          "Paria_F8",
#          "Paria_M12",
#          "Paria_M14",
#          "Paria_M2",
#          "Paria_M5",
#          "Paria_M9")
# UM<-c("UM1",
#       "UM10",
#       "UM11",
#       "UM12",
#       "UM13",
#       "UM14",
#       "UM15",
#       "UM2",
#       "UM3",
#       "UM5",
#       "UM6",
#       "UM7")
# 
# popmap <- list(paria,UM)
# names(popmap) <- c("Paria","UM")
# 
# josie_region1 <- genotype_plot(vcf="~/Exeter/VCFs/chr5.test.recode.vcf.gz",
#                                chr="chr5",
#                                start=9665770,
#                                end=17067807,
#                                popmap = popmap,
#                                cluster = FALSE)
# 
# plot_grid(josie_region1[[1]],
#           josie_region1[[2]],
#           rel_heights = c(1,5),ncol=1,
#           axis="tblr",align="v")


