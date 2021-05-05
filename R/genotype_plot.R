#' Genotype Plotting
#' 
#' Visualises genotypes within a VCF file located somewhere on the system, or loaded as a vcfR object
#'
#' @param vcf Path to a bgzipped vcf as a character string.
#' @param vcf_object A vcfR object to plot, specified as an alternative to vcf, chr, start and end.
#' @param chr Chr of scaf ID as a character string.
#' @param start Start of region as numeric/integer.
#' @param end End of region as numeric/integer.
#' @param popmap Two column data frame with column 1 for individual IDs as they appear in the VCF and column 2 for pop labels.
#' @param cluster Logical TRUE/FALSE whether to organise haplotypes by hclust clustering.
#' @param snp_label_size Physical distance between snp label markers as numeric/integer.
#' @param colour_scheme Character vector of colour values, must be length 3.
#' @param invariant_filter Logical TRUE/FALSE to remove invariants.
#' @param is.multi_allelic Logical TRUE/FALSE to handle multi-allelic VCFs.
#' @param is.haploid Logical TRUE/FALSE to handle haploid genotypes.
#' @param polarise_genotypes Character string that matches a population in the popmap, genotypes are polarised to the major allele in this population.
#' @param plot_allele_frequency Logical TRUE/FALSE to switch on/off plotting per-population allele frequencies.
#' @param plot_phased Logical TRUE/FALSE to determine whether we plot genotypes or phased haplotypes.
#' @param missingness Numeric between 0-1 to filter for missing genotypes.
#'
#' @import data.table
#' @import ggdendro
#' @import ggplot2
#' @importFrom glue glue
#' @import tidyr
#' @import vcfR
#' @import adegenet
#' @import pegas
#' @import viridis
#' @importFrom ade4 dudi.pca
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats na.omit
#' @importFrom reshape2 melt
#' @importFrom bedr check.binary
#' 
#' @return A list of plottable objects including genotypes, SNP labels, dendrogram
#' @export
genotype_plot<-function(vcf=NULL,
                        chr=NULL,
                        start=0,
                        end=NULL,
                        popmap=NULL,
                        cluster=FALSE,
                        snp_label_size=100000,
                        colour_scheme=c("#FCD225","#C92D59","#300060"),
                        vcf_object=NULL,
                        invariant_filter=TRUE,
                        is.multi_allelic=FALSE,
                        is.haploid=FALSE,
                        polarise_genotypes=NULL,
                        plot_allele_frequency=FALSE,
                        plot_phased=FALSE,
                        missingness=0.5){
  
  # Bind variables to function
  AF <- BP <- GEN_pos <- GT <- index <- pop_F <- snp <- snp_pos <- y <- y1 <- y2 <- NULL
  
  # Catch basic errors
  if(plot_allele_frequency & plot_phased){
    stop("Error: Cannot plot both allele frequencies and phased haplotypes")
  }
  if(plot_allele_frequency & cluster){
    stop("Error: Cannot plot allele frequencies and cluster")
  }
  
  # Reverse order of popmap if we are plotting AF to catch bug
  # Also do a check that we have >1 ind per pop
  if(plot_allele_frequency){
    popmap <- popmap[order(match(popmap[,2],rev(unique(popmap[,2])))),]
    if(any(table(popmap[,2]) < 2)){
      stop("Error: Cannot plot AFs with only one individual per pop. Check popmap.")
    }
  }
  
  # Read in our VCF if we need to
  if(is.null(vcf_object)){
    
    # First check that bcftools is accesible
    if(!(check.binary("bcftools"))){
      stop("Error: BCFtools is not visible to R's system() PATH, check installation or use vcf_object input")
    } else {
      
      # Subset the VCF on the command-line for the chr of interest
      vcf_tempfile <- tempfile(pattern = "gt_plot", fileext = '.vcf')
      
      on.exit({ unlink(vcf_tempfile) })
      if(tools::file_ext(vcf) == "gz"){
        # Make the popmap temp
        popmap_tmp <- system(glue::glue("bcftools query -l {vcf}"), wait=TRUE, intern=T)
        # Validate user popmap
        if(any(!(popmap[,1] %in% popmap_tmp))){
          missing <- paste(popmap[,1][!(popmap[,1] %in% popmap_tmp)],collapse = ", ")
          stop(paste0("ERROR The following inds are not in VCF: ",missing))
        } else {
          # Read in individuals in popmap order...
          inds_to_read <- paste(popmap[,1],collapse=",")
          system(glue::glue("bcftools view -r {chr}:{start}-{end} -s {inds_to_read} {vcf} > {vcf_tempfile}"), wait=TRUE)
          
          # If needs be, also co-erce multi-allelic SNPs to biallelic...
          if(is.multi_allelic){
            message("Coercing Multi-allelic VCF to Bi-allelic VCF with bcftools norm -m -")
            vcf_tempfile2 <- tempfile(pattern = "gt_plot_bi", fileext = '.vcf')
            on.exit({ unlink(vcf_tempfile2) })
            system(glue::glue("bcftools norm -m - {vcf_tempfile} > {vcf_tempfile2}"))
            vcf_in<-read.vcfR(vcf_tempfile2)
          } else {
            # Read in the subsetted VCF
            vcf_in<-read.vcfR(vcf_tempfile)
          }
        }
      } else {
        stop("VCF needs to be bgzipped, or path mis-specified")
      }
    }
  } else {
    # Allow VCF object to be given instead of outer subsetting
    vcf_in <- vcf_object
    if(length(unique(vcf_in@fix[,1])) != 1){
      stop(paste0("ERROR multiple chr/scaf detected, VCF should contain a single chr/scaf for plotting"))
    }
    chr <- unique(vcf_in@fix[,1])
    start <- min(as.integer(vcf_in@fix[,2]))
    end <- max(as.integer(vcf_in@fix[,2]))
  }
  
  # If we have included multi-allelic sites, we will need to rename the IDs
  if(is.multi_allelic){
    
    # Fetch the non-unique names
    duplicated_id <- vcf_in@fix[,3][duplicated(vcf_in@fix[,3])]
    duplicated_counts <- table(vcf_in@fix[,3])
    duplicated_counts <- duplicated_counts[duplicated_id]
    
    # Loop over and replace
    for(id in duplicated_id){
      new_ids <- paste0(id,"_",letters[1:duplicated_counts[id]])
      vcf_in@fix[vcf_in@fix[,3] == id,3] <- new_ids 
    }
  }
  
  # Filter for missingness
  per_site_missing <- apply(extract.gt(vcf_in), MARGIN = 1, function(x){ sum(is.na(x)) })
  per_site_missing <- per_site_missing/(ncol(vcf_in@gt)-1)
  keep_missing <- which(per_site_missing < missingness)
  message(paste0("Removing ",length(per_site_missing)-length(keep_missing)," SNPs with > ",missingness * 100,"% missing data"))
  vcf_in <- vcf_in[keep_missing,]
  
  # Tidy up the popmap and turn into a list
  popmap2 <- lapply(unique(popmap[,2]),function(pop){
    
    # Adjust the popmap in anticipation of phasing...
    if(plot_phased){
      popmap_phased <- data.frame(rbindlist(lapply(popmap[,1],function(ind){
        popmap_tmp <- data.frame(ind=paste0(ind,"_",c(0,1)),
                                 pop=popmap[popmap[,1]==ind,2])
      })))
      return(as.character(popmap_phased[popmap_phased[,2]==pop,1]))
    } else {
      return(as.character(popmap[popmap[,2]==pop,1]))
    }
  })
  names(popmap2)<-unique(popmap[,2])
  
  # Remove invariants following filtering
  if(invariant_filter){
    to_prune <- is.polymorphic(vcf_in,na.omit = T)
    if(!(table(to_prune)["TRUE"] == nrow(vcf_in@fix))){
      message(paste0(table(to_prune)["FALSE"]," invariants have been pruned"))
      vcf_in <- vcf_in[is.polymorphic(vcf_in, na.omit=TRUE),]
    }
  }
  
  ### Make the chromosome connecting lines
  SNP_pos<-data.frame(chr=chr,
                      pos=as.integer(vcf_in@fix[,2]))
  SNP_pos$y1 <- -0.5
  SNP_pos$y2 <- -1.5
  colnames(SNP_pos) <- c("chrom","BP","y1","y2")
  
  SNP_pos$index3<-seq(min(SNP_pos$BP),max(SNP_pos$BP),
                      by=(max(SNP_pos$BP)-min(SNP_pos$BP))/nrow(SNP_pos))[1:nrow(SNP_pos)]
  
  colnames(SNP_pos)<-c("chrom","BP","y1","y2","GEN_pos")
  
  # We also need a filter for lines if we're plotting whole chr...
  message("Plotting SNP label markers")
  if(max(SNP_pos$BP)-min(SNP_pos$BP) > snp_label_size){
    
    Mb_vec<-as.integer(seq(0,end+snp_label_size,by=snp_label_size))
    Mb_vec <- c(start,Mb_vec[Mb_vec > start])
    Mb_vec[length(Mb_vec)] <- end
    
    filtered_SNPs<-data.frame(rbindlist(lapply(Mb_vec,function(x){
      return(SNP_pos[abs(SNP_pos$BP-x) == min(abs(SNP_pos$BP-x)),])
    })))
    
    lines <- suppressMessages(
      ggplot(filtered_SNPs,aes(x=as.factor(BP),y=y))+
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
                           labels = round(Mb_vec/1000000,3))
    )
    
  } else {
    
    lines <- suppressMessages(
      ggplot(SNP_pos,aes(x=as.factor(BP),y=y))+
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
    )
  }
  
  # Now get our genotypes and transform
  vcf2<-t(extract.gt(vcf_in))
  
  # Convert to hom genotypes if the vcf is haploid
  if(is.haploid){
    if(any(names(table(is.biallelic(vcf_in))) == "FALSE")){
      stop("Error, haploid VCF contains multiallelic loci. Please remove or coerce to biallelic sites before plotting.")
    }
    vcf2[vcf2=="0"] <- "0/0"
    vcf2[vcf2=="1"] <- "1/1"
  }
  
  # If we are polarising the VCF we will do that here
  if(!(is.null(polarise_genotypes))){
    # Don't handle mutli-allelic
    # if(is.multi_allelic){
    #   stop("Currently allele frequency options are not supported for multiallelic VCFs")
    # }
    # Calculate the AF for the focal population
    message(paste0("Polarising genotypes to major alleles of ",polarise_genotypes))
    af <- by(vcfR::vcfR2loci(vcf_in[,c(1,which(colnames(vcf_in@gt) %in% popmap[popmap[,2] == polarise_genotypes,1]))]),
             popmap[popmap[,2] == polarise_genotypes,2])
    
    # Return the major allele based on counts
    major_alleles <- sapply(af,function(i){
      # Mask bum windows
      if(ncol(i) == 0){
        return(0)
      } else {
        if(length(colnames(i)) != 2){
          # What col do we want?
          new_col <- ifelse(colnames(i)=="0","1","0")
          i <- cbind(i,0)
          colnames(i) <- c(colnames(i)[1],new_col)
        }
        return(ifelse(i[,"0"] > i[,"1"],0,1))
      }
    })  
    # End polarising
  }
  
  # Reorder if needs be
  if(!(is.null(vcf_object))){
    if(any(!(popmap[,1] %in% rownames(vcf2)))){
      missing <- paste(popmap[,1][!(popmap[,1] %in% rownames(vcf2))],collapse = ", ")
      stop(paste0("ERROR The following inds are not in vcfR object: ",missing))
    } else {
      vcf2 <- vcf2[popmap[,1],]
    }
  }
  
  colnames(vcf2)<-as.character(SNP_pos$GEN_pos)
  
  # Now polarise if we need to...
  if(!(is.null(polarise_genotypes))){
    to_polarise <- which(major_alleles == 1)
    
    # Get the dividier
    geno_format <- strsplit(unique(na.omit(vcf2[1,])),"")[[1]][2]
    # Check this
    if(geno_format != "|" & plot_phased){
      stop("Error: VCF is not in phased format (e.g. 0|1)")
    }
    
    vcf2[,to_polarise][vcf2[,to_polarise] == paste0("0",geno_format,"0")]<-"holder"
    vcf2[,to_polarise][vcf2[,to_polarise] == paste0("1",geno_format,"1")]<-paste0("0",geno_format,"0")
    vcf2[,to_polarise][vcf2[,to_polarise] == "holder"]<-paste0("1",geno_format,"1")
  }
  
  # Calculate allele frequencies for populations here if needs be
  if(plot_allele_frequency){
    message("Calculating per population allele frequencies")
    # First calculate all population allele frequencies
    # Extract gt manually from VCF
    gt <- t(vcf2)
    
    # Get the dividier
    geno_format <- strsplit(unique(na.omit(vcf2[1,]))[1],"")[[1]][2]
    
    # Convert to numbers
    gt[gt == paste0("0",geno_format,"0")] <- 0
    gt[gt == paste0("0",geno_format,"1")] <- 0.5
    gt[gt == paste0("1",geno_format,"0")] <- 0.5
    gt[gt == paste0("1",geno_format,"1")] <- 1
    class(gt) <- "numeric"
    
    # Filter the genotypes based on popmap
    gt <- gt[,colnames(gt) %in% popmap[,1]]
    
    # For all populations, calculate the allele frequency of the ALT allele (it doesn't matter which we use)
    pops <- unique(unlist(popmap[,2]))
    
    pop_AF <- lapply(pops,function(pop){
      
      # Get only these pops...
      gt_tmp <- gt[,popmap[popmap[,2]==pop,1]]
      
      # Get the AFs
      allele_counts <- rowSums(gt_tmp,na.rm = T)
      
      # Get counts, and factor in missingness...
      total_counts <- ncol(gt_tmp) - apply(gt_tmp, 1, function(x) sum(is.na(x)))
      
      return(allele_counts/total_counts)
      
    })
    
    # Make a new matrix
    AF_mat <- matrix(ncol=length(pops),nrow=nrow(gt))
    for(i in 1:ncol(AF_mat)){
      AF_mat[,i] <- pop_AF[[i]]
    }
    colnames(AF_mat) <- pops
    
    # Get REF/polarised if needed
    AF_mat <- 1 - AF_mat
    
    # Clean
    rownames(AF_mat) <- colnames(vcf2)
    AF_mat_long <- reshape2::melt(AF_mat)
    colnames(AF_mat_long) <- c("snp_pos","pop","AF")
    
  }
  
  # Pull genos again
  if(!plot_allele_frequency){
    message("Converting genotypes for plotting...")
    genos<-reshape2::melt(vcf2)
    colnames(genos)<-c('index','snp','GT')
    genos<-genos %>%
      separate(GT, c("x1", "x2"), "[|/]")
    
    ### split the geno type by / into 2 columns
    if(plot_phased){
      message("Plotting phased haplotypes instead of genotypes")
      
      # Check this
      # Get the divider
      geno_format <- strsplit(unique(na.omit(vcf2[1,])),"")[[1]][2]
      if(geno_format != "|" & plot_phased){
        stop("Error: VCF is not in phased format (e.g. 0|1)")
      }
      
      # rbind each genotype to make haplotypes
      genos <- data.frame(rbindlist(lapply(c(1,2),function(hap){
        genos_tmp <- genos[,c("index","snp",paste0("x",hap))]
        colnames(genos_tmp)[3] <- "GT"
        genos_tmp$index <- paste0(genos_tmp$index,"_",hap-1)
        genos_tmp
      })))
      
    } else {
      
      genos$x1<-as.numeric(as.character(genos$x1))
      genos$x2<-as.numeric(as.character(genos$x2))
      genos$GT <- genos$x1 + genos$x2
    }
  }
  
  ### First, if we are clustering do that
  if(cluster){
    if(plot_allele_frequency){
      stop("Error: Can't cluster haplotypes based on allele frequencies")
    } else {
      
      # Remake numeric genotype matrix
      geno_matrix <- as.matrix(pivot_wider(genos[,c("index","snp","GT")], names_from = snp, values_from = GT))
      matrix_inds <- geno_matrix[,1]
      geno_matrix <- geno_matrix[,2:ncol(geno_matrix)]
      geno_matrix <- gsub(" ","",geno_matrix)
      colnames(geno_matrix) <- SNP_pos$BP
      rownames(geno_matrix) <- matrix_inds
      class(geno_matrix) <- "numeric"
      
      #-------------- Cluster based on PC space
      # If data is phased, we treat haplotypes as haploid individuals
      if(!plot_phased){
        pca_gen <- df2genind(geno_matrix,ncode = 1,ploidy = 1)
        x.gen<-tab(pca_gen,freq=TRUE,NA.method="mean")
      } else {
        geno_matrix2 <- geno_matrix
        geno_matrix2[geno_matrix2 == 0] <- "0/0"
        geno_matrix2[geno_matrix2 == 1] <- "0/1"
        geno_matrix2[geno_matrix2 == 2] <- "1/1"
        pca_gen <- df2genind(geno_matrix2,ncode = 1,ploidy = 2,sep="/")
        x.gen<-tab(pca_gen,freq=TRUE,NA.method="mean")
      }
      
      # Do PCA
      cluster_pca <- dudi.pca(x.gen,center=T,scale = F,scannf = FALSE, nf=10)
      
      # # Cluster with hclust
      dist_genos <- dist(cluster_pca$li)
      
      # Cluster filtered
      clust_genos <- hclust(dist_genos)
      # Plot dendrogram
      message("Plotting clustering dendrogram")
      dendro <- suppressMessages(
        ggdendrogram(clust_genos,rotate = T,labels = T)+
          theme_dendro()+
          scale_y_reverse()
      )
      
      dendro_labels <- clust_genos$labels[clust_genos$order]
    }
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
                           index=1:length(unlist(popmap2)[unlist(popmap2) %in% dendro_labels]))
    name_order<-name_order[order(name_order$names),]
  }
  
  # Apply our name_order indices
  if(!plot_allele_frequency){
    
    # Filter in case we have lost clustered inds
    genos <- genos[genos$index %in% name_order$names,]
    genos$index <- as.character(genos$index)
    
    # Store the new indices in a dummy variable
    index2 <- rep(NA,nrow(genos))
    for(i in 1:nrow(name_order)){
      #genos[genos$index == name_order$names[i],"index"] <- as.character(name_order$index[i])
      index2[which(genos$index == name_order$names[i])] <- as.character(name_order$index[i])
    }
    
    # And replace the index column with index2
    genos$index <- index2
    
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
    
    # # Catch bad labelling for individuals
    # if(length(unique(popmap[,2])) == nrow(popmap)){
    #   cutoffs[1,"labs2"] <- nrow(popmap)
    # }
    
    # Label again...
    # Catch bad labelling for individuals accounting for popmaps where first 2 rows are pops
    if(length(popmap2[[1]]) == 1 & length(popmap2[[1]]) == 1){
      cutoffs[1,"labs2"] <- nrow(popmap)
    }
    
    # Reformat genos
    genos$index <- as.numeric(genos$index)
  }
  
  # Plot genotypes based on clustering...
  if(!plot_allele_frequency){
    if(cluster == FALSE){
      if(plot_phased){
        message("Plotting haplotypes without clustering")
      } else {
        message("Plotting genotypes without clustering")
      }
      genotypes<-suppressMessages(
        ggplot(data=genos,aes(x=snp,y=index))+
          #geom_raster(aes(fill=factor(GT)))+
          geom_tile(aes(fill=factor(GT)))+
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
          scale_y_continuous(breaks = as.numeric(cutoffs$labs2),
                             labels = cutoffs$pops)+
          scale_x_continuous(expand = c(0, 0))
        
      )
      
      # Retun a null object for cluster_pca
      cluster_pca=NULL
      
      # Apply phasing...
      if(plot_phased){
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme[c(1,3)],name="Haplotype",
                                                                    breaks=c("0","1"),labels=c("REF","ALT")))
      } else if(is.haploid){
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme[c(1,3)],name="Haploid Genotype",
                                                                    breaks=c("0","2"),labels=c("REF","ALT")))                                                                    
      } else {
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme,name="Genotype",
                                                                    breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT")))
      }
      
    } else {
      if(plot_phased){
        message("Plotting haplotypes with clustering")
      } else {
        message("Plotting genotypes with clustering")
      }
      genotypes <- suppressMessages(
        ggplot(data=genos,aes(x=snp,y=index))+
          #geom_raster(aes(fill=factor(GT)))+
          geom_tile(aes(fill=factor(GT)))+
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
      )
      # Apply phasing...
      if(plot_phased){
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme[c(1,3)],name="Haplotype",
                                                                    breaks=c("0","1"),labels=c("REF","ALT")))
      } else if(is.haploid){
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme[c(1,3)],name="Haploid Genotype",
                                                                    breaks=c("0","2"),labels=c("REF","ALT")))                                                                    
      } else {
        genotypes <- suppressMessages(genotypes + scale_fill_manual(values=colour_scheme,name="Genotype",
                                                                    breaks=c("0","1","2"),labels=c("HOM REF","HET","HOM ALT")))
      }
      
    }
  } else {
    message("Plotting per-population allele frequencies")
    
    # Format for plotting
    AF_mat_long <- data.frame(AF_mat_long)
    AF_mat_long$pop_F <- factor(AF_mat_long$pop,levels=unique(popmap[,2]))
    
    # Plot
    genotypes <- suppressMessages(
      ggplot(data=data.frame(AF_mat_long),aes(x=snp_pos,y=pop_F,fill=AF))+
        #geom_raster(aes(fill=factor(GT)))+
        geom_tile()+
        scale_fill_viridis(option="B")+
        theme_bw()+
        theme(panel.grid = element_blank(), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(), 
              axis.ticks  = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size=16),
              legend.position = "bottom",
              legend.title=element_text(size=14,vjust=1),
              legend.text = element_text(size=14,angle=45,hjust=1),
              panel.border = element_blank())+
        scale_x_continuous(expand = c(0, 0))
    )
    
    # Retun a null object for cluster_pca
    cluster_pca=NULL
  }
  
  
  # Return everything for plots...
  output<-list(lines,genotypes,dendro,dendro_labels,cluster_pca)
  names(output) <- c("positions","genotypes","dendrogram","dendro_labels","cluster_pca")
  return(output)
  
}
