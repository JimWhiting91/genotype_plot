context("Example plot")

test_that("readme code works", {
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "LL1"),
                           stringsAsFactors = FALSE)
  new_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",
                                                 package = "GenotypePlot"),
                            chr    = 1,
                            start  = 11697036,
                            end    = 11836450,
                            popmap = our_popmap,
                            cluster        = FALSE,
                            snp_label_size = 100000,
                            colour_scheme=c("#d4b9da","#e7298a","#980043"))
  expect_equal(class(new_plot), "list", label = "correct type");
  expect_equal(names(new_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(new_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$dendrogram), c("NULL"), label = "null type")
  expect_equal(class(new_plot$dendro_labels), c("NULL"), label = "null type")
  expect_equal(class(new_plot$cluster_pca), c("NULL"), label = "null type")
})

test_that("clustering works",{
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "LL1"),
                           stringsAsFactors = FALSE)
  cluster_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",
                                                 package = "GenotypePlot"),
                            chr    = 1,
                            start  = 11697036,
                            end    = 11836450,
                            popmap = our_popmap,
                            cluster        = TRUE,
                            snp_label_size = 100000,
                            colour_scheme=c("#d4b9da","#e7298a","#980043"))  
  
  expect_equal(class(cluster_plot), "list", label = "correct type");
  expect_equal(names(cluster_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(cluster_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(cluster_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(cluster_plot$dendrogram), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(cluster_plot$dendro_labels), c("character"), label = "character type")
  expect_equal(class(cluster_plot$cluster_pca), c("pca", "dudi"), label = "pca type")
})

test_that("multi.allelic works",{
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "LL1"),
                           stringsAsFactors = FALSE)
  
  multi_plot <- genotype_plot(vcf    = system.file("example_multi.vcf.gz",
                                                     package = "GenotypePlot"),
                                chr    = 1,
                                start  = 11697036,
                                end    = 11836450,
                                popmap = our_popmap,
                                cluster        = FALSE,
                                snp_label_size = 100000,
                                is.multi_allelic = TRUE,
                                colour_scheme=c("#d4b9da","#e7298a","#980043"))  
  
  expect_equal(class(multi_plot), "list", label = "correct type");
  expect_equal(names(multi_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(multi_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(multi_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(multi_plot$dendrogram), c("NULL"), label = "null type")
  expect_equal(class(multi_plot$dendro_labels), c("NULL"), label = "null type")
  expect_equal(class(multi_plot$cluster_pca), c("NULL"), label = "null type")
})

test_that("allele frequencies works",{
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "XAA"),
                           stringsAsFactors = FALSE)
  
  AF_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",
                                                   package = "GenotypePlot"),
                              chr    = 1,
                              start  = 11697036,
                              end    = 11836450,
                              popmap = our_popmap,
                              cluster        = FALSE,
                              snp_label_size = 100000,
                              plot_allele_frequency = TRUE,
                              colour_scheme=c("#d4b9da","#e7298a","#980043"))  
  
  expect_equal(class(AF_plot), "list", label = "correct type");
  expect_equal(names(AF_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(AF_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(AF_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(AF_plot$dendrogram), c("NULL"), label = "null type")
  expect_equal(class(AF_plot$dendro_labels), c("NULL"), label = "null type")
  expect_equal(class(AF_plot$cluster_pca), c("NULL"), label = "null type")
})

test_that("phasing works",{
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "LL1"),
                           stringsAsFactors = FALSE)
  
  phase_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",
                                                package = "GenotypePlot"),
                           chr    = 1,
                           start  = 11697036,
                           end    = 11836450,
                           popmap = our_popmap,
                           cluster        = TRUE,
                           snp_label_size = 100000,
                           plot_phased = TRUE,
                           colour_scheme=c("#d4b9da","#e7298a","#980043"))  
  
  expect_equal(class(phase_plot), "list", label = "correct type");
  expect_equal(names(phase_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(phase_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(phase_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(phase_plot$dendrogram), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(phase_plot$dendro_labels), c("character"), label = "character type")
  expect_equal(class(phase_plot$cluster_pca), c("pca", "dudi"), label = "pca type")
  
  # Add a check that phasing has doubled haplotypes...
  expect_equal(nrow(our_popmap),length(phase_plot$dendro_labels)/2,label="popmap length")
})

test_that("combining plots works",{
our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                         pop = c(rep("CEU", 5), "XAA", "LL1"),
                         stringsAsFactors = FALSE)
new_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",
                                               package = "GenotypePlot"),
                          chr    = 1,
                          start  = 11697036,
                          end    = 11836450,
                          popmap = our_popmap,
                          cluster        = TRUE,
                          snp_label_size = 100000,
                          colour_scheme=c("#d4b9da","#e7298a","#980043"))

combined_plot <- combine_genotype_plot(new_plot)
expect_equal(class(combined_plot), c("gg", "ggplot"), label = "ggplot type")
})

test_that("plotting haploid genotypes works",{
  our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                           pop = c(rep("CEU", 5), "XAA", "LL1"),
                           stringsAsFactors = FALSE)
  
  vcf_in <- vcfR::read.vcfR(system.file("example.vcf.gz",
                                        package = "GenotypePlot"))
  
  # Make it haploid
  gts <- extract.gt(vcf_in)
  gts[gts=="0|0"] <- "0"
  gts[gts=="0|1"] <- "0"
  gts[gts=="1|0"] <- "0"
  gts[gts=="1|1"] <- "1"
  gts_haps <- data.frame(FORMAT="GT",
                         gts)
  vcf_in@gt <- as.matrix(gts_haps)
  new_plot <- genotype_plot(vcf_object = vcf_in,
                            popmap = our_popmap,
                            cluster        = TRUE,
                            snp_label_size = 100000,
                            is.haploid = T,
                            colour_scheme=c("#d4b9da","#e7298a","#980043"))
  
  combined_plot <- combine_genotype_plot(new_plot)
  expect_equal(class(new_plot), "list", label = "correct type");
  expect_equal(names(new_plot), c("positions", "genotypes", "dendrogram", "dendro_labels","cluster_pca"), label = "list has correct names")
  expect_equal(class(new_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$dendrogram), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$dendro_labels), c("character"), label = "character type")
  expect_equal(class(new_plot$cluster_pca), c("pca", "dudi"), label = "pca type")
  
})

