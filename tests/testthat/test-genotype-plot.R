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
  expect_equal(names(new_plot), c("positions", "genotypes", "dendrogram", "dendro_labels"), label = "list has correct names")
  expect_equal(class(new_plot$positions), c("gg", "ggplot"), label = "ggplot type")
  expect_equal(class(new_plot$genotypes), c("gg", "ggplot"), label = "ggplot type")
})