
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Genotype Plot

## Visualise and cluster genotypes based on a VCF

This function can be used to subset VCFs for regions and individuals of
interest and produce high-quality figures for publications.

The function depends on R being able to access a local installation of
bcftools in order to handle VCFs outside of R prior to reading in and
visualising. This functionality is only available therefore on UNIX
systems. However the function can be used on other systems by using the
`vcf_object` functionality (see below).

For a demo of the plotting options, see the shiny app:
<https://jimw91.shinyapps.io/genotype_plot_demo>

## Installation

``` r
install.packages("remotes")
remotes::install_github("JimWhiting91/genotype_plot")
```

## Tutorial

The package is really just a single function that handles everything. A
typical call requires the following coding (NOTE: you can just give the
full path to your vcf rather than use system.file):

``` r
library(GenotypePlot)

# popmap = two column data frame with column 1 for individual IDs as they appear in the VCF and column 2 for pop labels
our_popmap <- data.frame(ind = c("HG01914", "HG01985", "HG01986", "HG02013", "HG02051", "HG01879", "HG01880"),
                         pop = c(rep("CEU", 5), "XAA", "LL1"),
                         stringsAsFactors = FALSE)

# Make the genotype plot
new_plot <- genotype_plot(vcf    = system.file("example.vcf.gz",            # bgzipped VCF
                                               package = "GenotypePlot"),   
                          chr    = 1,                                       # chr or scaffold ID
                          start  = 11700000,                                # start of region
                          end    = 11800000,                                # end = end of region
                          popmap = our_popmap,                              # population membership
                          cluster        = FALSE,                           # whether to organise haplotypes by PCA clustering
                          snp_label_size = 10000,                          # breaks for position labels, eg. plot a position every 100,000 bp
                          colour_scheme=c("#FCD225","#C92D59","#300060"),   # character vector of colour values
                          invariant_filter = TRUE)                         # Filter any invariant sites before plotting
```

This can be minimally assembled into a plot for output.

``` r
combine_genotype_plot(new_plot)
```

This function uses `cowplot` to combine all the elements from the
genotype plot (genotypes and positions, or genotypes, positions, and
dendrogram if `cluster=TRUE`). The relative `widths` and `heights` flags
can be used to edit the width of the dendogram and height of the SNP
positions respectively.

Depending on the size of the region you’re looking at, these can take a
long time to plot within R, so I tend to plot them directly to a PDF,
e.g.

``` r
pdf("your_genotype_plot.pdf",width=10,height=8)
combine_genotype_plot(new_plot)
dev.off()
```

Genotypes are plotted according to the popmap order, so genotypes are
visualised within populations, labelled according to
`unique(popmap[,2])`, and are plotted in the order they appear in
`popmap[,1]`, unless `cluster=TRUE` (see below).

If `cluster=TRUE`, individuals are clustered according to PCA analysis
(`dudi.pca()`) and `hclust`. This is not designed to be an explicit test
of phylogeny, but can be useful to quickly visualise haplotype
relationships. If `cluster=TRUE`, haplotypes are no longer labelled
because the ordering is no longer user-defined but defined by the
clustering. New labels in clustered order are returned in the output as
`dendro_labels` (see Outputs). The clustering PCA output is also
returned as `cluster_pca` if needed, for example to plot PC scores
(`cluster_pca$li`)

The script first uses bcftools to subset the VCF based on the path given
and co-ordinates. This is simply so we’re not reading a VCF larger than
needs be into R.

``` r
# Subset our VCF
vcf_tempfile <- tempfile(pattern = "gt_plot", fileext = '.vcf')
system(glue("bcftools view -r {chr}:{start}-{end} {vcf} > {vcf_tempfile}"), wait=TRUE)
```

This command therefore writes a new file to the session temporary
directory, which is read in and then removed from the system using
`unlink`. The location of this directory is changed by modifying the
environment variable `TMPDIR`, usually in `~/.Renviron`.

After this, all plots are generated.

## No bcftools?

No problem. It is possible to use the package on vcfR objects directly
by reading any pre-made VCF into R with `my_vcf <-
read.vcfR("path/to/vcf")` and then parsing this object to
`genotype_plot()` with `vcf_object = my_vcf`. If doing this, the `chr`,
`start`, and `end` variables are set automatically using the chromosome
ID and min and max BP positions of the vcf object. For example:

``` r
new_plot <- genotype_plot(vcf_object  =  my_vcf,
                          popmap = our_popmap,                              
                          cluster        = FALSE,                           
                          snp_label_size = 10000,                          
                          colour_scheme=c("#d4b9da","#e7298a","#980043"))   
```

As for plotting VCFs generally, only vcfR objects with one chromosome or
scaffold ID are permitted.

This functionality can also be applied to more simply plotting VCFs
being manipulated/analysed in R with vcfR.

For more info on vcfR, see here:
<https://knausb.github.io/vcfR_documentation/>

## Inputs

The function handles the VCF outside of R using calls to `system()`, so
the input in terms of the VCF and region of interest are just character
strings, as described above.

#### Preparing the VCF

The VCF needs to be sorted, bgzipped and indexed prior to use. Sorting
can be done either `bcftools sort` or vcftool’s `vcf-sort` (if your VCF
was produced by stacks, you’ll need to use `vcf-sort`). To zip and
index, the following should work:

``` bash
bgzip -c your_data.vcf > your_data.vcf.gz
tabix -f -p vcf your_data.vcf.gz
```

For multi-allelic VCFs, these are coerced to bi-allelic VCFs with
BCFtools, such that multi-allelic variants are recoded as biallelic
variants at the same position. This is done based on the following:

``` bash
bcftools norm -m your_multi.vcf > your_bi.vcf
```

#### Popmap

The popmap should be a `data.frame` object with two columns: first
column = individual IDs as they appear in the VCF, and second column =
population label. The values from column 2 are used as labels in the
final plot. Column names are irrelevant, but the order must be column 1
for inds and column 2 for pop. This can either be made within R or read
in from a file with `read.table()`.

An example popmap should look like this:

``` r
> head(popmap)
     ind pop
1 LT_F16  TACHP
2 LT_F18  TACHP
3 LT_F19  TACHP
4 LT_F23  TACHP
5 LT_F24  TACHP
6  LT_M1  TACHP
```

The VCF is also filtered using the popmap, so you can read in a VCF with
many samples but only plot individuals in the popmap. When this is the
case, the VCF is also filtered for invariant sites between the remaining
individuals with a call to `vcfR::is.polymorphic()`.

If you want to produce a figure with a row per individual, rather than a
row with only a population label (as default), you can just give each
individual a unique value in the pop column for example to edit the
popmap above we can just do `popmap[,2] <- popmap[,1]` and produce a
popmap that looks like this:

``` r
> popmap[,2] <- popmap[,1]
> head(popmap)
      ind pop
1 LT_F16 LT_F16
2 LT_F18 LT_F18
3 LT_F19 LT_F19
4 LT_F23 LT_F23
5 LT_F24 LT_F24
6  LT_M1  LT_M1
```

## Alternative Plotting options

The easiest way to explore plotting options is to use the demo shiny app
at <https://jimw91.shinyapps.io/genotype_plot_demo>

#### Polarisation

Typically, genotypes in a VCF are polarised realtive to the reference
genome. However in some cases it can be desirable to re-polarise
genotypes to the major allele in a given population. This is done by
parsing the name of the focal population as a character string to
`polarise_genotypes` for eg. `polarise_genotypes="pop1"`.

#### Population Allele Frequencies

Per-population allele frequencies can be visualised by setting
`plot_allele_frequency=TRUE`. These are based on the frequency of
alleles within the populations defined in the popmap, providing a value
between 0 (absent) to 1 (fixed). This flag is compatible with
polarisation, so the frequency of the major allele within a focal
population can be visualised. This parameter is not compatible with
clustering, as the returning `genotypes` object only includes one row
per population, rather than one row per individual.

#### Phased haplotypes

If your data is phased, it can be desirable to visualise/cluster
haplotypes as opposed to genotypes. This can be done by setting
`plot_phased=TRUE`, but is contingent on genotypes being coded in phased
format (ie. 0|0, not 0/0). Phased genotypes are separated into two
haplotypes per individual, and are plotted in the first and last colour
given to `colour_scheme`. `plot_phased` is not compatible with
`plot_allele_frequency`.

## Outputs

The function returns a list where elements correspond to different parts
of the plots. As a standard, all plots return a `positions` and
`genotypes` element which correspond to the main genotype figure
(genotypes) and the genome position labels (positions).

If `cluster=TRUE`, the function also returns `dendrogram` and
`dendro_labels`, which correspond to a dendrogram of the clustered
haplotypes and the tip labels, respectively. The PCA output
`cluster_pca` used for clustering is also returned.

Each element is a ggplot object that can be modified as an individual
object in order for the user to modify any aspect of the plot as they
wish.

## Manipulating the dendrogram

The dendrogram is outputted with minimal formatting by default, but it
may be desirable to format this in such a way as to highlight
populations or individuals etc. The dendrogram object is just a ggplot
object made with the `ggdendro` package, so can be edited however you
wish (for examples, see
<http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning#ggdendro-package-ggplot2-and-dendrogram>).

For e.g., if we want to add back in the tip labels, we can plot as such

``` r
# Add dendrogram tips
dendro_with_tips <- new_plot$dendrogram +
                    geom_text(aes(x=1:length(new_plot$dendro_labels),
                    y=-2.5,
                    label=new_plot$dendro_labels))
```

Note here that `x` and `y` are inverted because the dendrogram has been
rotated 90 degrees. So here we are simply adding the
`new_plot$dendro_labels` back in at inverted `x` positions of 1 to
however many tips we have, and plotting them at an inverted `y` position
of -2.5 so they don’t overlap with the plot. In the example figure
below, where individuals are represented by points, this is done by
using the `new_plot$dendro_labels` to build a metadata `data.frame` in
which individuals have a Predation and River label that is then added in
a similar way to the above but with `geom_point()`.

Another way to sort of add tip labels to the dendrogram would be to run
`genotype_plot` twice, the first time with `cluster=TRUE` and use the
`output1$dendro_labels` to build a new popmap, e.g. `popmap2 <-
data.frame(ind=output1$dendro_labels,pop=output1$dendro_labels)`. Then
run again with the cluster-ordered `popmap2` in which individuals are
labelled in the pop column as individuals and set `cluster=FALSE`. You
can then plot the `output1$dendrogram` with `output2$genotypes` as
below.

Like all outputs from `genotype_plot()`, you can simply overwrite the
dendrogram slot with the new edited version and combine the plots as
normal.

## Change Log

    #### v0.2.1
      * is.haploid flag for handling VCFs with haploid genotypes
      * Bug fix for cases where VCF IDs are numbers rather than characters
      
    #### v0.2
      * Includes plotting options: phasing, allele frequencies, polarisation
      * Handling for multi-allelic VCFs
      * Missing data and invariant filtering
      * Update to clustering method based on PCA of genotype matrix
      * Updated documentation and shiny app demo
      
    #### v0.1
      * Individual filtering/reordering performed with bcftools, resolves issues with filtering vcf with vcfR.
      * Popmap is checked against the VCF prior to reading in, and errors are caught.
      * `vcf_object` flag added, can be used to read vcfR objected already in R and allows package to be used with windows
      * Error messages and general reporting fixes
      * `invariant_filter` flag added, so invariant filtering now optional but still default.
      * SNP marker labels improved so are added at regular intervals
    
    #### v0.0.9
    Original R package
