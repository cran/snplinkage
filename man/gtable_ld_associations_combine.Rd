% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtable_associations.R
\name{gtable_ld_associations_combine}
\alias{gtable_ld_associations_combine}
\title{Build gtable by combining ggplots}
\usage{
gtable_ld_associations_combine(ggplots, diamonds)
}
\arguments{
\item{ggplots}{List of ggplots}

\item{diamonds}{Does the LD visualization use diamond-type layout}
}
\value{
gtable of ggplots
}
\description{
Build gtable by combining ggplots
}
\examples{

library(snplinkage)

# example rnaseq data frame, 20 variables of 20 patients
m_rna = matrix(runif(20 ^ 2), nrow = 20)

# pair-wise correlation matrix
m_ld = cor(m_rna) ^ 2

# keep only upper triangle and reshape to data frame
m_ld[lower.tri(m_ld, diag = TRUE)] = NA
df_ld = reshape2::melt(m_ld) |> na.omit()

# rename for SNPLinkage
names(df_ld) = c('SNP_A', 'SNP_B', 'R2')

# visualize with ggplot_ld
gg_ld = ggplot_ld(df_ld)

# let's imagine the 20 variables came from 3 physically close regions
positions = c(runif(7, 10e5, 15e5), runif(6, 25e5, 30e5),
              runif(7, 45e5, 50e5)) |> sort()

# build the dataframe
df_snp_pos = data.frame(position = positions)
df_snp_pos$label = c(rep('HLA-A', 7), rep('HLA-B', 6), rep('HLA-C', 7))

gg_pos_biplot = ggplot_snp_pos(df_snp_pos, labels_colname = 'label',
                               upper_subset = TRUE)

# let's assume HLA-B is more associated with the outcome than the other genes
pvalues = c(runif(7, 1e-3, 1e-2), runif(6, 1e-8, 1e-6), runif(7, 1e-3, 1e-2))
log10_pvals = -log10(pvalues)

# we can reuse the df_snp_pos object
df_snp_pos$pvalues = log10_pvals

# add the chromosome column
df_snp_pos$chromosome = 6

gg_assocs = ggplot_associations(df_snp_pos, labels_colname = 'label',
                                linked_area = TRUE, nudge = c(0, 0.5),
                                n_labels = 12)

l_ggs = list(pos = gg_pos_biplot, ld = gg_ld, pval = gg_assocs)
gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)
grid::grid.draw(gt_ld)

}
