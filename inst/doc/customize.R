## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  fig.width = 10,
  fig.height = 10,
  fig.align = 'center',
  dpi = 120,
  out.width = "100%",
  out.height = "100%",
  comment = "#>"
)

## -----------------------------------------------------------------------------
  library(snplinkage)

  # example rnaseq data matrix, 20 variables of 20 patients
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
  gg_ld

## -----------------------------------------------------------------------------
  # let's imagine the 20 variables came from 3 physically close regions
  positions = c(runif(7, 31e6, 31.5e6), runif(6, 32e6, 32.5e6),
                runif(7, 33e6, 33.5e6)) |> sort()

  # build the dataframe
  df_snp_pos = data.frame(position = positions)

  # minimal call
  gg_snp_pos = ggplot_snp_pos(df_snp_pos)

## -----------------------------------------------------------------------------
  df_snp_pos$label = c(rep('HLA-A', 7), rep('HLA-B', 6), rep('HLA-C', 7))
  gg_snp_pos = ggplot_snp_pos(df_snp_pos, labels_colname = 'label')

## -----------------------------------------------------------------------------
  l_ggs = list(snp_pos = gg_snp_pos, ld = gg_ld)
  gt_ld = gtable_ld_grobs(l_ggs, labels_colname = TRUE,
                          title = 'RNASeq correlations')
  grid::grid.draw(gt_ld)

## -----------------------------------------------------------------------------
  # let's imagine the middle region, HLA-B, is more associated with the outcome
  pvalues = c(runif(7, 1e-3, 1e-2), runif(6, 1e-8, 1e-6), runif(7, 1e-3, 1e-2))
  log10_pvals = -log10(pvalues)

  # we can reuse the df_snp_pos object
  df_snp_pos$pvalues = log10_pvals
  
  # add the chromosome column
  df_snp_pos$chromosome = 6

  gg_assocs = ggplot_associations(df_snp_pos, labels_colname = 'label',
                                  linked_area = TRUE, nudge = c(0, 0.5),
                                  n_labels = 12)

## -----------------------------------------------------------------------------
  gg_pos_biplot = ggplot_snp_pos(df_snp_pos, labels_colname = 'label',
                                 upper_subset = TRUE)

  # let's also say the middle region HLA-B is particularly correlated
  df_ld$R2[df_ld$SNP_A %in% 8:13 & df_ld$SNP_B %in% 8:13] = runif(15, 0.7, 0.9)
  gg_ld = ggplot_ld(df_ld)

  l_ggs = list(pos = gg_pos_biplot, ld = gg_ld, pval = gg_assocs)
  gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)

## -----------------------------------------------------------------------------
  library(ggplot2)
  gg_assocs <- gg_assocs + theme(axis.text.x = element_blank())
  title <- gg_assocs$labels$x %>% gsub(' (Mbp)', '', ., fixed = TRUE) %>%
    paste('-', nrow(df_snp_pos), 'SNPs')
  gg_assocs <- gg_assocs + labs(title = title, x = NULL)
  
  l_ggs$pval = gg_assocs
  gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)
  grid::grid.draw(gt_ld)

## -----------------------------------------------------------------------------
  gg_assocs$layers

## -----------------------------------------------------------------------------
  gg_assocs$layers[[1]]$aes_params$fill = "#0147ab"

## -----------------------------------------------------------------------------
  l_ggs$pval = gg_assocs
  gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)

## -----------------------------------------------------------------------------
  gg_pos_biplot = ggplot_snp_pos(df_snp_pos, labels_colname = 'label',
                                 upper_subset = TRUE, colors = '#101d6b')

  gg_assocs = ggplot_associations(df_snp_pos, labels_colname = 'label',
                                  linked_area = TRUE, nudge = c(0, 0.5),
                                  n_labels = 12, colors = '#101d6b')

  # extract title
  gg_assocs <- gg_assocs + theme(axis.text.x = element_blank())
  title <- gg_assocs$labels$x %>% gsub(' (Mbp)', '', ., fixed = TRUE) %>%
    paste('-', nrow(df_snp_pos), 'SNPs')
  gg_assocs <- gg_assocs + labs(title = title, x = NULL)

  # replace area color
  gg_assocs$layers[[1]]$aes_params$fill = "#0147ab"

  # rebuild
  l_ggs = list(pos = gg_pos_biplot, ld = gg_ld, pval = gg_assocs)
  gt_ld = gtable_ld_associations_combine(l_ggs, diamonds = TRUE)
  grid::grid.draw(gt_ld)

## -----------------------------------------------------------------------------
  data('crohn')
  m_hla = crohn[, -(1:6)]
  m_ld = cor(m_hla) ^ 2

  # keep only upper triangle and reshape to data frame
  m_ld[lower.tri(m_ld, diag = TRUE)] = NA
  df_ld = reshape2::melt(m_ld) |> na.omit()

  # rename for SNPLinkage
  names(df_ld) = c('SNP_A', 'SNP_B', 'R2')

  # visualize with ggplot_ld
  gg_ld = ggplot_ld(df_ld)

## -----------------------------------------------------------------------------
  mlog10_pvals = chisq_pvalues(m_hla, crohn[, 'crohn'])
  df_pos = data.frame(probe_id = colnames(m_hla), pvalues = mlog10_pvals,
                      chromosome = 5)

  # if we don't have positions we can use byindex = TRUE
  gg_assocs = ggplot_associations(df_pos, byindex = TRUE, nudge = c(0, 0.5))

## -----------------------------------------------------------------------------
  cowplot::plot_grid(gg_assocs, gg_ld, nrow = 2)

## -----------------------------------------------------------------------------
  df_top_assocs = subset(df_pos, pvalues > quantile(pvalues, 0.9))
  gg_assocs = ggplot_associations(df_top_assocs, linked_area = TRUE,
                                  nudge = c(0, 0.5))

  df_ld = subset(df_ld, SNP_A %in% df_top_assocs$probe_id &
                        SNP_B %in% df_top_assocs$probe_id)

  gg_ld = ggplot_ld(df_ld)
   
  cowplot::plot_grid(gg_assocs, gg_ld, nrow = 2)

## -----------------------------------------------------------------------------
  sessionInfo()

