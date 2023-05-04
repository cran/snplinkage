### R code from vignette source 'snplinkage.Rnw'

###################################################
### code chunk number 1: snplinkage.Rnw:45-46
###################################################
  options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: qc
###################################################
  library('snplinkage')

  gds_path <- save_hgdp_as_gds()
  gdata <- load_gds_as_genotype_data(gds_path)
  qc <- snprelate_qc(gdata, tagsnp = .99)
  print_qc_as_tex_table(qc)


###################################################
### code chunk number 3: 1p13
###################################################
  snp_idxs_1p13 <- select_region_idxs(qc$gdata,
    chromosome = 1, position_min = 114.4e6, n_snps = 20, offset = 12)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13, labels_colname = 'probe_id')
  grid::grid.draw(plt)


###################################################
### code chunk number 4: 1p13_large
###################################################
  snp_idxs_1p13_large <- select_region_idxs(qc$gdata, chromosome = 1,
    position_min = 114e6, n_snps = 100)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13_large)
  grid::grid.draw(plt)


###################################################
### code chunk number 5: 8p23
###################################################
  snp_idxs_8p23 <- select_region_idxs(qc$gdata, chromosome = 8,
    position_min = 11e6, position_max = 12e6)

  df_ld <- snprelate_ld(qc$gdata, snps_idx = snp_idxs_8p23, quiet = TRUE)
  plt <- gtable_ld(df_ld, df_snp = gdata_snps_annots(qc$gdata))
  grid::grid.draw(plt)


###################################################
### code chunk number 6: hla
###################################################
  snp_idxs_hla <- select_region_idxs(qc$gdata,
    chromosome = 6, position_min = 32.2e6, position_max = 32.8e6)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_hla)
  grid::grid.draw(plt)


###################################################
### code chunk number 7: hladr
###################################################
  snp_idxs_hladr <- select_region_idxs(qc$gdata,
    chromosome = 6, position_min = 32.5e6, n_snps = 20, offset = 9)

  # qc$gdata <- gdata_add_gene_annots(qc$gdata, snp_idxs_hladr)
  qc$gdata <- gdata_add_gene_annots_hladr_example(qc$gdata, snp_idxs_hladr)

  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_hladr, labels_colname = 'gene')
  grid::grid.draw(plt)


###################################################
### code chunk number 8: aim_assocs
###################################################
  snp_idxs_mhc <- select_region_idxs(qc$gdata,
    chromosome = 6, position_min = 29e6, position_max = 33e6)
  df_assocs <- chisq_pvalues_gdata(qc$gdata, snp_idxs_mhc)

  df_top_aim <- subset(df_assocs, rank(-pvalues, ties.method = 'first') <= 20)

  #qc$gdata <- gdata_add_gene_annots(qc$gdata, rownames(df_top_aim))
  qc$gdata <- gdata_add_gene_annots_aim_example(qc$gdata, rownames(df_top_aim))

  plt <- gtable_ld_associations_gdata(df_top_aim, qc$gdata,
    labels_colname = 'gene')
  grid::grid.draw(plt)


###################################################
### code chunk number 9: aim_assocs_large
###################################################
  plt <- gtable_ld_associations_gdata(df_assocs, qc$gdata,
    labels_colname = 'gene')
  grid::grid.draw(plt)


###################################################
### code chunk number 10: tagsnp
###################################################
  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13, r2 = 0.8)
  grid::grid.draw(plt)


###################################################
### code chunk number 11: tagsnp_large
###################################################
  plt <- gtable_ld_gdata(qc$gdata, snp_idxs_1p13_large, r2 = 0.8)
  grid::grid.draw(plt)


