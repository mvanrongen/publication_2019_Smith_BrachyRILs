#' Genome scan using custom lmer formula
#'
#' @param genoprobs Genotype probabilities as calculated by calc_genoprob().
#' @param pheno A data.frame (does not have to be a numeric matrix).
#' @param map (optional) Genetic map of markers.
#' @param formula A formula of the null model (no genotype) using `lme4::lmer()`'s formula syntax.
#' @param null (optional) A formula of the null model to be fit. If provided `formula` is ignored. See details for how to specify
#' @param alternative (optional) A formula of the alternative model to be fit. If provided `formula` is ignored.
scan1_lmer <- function(cross, map, pheno, idcol, m0 = NULL, m1 = NULL){
  
  # Subset to retain only individuals with genotype data
  pheno <- pheno[pheno[[idcol]] %in% qtl2::ind_ids(cross), ]
  
  # define alternative model
  if(is.null(m1)){
    m1 <- paste0(format(terms(f)), " + GEN")
  }
  
  # fit null model
  fit_null <- lme4::lmer(m0, data = pheno)
  
  # get genotypes
  impute_genos <- viterbi(cross, map, error_prob = 0.01)
  all_gen <- impute_genos[pheno[[idcol]], ]
  
  # fit the model for each marker
  lod_scores <- lapply(all_gen, function(chr){
    lods <- lapply(dimnames(chr)[[2]], function(marker){
      
      # Get genotype matrix (as a factor)
      GEN <- factor(chr[, marker])
      
      tryCatch(
        {
          # Fit the genetic model
          fit_full <- lme4::lmer(m1, data = pheno)
          
          # Calculate LOD score
          LOD = (logLik(m1) - logLik(m0))/log(10)
          
        }, error=function(err) LOD = NA)
    })
    
    # turn to vector
    lods <- unlist(lods)
    names(lods) <- dimnames(chr)[[2]]
    
    return(lods)
  })
  
}
