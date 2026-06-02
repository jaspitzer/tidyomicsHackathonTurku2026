adjust_effect_sizes <- function(SE, contrast = NULL){
  deseq_object <- metadata(SE)$tidybulk$DESeq2_object #legibility
  formula_terms <- tidybulk:::parse_formula(BiocGenerics::design(deseq_object)) # legibility
  
  if(!is.null(contrast)){
    .coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
  }else if(sum(stringr::str_detect(names(rowData(SE)), paste0("___", formula_terms, " .+\\-"))) > 1){
    
  }else if(is.null(contrast)){
    factor_levels <- SE@colData[, formula_terms[1]] |> 
      as.factor() |> 
      levels()
    contrast = c(formula_terms[1], 
                 factor_levels[2], 
                 factor_levels[1])
    .coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
  }
  
  adjust <- DESeq2::lfcShrink(dds = deseq_object, 
                              res  = DESeq2::DESeqResults(BiocGenerics::as.data.frame(BiocGenerics::metadata(SE)$tidybulk$DESeq2_fit)), 
                              coef = .coef,
                              apeAdapt = F
  )
  return(adjust)
}