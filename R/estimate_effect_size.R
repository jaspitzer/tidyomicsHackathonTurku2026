#' Title
#'
#' @param SE A summarized experiment with calculated DE genes (DESseq2 only at the moment)
#' @param contrast A 3 part vector, analogous to contrast specification in test_differential_expression()
#'
#' @returns A DESeqResults instance with adjusted lfc fold changes
#' @export
#'
#' @examples
adjust_effect_sizes <- function(SE, contrast = NULL, prefix = NULL, method = "apeglm"){
  # currently this works for apeglm, but not apeAdapt, as that requires to pass
  # a check in Deseq2::lfcShrink that checks for the coefficient embedded in the
  # results meta data (not in object)
  deseq_object = S4Vectors::metadata(SE)$tidybulk$DESeq2_object #legibility
  formula_terms = tidybulk:::parse_formula(BiocGenerics::design(deseq_object)) # legibility
  
  # big issue: no clean way to distinguish various results; there is the prefix,
  # but that is not generally set. sensible default would be the most recent 
  # 
  
  if(!is.null(contrast)){
    .coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
  # }else if(sum(stringr::str_detect(names(rowData(SE)), paste0("___", formula_terms, " .+\\-"))) > 1){
  #   
  }else if(is.null(contrast)){
    factor_levels = SE@colData[, formula_terms[1]] |> 
      as.factor() |> 
      levels()
    contrast = c(formula_terms[1], 
                 factor_levels[2], 
                 factor_levels[1])
    .coef = paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
  }
  
  if(!is.null(prefix)){
    res_df = S4Vectors::metadata(SE)$tidybulk$DESeq2_fit %>% 
      dplyr::select(transcript, contains(prefix)) %>% 
      rename_with(\(x) str_remove(x, paste0("___", prefix)))
    results_object = res_df %>% 
      BiocGenerics::as.data.frame() %>% 
      DESeq2::DESeqResults()
  }else{
    results_object = DESeq2::DESeqResults(BiocGenerics::as.data.frame(S4Vectors::metadata(SE)$tidybulk$DESeq2_fit)) 
  }
  
  if(method == "apeglm"){
    adjust <- DESeq2::lfcShrink(dds = deseq_object, 
                                res  = results_object, 
                                coef = .coef,
                                apeAdapt = F)
  }else if(method == "ashr"){
    adjust <- DESeq2::lfcShrink(dds = deseq_object, 
                                res  = results_object, 
                                apeAdapt = F)
  }
  
  
  # right now it returns the object as is, goal would be to insert it into the 
  # rowdata and the results tibble in the meta data with an alternative column
  # name like "adjusted fold changes" or "estimated effects"
  return(adjust)
}