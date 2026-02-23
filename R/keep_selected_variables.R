# R/keep_selected_variables.R

#' @noRd
#'
keep_selected_variables = function(da, vars){
  da1 = da
  sel = which(colnames(da) %in% vars)
  nsel = length(sel)
  if(nsel == 0) da1 = NULL
  if(nsel == 1){
    da1 = data.frame(da[, sel])
    colnames(da1) = colnames(da)[sel]
  }
  if(nsel > 1){
    da1 = da[, sel]
  }
  return(da1)
}
