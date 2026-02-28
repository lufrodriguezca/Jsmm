# R/merge_rc_data.R

#' Merge releases and capture data
#'
#' Sets releases and capture data.
#'
#' @param m Named list that contains all the model components for running Joint species movement modelling (Jsmm) analyses.
#' @return A modified Jsmm model m.
#' @export
#'

merge_rc_data = function(m){
  releases = m[["releases"]]
  captures = m[["captures"]]
  merged.releases = c()
  merged.captures = c()
  mult = c()
  j = 0
  for(i in 1:nrow(releases)){
    re = releases[i, ]
    ca = matrix(captures[captures[,1] == i, -1], ncol = 2)
    new.rc.history = TRUE
    same.release = which(re[1] == merged.releases[, 1] & re[2] == merged.releases[, 2])
    for(ind in same.release){
      if(!is.null(merged.captures)){
        old.ca = matrix(merged.captures[merged.captures[, 1] == ind, -1], ncol = 2)
      } else {
        old.ca = matrix(nrow = 0, ncol = 2)
      }
      if(nrow(old.ca) == nrow(ca)){
        if(all(old.ca == ca)){
          new.rc.history = FALSE
          mult[ind] = mult[ind] + 1
        }
      }
    }
    if(new.rc.history){
      j = j + 1
      mult[j] = 1
      merged.releases = rbind(merged.releases, re)
      if(nrow(ca) > 0) merged.captures = rbind(merged.captures, cbind(rep(j, nrow(ca)), ca))
    }
  }
  merged.releases = cbind(merged.releases, mult)
  rownames(merged.releases) = 1:nrow(merged.releases)
  colnames(merged.captures) = colnames(captures)
  m[["merged.releases"]] = merged.releases
  m[["merged.captures"]] = merged.captures
  class(m) = "class_m"
  return(m)
}
