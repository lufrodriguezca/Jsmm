# R/model_summary.R

#' Show summary of the Jsmm model.
#'
#' Generates a summary of all the Jsmm components.
#'
#' @param m named list corresponding to the Jsmm model.
#'
#' @return Produces a visualization on screen of the Jsmm main components.
#' @export
#'
model_summary = function(m){
  cat(paste("Domain: nele =", nrow(m[["domain"]][["triangulation"]][["ele"]]), "| nnod =", nrow(m[["domain"]][["triangulation"]][["node"]])), "\n")
  cat(paste("Capture method:",m[["observation_effort"]][["method"]], "| Number release sites:", length(unique(m[["observation_effort"]][["releases"]][["location"]])),
        "| Number capture sites:", length(unique(m[["observation_effort"]][["captures"]][["location"]]))), "\n")
  cat(paste("Number of time partitions:", length(m[["times"]])), "\n")
  cat("------------------------------------- \n")
  cat(paste("Model formulas:"),"\n")
  cat(paste("    Diffusion ~", as.character(m[["model_formula"]][["diffusion"]])[2], "| nc =", m[["nc"]][1]  ), "\n")
  if(!is.null(m[["model_formula"]][["advection"]][["b1"]])){
    cat(paste("    Advection b1 ~", as.character(m[["model_formula"]][["advection"]][["b1"]])[2], "| nc =", m[["nc"]][2] ), "\n")
  }
  if(!is.null(m[["model_formula"]][["advection"]][["b2"]])){
    cat(paste("    Advection b2 ~", as.character(m[["model_formula"]][["advection"]][["b2"]])[2], "| nc =", m[["nc"]][3]), "\n")
  }
  if(!is.null(m[["model_formula"]][["mortality"]])){
    cat(paste("    Mortality ~", as.character(m[["model_formula"]][["mortality"]])[2], "| nc =", m[["nc"]][4]), "\n")
  }
  if(!is.null(m[["model_formula"]][["habitat_preference"]])){
    cat(paste("    Habitat preference ~", as.character(m[["model_formula"]][["habitat_preference"]])[2], "| nc =", m[["nc"]][5]), "\n")
  }
  if(!is.null(m[["model_formula"]][["traits"]])){
    cat(paste("    Traits ~", as.character(m[["model_formula"]][["traits"]])[2] ), "\n")
  }
  if(!is.null(m[["model_formula"]][["observation"]])){
    cat(paste("    Observation ~", as.character(m[["model_formula"]][["observation"]])[2], "| nc =", m[["nc"]][6] ), "\n")
  }
  cat("------------------------------------- \n")
  cat(paste("Number of species:", m[["ns"]]),"\n")
  if(!is.null(m[["releases"]])){
    cat("Releases Data:\n")
    cat(paste("    ", nrow(m[["releases"]]), "unique releases" ),"\n")
  }
  if(!is.null(m[["captures"]])){
    cat("Captures Data:\n")
    cat(paste("    ", nrow(m[["captures"]]), "captures" ), "\n")
  }
}
