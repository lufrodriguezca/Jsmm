# R/print.class_m.R

#' @noRd
#'

print.class_m = function(x){
  cat(paste("Domain: nele =", nrow(x[["domain"]][["triangulation"]][["ele"]]), "| nnod =", nrow(x[["domain"]][["triangulation"]][["node"]])), "\n")
  cat(paste("Capture method:", x[["observation_effort"]][["method"]], "| Number release sites:", length(unique(x[["observation_effort"]][["releases"]][["location"]])),
            "| Number capture sites:", length(unique(x[["observation_effort"]][["captures"]][["location"]]))), "\n")
  cat(paste("Number of time partitions:", length(x[["times"]])), "\n")
  cat("------------------------------------- \n")
  cat(paste("Model formulas:"),"\n")
  cat(paste("    Diffusion ~", as.character(x[["model_formula"]][["diffusion"]])[2], "| nc =", x[["nc"]][1]  ), "\n")
  if(!is.null(x[["model_formula"]][["advection"]][["b1"]])){
    cat(paste("    Advection b1 ~", as.character(x[["model_formula"]][["advection"]][["b1"]])[2], "| nc =", x[["nc"]][2] ), "\n")
  }
  if(!is.null(x[["model_formula"]][["advection"]][["b2"]])){
    cat(paste("    Advection b2 ~", as.character(x[["model_formula"]][["advection"]][["b2"]])[2], "| nc =", x[["nc"]][3]), "\n")
  }
  if(!is.null(x[["model_formula"]][["mortality"]])){
    cat(paste("    Mortality ~", as.character(x[["model_formula"]][["mortality"]])[2], "| nc =", x[["nc"]][4]), "\n")
  }
  if(!is.null(x[["model_formula"]][["habitat_preference"]])){
    cat(paste("    Habitat preference ~", as.character(x[["model_formula"]][["habitat_preference"]])[2], "| nc =", x[["nc"]][5]), "\n")
  }
  if(!is.null(x[["model_formula"]][["traits"]])){
    cat(paste("    Traits ~", as.character(x[["model_formula"]][["traits"]])[2] ), "\n")
  }
  if(!is.null(x[["model_formula"]][["observation"]])){
    cat(paste("    Observation ~", as.character(x[["model_formula"]][["observation"]])[2], "| nc =", x[["nc"]][6] ), "\n")
  }
  cat("------------------------------------- \n")
  cat(paste("Number of species:", x[["ns"]]),"\n")
  if(!is.null(x[["releases"]])){
    cat("Releases Data:\n")
    cat(paste("    ", nrow(x[["releases"]]), "unique releases" ),"\n")
  }
  if(!is.null(x[["captures"]])){
    cat("Captures Data:\n")
    cat(paste("    ", nrow(x[["captures"]]), "captures" ), "\n")
  }
  #invisible(x)
}
