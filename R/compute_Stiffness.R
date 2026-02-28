# R/compute_Stiffness.R

#' @noRd
#'
compute_Stiffness = function(m, XBeta, pre){
  included = (m[["nc"]] > 0)
  detB  = m[["f_data"]][["det_B"]]
  el    = m[["domain"]][["triangulation"]][["ele"]]
  x     = m[["domain"]][["triangulation"]][["node"]]
  nele  = nrow(el)

  Bj = m[["f_data"]][["B_j"]]
  C  = m[["f_data"]][["fix_vals"]][, 5]
  metric = m[["f_data"]][["metric"]]

  MT   = pre[["MT"]]
  xPos = pre[["xPos"]]
  ML   = pre[["ML"]]
  vals = pre[["vals"]]

  if(included["diffusion"]){
    XBeta_diffusion = XBeta[["diffusion"]]
  }else{
    XBeta_diffusion = XBeta[["diffusion"]]*0
  }
  if(included["advection_1"]){
    XBeta_advection_1 = XBeta[["advection_1"]]
  }else{
    XBeta_advection_1 = XBeta[["diffusion"]]*0
  }
  if(included["advection_2"]){
    XBeta_advection_2 = XBeta[["advection_2"]]
  }else{
    XBeta_advection_2 = XBeta[["diffusion"]]*0
  }
  if(included["mortality"]){
    XBeta_mortality = XBeta[["mortality"]]
  }else{
    XBeta_mortality = XBeta[["diffusion"]]*0
  }
  if(m[["observation_effort"]][["method"]] == "CCP" && included["observation"]){
    XBeta_mortality = XBeta_mortality + XBeta[["observation"]]
  }

  Ks = list()

  for(sp in 1:m[["ns"]]){
    vals = vals*0

    coefs_new_new = assemble_coefs(pqx = m[["f_data"]][["fix_vals"]][, 1],
                                  rqx = m[["f_data"]][["fix_vals"]][, 3],
                                  pqy = m[["f_data"]][["fix_vals"]][, 2],
                                  rqy = m[["f_data"]][["fix_vals"]][, 4],
                                  C =  m[["f_data"]][["fix_vals"]][, 5],
                                  method = m[["observation_effort"]][["method"]],
                                  XBeta_diffusion[, sp], XBeta_advection_1[, sp], XBeta_advection_2[, sp], XBeta_mortality[, sp],
                                  nele, node_x = m[["domain"]][["triangulation"]][["node"]][, 1],
                                  node_y = m[["domain"]][["triangulation"]][["node"]][, 2],
                                  el_1 = el[, 1])

    k = XBeta[["habitat_preference"]][, sp]

    vals = assemble_stiffness(k, ML, detB, vals, nele, xPos, coef= coefs_new_new, Bj= c(t(Bj)),
                              C,  node_x = m[["domain"]][["triangulation"]][["node"]][, 1],
                              node_y = m[["domain"]][["triangulation"]][["node"]][, 2], el_1 = el[, 1],
                              pre11_01 = pre[["pre11_01"]], pre11_02 = pre[["pre11_02"]], pre11_03 = pre[["pre11_03"]], pre11_AC = pre[["pre11_AC"]],
                              pre22_01 = pre[["pre22_01"]], pre22_02 = pre[["pre22_02"]], pre22_03 = pre[["pre22_03"]], pre22_BC = pre[["pre22_BC"]],
                              preB1_01 = pre[["preB1_01"]],  preB1_02 = pre[["preB1_02"]],  preB1_03 = pre[["preB1_03"]],
                              preB2_01 = pre[["preB2_01"]],  preB2_02 = pre[["preB2_02"]],  preB2_03 = pre[["preB2_03"]],
                              Kxpp = pre[["Kxpp"]],  Kypp = pre[["Kypp"]],  Kpp = pre[["Kpp"]])
    MT@x = vals
    Ks[[sp]] = MT
  }
  return (Ks)
}
