# R/convert_to_pslg.R

#' @noRd
#'
convert_to_pslg = function(lands_p){
  pol_nodes = matrix(unlist(lands_p[which(lands_p[["id"]] == 1), ][["geometry"]][[1]][1]), ncol = 2, byrow = FALSE)
  nnpol = nrow(pol_nodes)
  pol_nodes = pol_nodes[-nnpol, ]
  nnpol = nnpol - 1
  c = nnpol

  s_0 = 1:(nnpol)
  s_1 = 2:(nnpol + 1)
  s_1[nnpol] = 1

  v1_x = pol_nodes[, 1]
  v1_y = pol_nodes[, 2]

  npoly =  nrow(lands_p)

  for(k in 2:npoly){
    pol_nodes_a = matrix(unlist(lands_p[["geometry"]][k][[1]][1]), ncol = 2, byrow = FALSE)
    pol_nodes_a = pol_nodes_a[-nrow(pol_nodes_a), ]

    v2_x = pol_nodes_a[, 1]
    v2_y = pol_nodes_a[, 2]

    s2 = c()
    c2 = 0
    nnpol_a = nrow(pol_nodes_a)

    for(j in 1:nnpol_a){
      c2 = c2 + 1
      seg_x_ind = which(v1_x == pol_nodes_a[j, 1])
      seg_y_ind = which(v1_y == pol_nodes_a[j, 2])

      intl = length(intersect(seg_x_ind, seg_y_ind))

      if(intl == 1){
        s2[c2] = intersect(seg_x_ind, seg_y_ind)
      }else{
        c = c + 1
        s2[c2] = c
      }
    }

    seg_aux = matrix(NA, ncol = 2, nrow = (nnpol_a - 1))

    seg_aux[, 1] = s2[1:(nnpol_a - 1)]
    seg_aux[, 2] = s2[2:nnpol_a]

    seg_aux = rbind(seg_aux, c(s2[nnpol_a], s2[1]))

    sa_0 = seg_aux[, 1]
    sa_1 = seg_aux[, 2]

    s_0 = c( matrix(c(s_0,sa_0), nrow = 1, byrow = TRUE))
    s_1 = c( matrix(c(s_1,sa_1), nrow = 1, byrow = TRUE))

    for(i in 1:length(v1_x)){
      px = which(v2_x == v1_x[i])
      py = which(v2_y == v1_y[i])

      index = intersect(px, py)

      if(length(index) == 1){
        v2_x = v2_x[-index]
        v2_y = v2_y[-index]
      }
    }

    v1_x = c(v1_x, v2_x)
    v1_y = c(v1_y, v2_y)
    c = length(v1_x)
  }

  nrs = length(s_0)

  S = matrix(NA, nrow = nrs, ncol = 2)
  S[, 1] = s_0
  S[, 2] = s_1

  P = matrix(NA, nrow = length(v1_x), ncol = 2)
  P[, 1] = v1_x
  P[, 2] = v1_y

  SB = rep(0, nrs)

  for(i in 1:nrs){
    for(j in 1:nnpol){
      if(S[j, 1] == S[i, 1] && S[j, 2] == S[i, 2]){
        SB[i] = 1
      }
    }
  }

  PB = rep(0, length(v1_x))

  for(i in 1:length(v1_x)){
    if(i < (nnpol + 1)){
      PB[i] = 1
    }
  }

  pslg_format = RTriangle::pslg(P = P, PB = PB, S = S, SB = SB)

  return(pslg_format)
}
