# R/precomputation_KM.R

#' @noRd
#'
precomputation_KM = function(m){
  included = (m[["nc"]]>0)
  detB = m[["f_data"]][["det_B"]]
  el = m[["domain"]][["triangulation"]][["ele"]]
  x = m[["domain"]][["triangulation"]][["node"]]
  nnod = nrow(x)
  nele = nrow(el)

  metric = m[["f_data"]][["metric"]]

  ML = c(2,1,1, 1,2,1, 1,1,2)/24
  KxXX1 = c(1,-1,0, -1,1,0, 0,0,0)/6
  KxXY1 = c(1,0,-1, -1,0,1, 0,0,0)/6
  KxYY1 = c(1,0,-1, 0,0,0, -1,0,1)/6
  KyXX1 = c(1,-1,0, -1,1,0, 0,0,0)/6
  KyXY1 = c(1,0,-1, -1,0,1, 0,0,0)/6
  KyYY1 = c(1,0,-1, 0,0,0, -1,0,1)/6
  KXX1 = c(1,-1,0, -1,1,0, 0,0,0)/2
  KXY1 = c(1,0,-1, -1,0,1, 0,0,0)/2
  KYY1 = c(1,0,-1, 0,0,0, -1,0,1)/2
  KpX1 = c(-1,1,0, -1,1,0, -1,1,0)/6
  KpY1 = c(-1,0,1, -1,0,1, -1,0,1)/6
  KxpX1 = c(-1,1,0, -2,2,0, -1,1,0)/24
  KxpY1 = c(-1,0,1, -2,0,2, -1,0,1)/24
  KypX1 = c(-1,1,0, -1,1,0, -2,2,0)/24
  KypY1 = c(-1,0,1, -1,0,1, -2,0,2)/24
  Kxpp1 = c(2,2,1, 2,6,2, 1,2,2)/120
  Kypp1 = c(2,1,2, 1,2,2, 2,2,6)/120
  Kpp1 = c( 2,1,1, 1,2,1, 1,1,2)/24

  MT = matrix(0, ncol = nnod, nrow = nnod)
  xPos = matrix(0, ncol = 9, nrow = nele)

  for(e in 1:nele){
    ind = el[e, ]
    MT[ind, ind] = 1
  }
  MT = as(MT, "dgCMatrix")
  MTpos = which(matrix(MT, nrow = 1) > 0)

  pre11_01 = pre11_02 = pre11_03 = pre11_AC = pre22_01 = pre22_02 = pre22_03 =
  pre22_BC = preB1_01 = preB1_02 = preB1_03 = preB2_01 = preB2_02 = preB2_03 = rep(0,nele*9)

  for (e in 1:nele){
    ind = el[e, ]
    pos = c((ind - 1)*nnod + (ind[1] - 1), (ind - 1)*nnod + (ind[2] - 1), (ind - 1)*nnod + (ind[3] - 1)) + 1
    for(i in 1:9){
      xPos[e,i] =  which(MTpos == pos[i])
    }
    pre11_01[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,1]*metric[e,1]*KxXX1 + metric[e,1]*metric[e,2]*(KxXY1 + KxXY1[c(1,4,7,2,5,8,3,6,9)]) + metric[e,2]*metric[e,2]*KxYY1)
    pre11_02[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,1]*metric[e,1]*KyXX1 + metric[e,1]*metric[e,2]*(KyXY1 + KyXY1[c(1,4,7,2,5,8,3,6,9)]) + metric[e,2]*metric[e,2]*KyYY1)
    pre11_03[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,1]*metric[e,1]*KXX1 + metric[e,1]*metric[e,2]*(KXY1 + KXY1[c(1,4,7,2,5,8,3,6,9)] ) + metric[e,2]*metric[e,2]*KYY1)
    pre11_AC[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,1]*KpX1 + metric[e,2]*KpY1
    pre22_01[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,3]*metric[e,3]*KxXX1 + metric[e,3]*metric[e,4]*(KxXY1 + KxXY1[c(1,4,7,2,5,8,3,6,9)]) + metric[e,4]*metric[e,4]*KxYY1)
    pre22_02[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,3]*metric[e,3]*KyXX1 + metric[e,3]*metric[e,4]*(KyXY1 + KyXY1[c(1,4,7,2,5,8,3,6,9)]) + metric[e,4]*metric[e,4]*KyYY1)
    pre22_03[(9*(e-1)+1):(9*(e-1)+9)] = (metric[e,3]*metric[e,3]*KXX1 + metric[e,3]*metric[e,4]*(KXY1 + KXY1[c(1,4,7,2,5,8,3,6,9)]) + metric[e,4]*metric[e,4]*KYY1)
    pre22_BC[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,3]*KpX1 + metric[e,4]*KpY1
    preB1_01[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,1]*KxpX1 + metric[e,2]*KxpY1
    preB1_02[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,1]*KypX1 + metric[e,2]*KypY1
    preB1_03[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,1]*KpX1 + metric[e,2]*KpY1
    preB2_01[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,3]*KxpX1 + metric[e,4]*KxpY1
    preB2_02[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,3]*KypX1 + metric[e,4]*KypY1
    preB2_03[(9*(e-1)+1):(9*(e-1)+9)] = metric[e,3]*KpX1 + metric[e,4]*KpY1
  }

  vals = rep(0, length(MTpos))

  pre = list()
  pre[["xPos"]] =c(xPos)
  pre[["vals"]] =vals
  pre[["MTpos"]] = MTpos
  pre[["MT"]] = MT
  pre[["ML"]] = ML

  pre[["pre11_01"]] = pre11_01
  pre[["pre11_02"]] = pre11_02
  pre[["pre11_03"]] = pre11_03
  pre[["pre11_AC"]] = pre11_AC
  pre[["pre22_01"]] = pre22_01
  pre[["pre22_02"]] = pre22_02
  pre[["pre22_03"]] = pre22_03
  pre[["pre22_BC"]] = pre22_BC
  pre[["preB1_01"]] = preB1_01
  pre[["preB1_02"]] = preB1_02
  pre[["preB1_03"]] = preB1_03
  pre[["preB2_01"]] = preB2_01
  pre[["preB2_02"]] = preB2_02
  pre[["preB2_03"]] = preB2_03

  pre[["Kxpp"]] = Kxpp1
  pre[["Kypp"]] = Kypp1
  pre[["Kpp"]] = Kpp1

  return(pre)
}
