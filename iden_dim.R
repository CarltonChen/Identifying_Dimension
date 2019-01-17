iden_dim = function(eig, d = 10, n, alpha = 0.1, pic = T){
  m1 = max(which(eig >= 0))
  # m = min(Tu, m1, n, d)
  m = min(Tu/2, m1, n/2)
  
  if(T){
    ratio = cumsum(eig[1:d]) / sum(eig[eig > 0])
    d1 = min(which(ratio >= 1-alpha))
  }
  
  if(T){
    # Eigen Ratio: Lam, C. and Yao, Q., 2012
    eig0 = eig[1:m]
    r_er = eig0[2:m] / eig0[1:(m - 1)]
    d2 = which.min(r_er[1:d])
  }
  
  if(T){
    # Growth Ratio: Ahn and Horenstein 2013
    s_eig = rev(cumsum(rev(eig[1:m])))
    gr1 = log(s_eig[2:(m-1)] / s_eig[3:m])
    gr2 = log(s_eig[1:(m-2)] / s_eig[2:(m-1)])
    r_gr = gr1 / gr2
    d3 = which.min(r_gr[1:d])
  }
  
  if(T){
    # Contribution Ratio: Xia, Liang, et.al 2018
    s_eig = rev(cumsum(rev(eig[1:m])))
    cr1 = eig[2:m] / s_eig[2:m]
    cr2 = eig[1:(m - 1)] / s_eig[1:(m - 1)]
    r_cr = cr1 / cr2
    d4 = which.min(r_gr[1:d])
  }
  
  if(T){
    # Information Criteria: Cho, H. et al. (2016) 
    
    iden_dim_ic = function(lambdas, d, n, pic = T){
      # Cho, H. et al. (2016) 
      
      lbd0 = lambdas
      
      tau_bot = 0
      tau_top = 1e8
      q0 = 0:(d-1) ############### d - 1
      lq = length(q0)
      lt = 101
      
      makeIC2 = function(lbd, d, c_star = 0, n) {
        IC2 = function(tau, q) {
          ic = log(c_star + sum(lbd[(q + 1):d]) / (d ^ 2)) + tau * q / sqrt(n)
          return(ic)
        }
      }
      
      # makeIC2 = function(lbd, d, n) {
      #   IC2 = function(tau, q) {
      #     ic = (sum(lbd[(q + 1):d]) / (d ^ 2)) + tau * q / sqrt(n)
      #     return(ic)
      #   }
      # }
      
      kk = 0
      tau_bot0 = tau_top0 = rep(NA, 10)
      while(kk <= 10){
        kk = kk + 1
        tau0 = seq(tau_bot, tau_top, length.out = lt)
        # lt = length(tau0)
        prs = expand.grid(tau0, q0)
        
        IC2 = makeIC2(lbd = lbd0, d = d, n = n)
        ics = mapply(IC2, prs$Var1, prs$Var2)
        ics = matrix(ics, lt)
        
        t_bot = max(which(apply(ics, 1, which.min) == d))
        t_top = min(which(apply(ics, 1, which.min) == 1))
        tau_bot0[kk] = tau0[t_bot]
        tau_top0[kk] = tau0[t_top]
        tau_bot = tau_bot0[kk] * 0.99 
        tau_top = tau_top0[kk] * 1.01
      }
      
      tau0 = seq(max(tau_bot0), min(tau_top0), length.out = lt)
      # lt = length(tau0)
      prs = expand.grid(tau0, q0)
      ics = mapply(IC2, prs$Var1, prs$Var2)
      ics = matrix(ics, lt)
      
      r = which.max(table(apply(ics, 1, which.min)))
      r = as.numeric(names(r)) - 1
      names(r) = NULL
      
      if(pic){
        plot(range(q0), range(ics), type = "n", xlab = "q", ylab = "IC(q)")
        for(i in 1:lt){
          if(which.min(ics[i,]) == r+1){
            col0 = "red"
          } else {
            col0 = 1
          }
          lines(q0, ics[i,], type = "o", pch = "+", col = col0, lty = 3)
        }
        abline(v = r, lty = 2, col = "red")
      }
      
      return(r)
      
    }
    d5 = iden_dim_ic(lambdas = eig, pic = pic)
  }
  
  return(c(d1, d2, d3, d4, d5))
}
