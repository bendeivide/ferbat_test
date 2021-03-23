# Cenario 3: ET1 e Poder sob H0 parcial heterogeneo


# Simulacao do DIC em cinco cenarios
SimulateData <- function(n, r, delta, mu = 100, cvp = 0.10, cenario){
  sig <- cvp * mu
  trat <- 1:n
  trat <- as.factor(rep(trat, each = r))
  if (cenario == 1) { # H0 completa
    taui = 0 
    y <- mu + taui + rnorm(n * r, 0, sig)
  } else if (cenario == 2) { # H1 completa heterogeneo
    stderrdif <- sig / sqrt(r)
    ct <- 1:(n - 1)
    tau1 = 0
    tau <- c(tau1,tau1 + ct*delta*stderrdif)
    taui <- rep(tau, each = r)
    y <- mu + taui + rnorm(n * r, 0, sig)
  }   else if (cenario == 3) { # H0 parcial heterogeneo
    stderrdif <- sig / sqrt(r)
    k1 <- n %/% 2
    k2 <- n - k1        
    tau1 <- rep(c(0), times = k1)
    tau2 <- rep(c(delta*stderrdif), times = k2)
    taui <- rep(c(tau1, tau2), each = r)
    y <- mu + taui + rnorm(n * r, 0, sig)
  } else if (cenario == 4) { # H1 completa homogeneo
    stderrdif <- sig / sqrt(r)
    ct <- rep(delta, (n - 1))
    tau1 = 0
    tau <- c(tau1,tau1 + ct*stderrdif)
    taui <- rep(tau, each = r)
    y <- mu + taui + rnorm(n * r, 0, sig)
  } else if (cenario == 5) { # H0 parcial homogeneo
    stderrdif <- sig / sqrt(r)
    k1 <- n %/% 2
    k2 <- n - k1        
    ct <- rep(delta, k2)
    tau1 <- rep(c(0), times = k1)
    tau2 <- ct*stderrdif
    taui <- rep(c(tau1, tau2), each = r)
    y <- mu + taui + rnorm(n * r, 0, sig)
  } 
  y <- as.numeric(y)
  exper <- data.frame(trat, y)
  return(exper)
}   

# Analise de variancia
anava <- function(Y, controle = 1){    
  trat <- Y[,1] 
  av <- aov(Y[,2] ~ Y[,1])
  nu <- av$df.residual    
  QME <- deviance(av)/nu
  r <- 1/mean(1/table(trat)) 
  meanc <- mean(Y[,2][trat == controle])
  aux <- as.data.frame(subset(Y, Y[,1] != controle))
  meantrat <- tapply(aux[,2], aux[,1], mean)
  pos <- which(names(meantrat) == controle)
  meantrat <- meantrat[-pos]
  return(list(QME=QME, gle = nu, rh = r, meanc = meanc, meantrat = meantrat))
}  



Cen3Dun <- function(n, r, delta, k, alpha, prob = "ET1"){  
  
  # Simplificando a simulacao
  naux <- c(5, 10, 20, 40, 100)
  raux <- c(4, 10, 20)
  posn <- which(naux == n)
  posr <- which(raux == r)
  

  # Tam grupo I
  nk1 <- (n %/% 2)
  
  # Tm grupo II
  nk2 <- n - (n %/% 2)
  
  # Simulacao do experimento DIC:
  Y <- SimulateData(n = n, r = r, delta = delta, cenario = 3)
  av <- anava(Y, controle = 1)
  qme <- av$QME
  nu <- av$gle
  meanc <- av$meanc
  p <- 1 - alpha
  ri <- rep(av$rh, n)
  sd <- (2*qme/ri[1])^0.5 # Consideramos o mesmo numero de rep
  
  # # Dunnett
  library(rpgm)
  # q1 <- qDun(p, n, r, two.sided = TRUE)
  
  teste1 <- abs(av$meanc - av$meantrat)/sd # Teste Dunnett
  
  # Erro Tipo I
  if (prob == "ET1") {
    Decision <- any(teste1[1:(nk1 - 1)] >= q1[posn, posr])
  }
  # Poder
  if (prob == "POD") {
    teste1 <- teste1[nk1:(n - 1)]
    Decision <- as.vector(teste1[k] >= q1[posn, posr])
  }
  return(list(Decicion = Decision))
}

SimMCCen3 <- function(n, r, delta=1, mu =100, cvp=0.10, alpha=0.05, N=1000)
{  
  # k : números de delta's! 
  k <- c(1, 2, 4, 8, 16, 32)
  k <- k[k <= n - (n %/% 2)]
  nk <- length(k)
  Result <- matrix(0, nk + 1, 2)
  rownames(Result) <- paste(c("alpha", rep("k", nk)), as.character(c(alpha, k)),sep = " : ")  
  colnames(Result) <- c("DunExato","DunBoot") 
  add <- 1 / N

  rDun <- function(ns, n, r, two.sided = TRUE){
    nu <- n * (r - 1)
    X <- foreach(i = 1:n, .combine = 'rbind') %dopar% rnorm(ns, 0, (1/r)^0.5)
    X <- t(X)
    X[,1:(n - 1)] <- (X[,1:(n - 1)] - X[,n])/(2 / r)^0.5 #- 2*meanc/(qme/ri[1:(rp1 - 1)] + qme/ri[rp1])^0.5
    X <- if (two.sided == TRUE) abs(X[,-n]) else X[,-n]
    Max <- rowMaxs(X) / (rchisq(ns,nu)/nu)^0.5
    return(Max)
  } 
  qDun <- function(p, n, r, two.sided = TRUE){
    x <- rDun(ns = 1000000, n, r, two.sided) # mid
    qMax <- quantile(x, p) 
    return(qMax)
  }
  
  funaux <- function(teste1, anav, Result, q1, k, nk){
    #############
    # Erro tipo I
    #############
    if (any(teste1[1:(n %/% 2 - 1)] >= q1)) {
      Result[1,2] <- Result[1,2] + add
    } else Result[1,2] <- 0
    library(multcomp)
    t.dunnet <- glht(anav, linfct = mcp(trat = "Dunnet"))
    pvalue <- summary(t.dunnet)$test$pvalues
    if (any(pvalue[1:(n %/% 2 - 1)] <= alpha)) {
      Result[1,1] <- Result[1,1] + add
    } else Result[1,1] <- 0
    
    #######
    # Poder
    #######
    funaux2 <- function(x, q1){
      if (x >= q1) add else 0
    }
    funaux3 <- function(x, alpha){
      if (x <= alpha) add else 0
    }
    
    # Criterio MR Dunnett
    teste <- teste1[(n %/% 2):(n - 1)]
    teste <- teste[k]
    Result[2:(nk + 1),2] <- apply(teste, 1, funaux2, q1)
    
    # Criterio Dunnett
    library(multcomp)
    t.dunnet <- glht(anav, linfct = mcp(trat = "Dunnet"))
    pvalue <- pvalue[((n %/% 2)):(n - 1)]
    pvalue <- array(pvalue[k])
    Result[2:(nk + 1),1] <- apply(pvalue, 1, funaux3, alpha)
    
    # Resultado da simulacao
    return(Result)
  }
  SimulateData <- function(n, r, delta, mu, cvp, cenario){
    sigma <- cvp * mu
    trat <- 1:n
    trat <- as.factor(rep(trat, each = r))
    if (cenario == 1) 
    {
      taui = 0 
      y <- mu + taui + rnorm(n * r, 0, sigma)
    } else if (cenario == 2) 
    {
      stderrdif <- sigma / sqrt(r)
      ct <- 1:(n - 1)
      tau1 = 0
      tau <- c(tau1,tau1 + ct*delta*stderrdif)
      taui <- rep(tau, each = r)
      y <- mu + taui + rnorm(n * r, 0, sigma)
    }   else if (cenario == 3) 
    {
      stderrdif <- sigma / sqrt(r)
      k1 <- n %/% 2
      k2 <- n - k1        
      tau1 <- rep(c(0), times = k1)
      tau2 <- rep(c(delta*stderrdif), times = k2)
      taui <- rep(c(tau1, tau2), each = r)
      y <- mu + taui + rnorm(n * r, 0, sigma)
    } else if (cenario == 4) 
    {
      stderrdif <- sigma / sqrt(r)
      ct <- rep(delta, (n - 1))
      tau1 = 0
      tau <- c(tau1,tau1 + ct*stderrdif)
      taui <- rep(tau, each = r)
      y <- mu + taui + rnorm(n * r, 0, sigma)
    }  
    y <- as.numeric(y)
    exper <- data.frame(trat, y)
    return(exper)
  }  
  anava <- function(Y, controle = 1){    
    trat <- Y[,1] 
    av <- aov(Y[,2] ~ Y[,1])
    nu <- av$df.residual    
    QME <- deviance(av)/nu
    r <- 1/mean(1/table(trat)) 
    meanc <- mean(Y[,2][trat == controle])
    aux <- as.data.frame(subset(Y, Y[,1] != controle))
    meantrat <- tapply(aux[,2], aux[,1], mean)
    pos <- which(names(meantrat) == controle)
    meantrat <- meantrat[-pos]
    return(list(QME=QME, gle = nu, rh = r, meanc = meanc, meantrat = meantrat))
  }  
  aux <- foreach(i = 1:N, .combine = "+") %dopar% {
    Y <- SimulateData(n, r, delta, mu,cvp, 3)
    av <- anava(Y, controle = 1)
    qme <- av$QME
    nu <- av$gle
    meanc <- av$meanc
    p <- 1 - alpha
    ri <- rep(av$rh, n)
    sd <- (2*qme/ri[1])^0.5 # Consideramos o mesmo numero de rep
    
    # Dunnett bootstrap
    library(rpgm)
    library(doParallel)
    q1 <- qDun(p, n, r)
    teste1 <- abs(meanc - av$meantrat)/sd
    
    # Dunnet original
    anav <- aov(y~trat, Y)
    
    # Resultado
    funaux(teste1, anav, Result, q1, k, nk)
  }  
  return(aux)
}  
