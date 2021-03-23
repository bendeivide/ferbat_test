# Cenario 2: Poder do teste em situacao heterogenea

SimMCCen2 <- function(n, r, delta=1, alpha=0.05, N=1000, )
{  
  # k : numeros de delta's! 
  k <- c(1, 2, 4, 8, 16, 32)
  k <- k[k <= (n - 1)]
  nk <- length(k)
  Result <- matrix(0, nk, 3)
  rownames(Result) <- paste("k", as.character(k),sep = " : ")
  colnames(Result) <- c("Dunnett_glht","Dunnett", "Ferbat")
  add <- 1 / N
  # Teste Dunnett 
  rDun <- function(ns, n, r, two.sided){
    nu <- n * (r - 1)
    X <- rpgm.rnorm(n * ns, 0, (1/r)^0.5, nthreads = maxthreads())
    X <- matrix(X, ns, n)    
    X[,1:(n - 1)] <- (X[,1:(n - 1)] - X[,n])/(2 / r)^0.5 
    X <- if (two.sided == TRUE) abs(X[,-n]) else X[,-n]
    Max <- rowMaxs(X) / (rchisq(ns,nu)/nu)^0.5
    return(Max)
  } 
  qDun <- function(p, n, r, two.sided){
    x <- rDun(ns = 100000, n, r, two.sided) # mid
    qMax <- quantile(x, p) 
    return(qMax)
  }
  
  # Distribuicao do teste Ferreira&Batista
  dn <- function(r, np = 50){
    GL <- function(size) {
      size <- as.integer(size)
      if (size < 0) 
        stop("Must be a non-negative number of nodes!")
      if (size == 0) 
        return(list(x = numeric(0), w = numeric(0)))
      i  <- 1:size
      j   <- 1:(size - 1)
      mu0 <- 2
      b <- j / (4 * j^2 - 1)^0.5
      A <- rep(0, size * size)
      A[(size + 1) * (j - 1) + 2] <- b
      A[(size + 1) * j] <- b
      dim(A) <- c(size, size)
      sd <- eigen(A, symmetric = TRUE)
      w <- rev(as.vector(sd$vectors[1, ]))
      w <- mu0 * w^2
      x <- rev(sd$values)
      return(list(nodes = x, weights = w))
    }
    
    xx <- GL(np)
    x  <- xx$nodes
    w  <- xx$weights
    
    y <- x / (1 - x^2)
    aux <- y * dnorm(y) * pnorm(y)^(r - 1)
    aux[aux <= 0] <- .Machine$double.eps^19 
    
    aux <- r * (1 + x)^2 / (1 - x^2)^2 * aux
    return(sum(aux * w))
  }
  rFerbat <- function(n, r, ns, two.sided){
    X <- rpgm.rnorm(n * r * ns, nthreads = maxthreads())
    X <- matrix(X, ns * n, r)
    # Calculo da media das ranges
    range <- rpgm::rowMaxs(X) - rpgm::rowMins(X)
    range <- matrix(range, n, ns)
    meanrange <- colMeans(range)
    
    # Calculo do maximo da diferenca entre as medias
    # e o controle
    media <- rowMeans(X)
    media <- matrix(media, n, ns)
    rmedia <- media[1:(n - 1),] - media[n,] 
    if (two.sided == TRUE) rmedia <- abs(rmedia) else rmedia
    maxrmedia <- colMaxs(rmedia)
    
    # Sdtudentizacao
    sd <- (meanrange/dn(r) * (2 / r)^0.5)
    ru <- maxrmedia / sd
    
    return(ru)  
  }
  qFerbat <- function(p, n, r, two.sided){
    x <- rFerbat(n, r, ns = 100000, two.sided) # mid
    qut <- quantile(x, p) 
    return(qut)
  }
  
  funaux2 <- function(x, q){
    if (x >= q) add else 0
  }
  funaux3 <- function(x, alpha){
    if (x <= alpha) add else 0
  }
  funaux <- function(teste1, teste2, anav, Result, q1, q2, k, nk){
    # Criterio do teste Dunnett
    teste <- teste1[k]
    Result[1:nk,2] <- apply(teste, 1, funaux2, q = q1)
    
    # Criterio do teste Dunnett
    teste <- teste2[k]
    Result[1:nk,3] <- apply(teste, 1, funaux2, q = q2)
    
    # Criterio Dunnett
    # library(multcomp)
    # t.dunnet <- glht(anav, linfct = mcp(trat = "Dunnet"))
    # pvalue <- as.numeric(summary(t.dunnet)$test$pvalues)[k]
    # pvalue <- array(pvalue)
    # Result[1:nk,1] <- apply(pvalue, 1, funaux3, alpha)
    
  return(Result)
  }
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
    return(list(QME = QME, gle = nu, rh = r, meanc = meanc, meantrat = meantrat))
  }  
  aux <- foreach(i = 1:N, .combine = "+") %dopar% {
    Y <- SimulateData(n = n, r = r, delta = delta, cenario = 2)
    av <- anava(Y, controle = 1)
    qme <- av$QME
    nu <- av$gle
    meanc <- av$meanc
    p <- 1 - alpha
    ri <- rep(av$rh, n)
    sd <- (2*qme/ri[1])^0.5# Teste 3
    
    # Dunnett
    library(rpgm)
    q1 <- qDun(p, n, r, two.sided = TRUE)
    teste1 <- abs(meanc - av$meantrat)/sd # Teste Dunnett
    
    
    # Anava
    anav <- aov(y~trat, Y)
    
    # Teste Ferbat
    q2 <- qFerbat(p, n, r, two.sided = TRUE)
    mataux <- matrix(Y$y, n, r, byrow = TRUE)
    range <- rowMaxs(mataux) - rowMins(mataux)
    meanrange <- mean(range)
    sdrange <- meanrange/dn(r) * (2 / r)^0.5
    teste2 <- (abs(meanc - av$meantrat) / sdrange) 
    
    # Resultado
    funaux(teste1, teste2, anav, Result, q1, q2, k, nk)
  }  
  return(aux)
}  




