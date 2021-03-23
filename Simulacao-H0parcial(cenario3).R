# Cenario 3: ET1 e Poder sob H0 parcial heterogeneo

SimMCCen3 <- function(n, r, delta, alpha, N, two.sided){  
  # k: numero multiplicado de delta
  k <- c(1, 2, 4, 8, 16, 32)
  k <- k[k <= n - (n %/% 2)]
  nk <- length(k)
  
  # Objeto que armazena os resultados da simulacap ET1 e poder
  Result <- matrix(0, nk + 1, 2)
  rownames(Result) <- paste(c("alpha", rep("k", nk)), as.character(c(alpha, k)),sep = " : ")  
  colnames(Result) <- c("Dunnett","Ferbat") 
  
  # Probabilidade do erro tipo I por experimento e do poder
  add <- 1 / N
  
  # Distribuicao do Teste Dunnett 
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
  
  # Distribuicao do teste Ferreira&Batista (Ferbat)
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
  
  # Funcao que calcular a probabilidade de ET1 por experimento e Poder em cada experimento
  funaux <- function(teste1, teste2, Result, q1, q2, k, nk){
    #funaux <- function(teste1, teste2, anav, Result, q1, q2, k, nk){
    #############
    # Erro tipo I
    #############
    # Teste Dunnett
    if (any(teste1[1:(n %/% 2 - 1)] >= q1)) {
      Result[1,1] <- Result[1,1] + add
    } else Result[1,2] <- 0
    
    # Teste Ferbat
    if (any(teste2[1:(n %/% 2 - 1)] >= q2)) {
      Result[1,2] <- Result[1,2] + add
    } else Result[1,2] <- 0
    
    #Teste Dunnett (usando o multcomp) # Foi usado apenas para comparar com as minhas rotinas
    # library(multcomp)
    # t.dunnet <- glht(anav, linfct = mcp(trat = "Dunnet"))
    # pvalue <- summary(t.dunnet)$test$pvalues
    # if (any(pvalue[1:(n %/% 2 - 1)] <= alpha)) {
    #   Result[1,1] <- Result[1,1] + add
    # } else Result[1,1] <- 0
    
    #######
    # Poder
    #######
    funaux2 <- function(x, q){
      if (x >= q) add else 0
    }
    # funaux3 <- function(x, alpha){
    #   if (x <= alpha) add else 0
    # }

    # Criterio Ferbat
    teste <- teste2[(n %/% 2):(n - 1)]
    teste <- teste[k]
    Result[2:(nk + 1),2] <- apply(teste, 1, funaux2, q = q2)
    
    # Criterio Dunnnett
    teste <- teste1[(n %/% 2):(n - 1)]
    teste <- teste[k]
    Result[2:(nk + 1),1] <- apply(teste, 1, funaux2, q = q1)
      
    # Criterio Dunnett
    # library(multcomp)
    # t.dunnet <- glht(anav, linfct = mcp(trat = "Dunnet"))
    # pvalue <- pvalue[((n %/% 2)):(n - 1)]
    # pvalue <- array(pvalue[k])
    # Result[2:(nk + 1),1] <- apply(pvalue, 1, funaux3, alpha)
    
    # Resultado da simulacao
    return(Result)
  }
  
  # Simulacao do DIC em cinco cenarios
  SimulateData <- function(n, r, delta, mu = 100, cvp = 0.10, cenario){
    sigma <- cvp * mu
    trat <- 1:n
    trat <- as.factor(rep(trat, each = r))
    if (cenario == 1) 
    {
      taui = 0 
      y <- mu + taui + rnorm(n * r, 0, sigma)
    } else if (cenario == 2) {
      stderrdif <- sigma / sqrt(r)
      ct <- 1:(n - 1)
      tau1 = 0
      tau <- c(tau1,tau1 + ct*delta*stderrdif)
      taui <- rep(tau, each = r)
      y <- mu + taui + rnorm(n * r, 0, sigma)
    }   else if (cenario == 3) {
      stderrdif <- sigma / sqrt(r)
      k1 <- n %/% 2
      k2 <- n - k1        
      tau1 <- rep(c(0), times = k1)
      tau2 <- rep(c(delta*stderrdif), times = k2)
      taui <- rep(c(tau1, tau2), each = r)
      y <- mu + taui + rnorm(n * r, 0, sigma)
    } else if (cenario == 4) {
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
  
  # Simulacao do ET1 por experimento e poder em paralelo
  aux <- foreach(i = 1:N, .combine = "+") %dopar% {
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
    q1 <- qDun(p, n, r, two.sided)
    teste1 <- abs(av$meanc - av$meantrat)/sd # Teste Dunnett

    # # Anava
    anav <- aov(y~trat, Y)
    
    # Teste Ferbat
    q2 <- qFerbat(p, n, r, two.sided)
    
    mataux <- matrix(Y$y, n, r, byrow = TRUE)
    range <- rowMaxs(mataux) - rowMins(mataux)
    meanrange <- mean(range)
    sdrange <- (meanrange * (2 / r)^0.5) / dn(r + 0.55)
    
    teste2 <- (abs(av$meanc - av$meantrat) / sdrange) 
    
    # Resultado
    funaux(teste1, teste2, Result, q1, q2, k, nk)
  }  
  return(aux)
}  




