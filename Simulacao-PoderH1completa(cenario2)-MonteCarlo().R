

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
  aux <- 2 * r * (1 + x^2) / (1 - x^2)^2 * aux
  return(sum(aux * w))
}



# Simulacao do experimento
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
anava <- function(Y, controle){    
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

# Distribuicao Ferbat
rFerbat <- function(n, r, ns, two.sided = TRUE){
  dn <- function(r, np = 64){
    #r <- r - 1
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
    aux <- 2 * r * (1 + x^2) / (1 - x^2)^2 * aux
    return(sum(aux * w))
  }

  X <- rpgm.rnorm(n * r * ns, nthreads = maxthreads())
  X <- matrix(X, ns * n, r)
  # Calculo da media das ranges
  range <- rowMaxs(X) - rowMins(X)
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
  sd <- (meanrange * (2 / r)^0.5) / dn(r)
  ru <- maxrmedia / sd
  
  return(ru)  
}

qFerbat <- function(p, n, r, ns = 50000, two.sided = TRUE){
  x <- rFerbat(n, r, ns = ns, two.sided) # mid
  qut <- quantile(x, p) 
  return(qut)
}


# Cenario de simulacao do teste Ferbat
Cen2Ferbat <- function(n, r, delta = 1, k = 1, alpha=0.05){  

    
    Y <- SimulateData(n = n, r = r, mu = 100, cvp = 0.10, delta = delta, cenario = 2)
    av <- anava(Y, controle = 1)
    qme <- av$QME
    ri <- rep(av$rh, n)
    
    # Teste Ferbat
    
    mataux <- matrix(Y$y, n, r, byrow = TRUE)
    range <- rowMaxs(mataux) - rowMins(mataux)
    meanrange <- mean(range)
    if (r > 10) {
      sdrange <- (meanrange * (2 / r)^0.5) / dn(r + 0.23 * r)
    }
    if (r <= 10) {
      sdrange <- (meanrange * (2 / r)^0.5) / dn(r + 0.1 * r)
    }
    teste2 <- (abs(av$meanc - av$meantrat) / sdrange) 
    
    # Decisao
    Decisao <- as.vector(teste2[k] > qFerbat(p = 1 - alpha, n, r))
  return(list(Decisao = Decisao))
}  

# Cenario de simulacao do teste Dunnett
Cen2Dun <- function(n, r, delta = 1, k = 1, alpha=0.05){  
  # Simplificando a simulacao
  naux <- c(5, 10, 20, 40, 100)
  raux <- c(4, 10, 20)
  posn <- which(naux == n)
  posr <- which(raux == r)
  
  if (alpha == 0.05) {
    #Quantis Ferbat
    q2 <- c(2.69, 2.52, 2.50, 2.87, 2.76, 2.75, 3.03, 2.96, 2.97, 3.19, 3.15, 
            3.17, 3.41, 3.39, 3.40)
    q2 <- matrix(q2, 5, 3, byrow = TRUE)
    #Quantis Dunnett
    q1 <- c(2.73, 2.53, 2.48, 2.86, 2.74, 2.71, 2.99, 2.93, 2.91, 3.13, 3.09, 
            3.08, 3.32, 3.30, 3.30)
    q1 <- matrix(q1, 5, 3, byrow = TRUE)
  }
  if (alpha == 0.01) {
    # Quantis Ferbat (bilaterais) alpha = 0.01
    q2 <- c(3.48, 3.13, 3.07, 3.48, 3.29, 3.26, 3.56, 3.45, 3.44, 3.66, 3.59,
            3.61, 3.82, 3.80, 3.82)
    q2 <- matrix(q2, 5, 3, byrow = TRUE)
    
    # Quantis  Dunnett bilaterais (alpha = 0.01)
    q1 <- c(3.55, 3.17, 3.07, 3.53, 3.31, 3.26, 3.58, 3.46, 3.43, 3.67, 3.61, 
            3.59, 3.83, 3.80, 3.79)
    q1 <- matrix(q1, 5, 3, byrow = TRUE)
    
  }
    
    Y <- SimulateData(n = n, r = r, mu = 100, cvp = 0.10, delta = delta, cenario = 2)
    av <- anava(Y, controle = 1)
    qme <- av$QME
    nu <- av$gle
    meanc <- av$meanc
    p <- 1 - alpha
    ri <- rep(av$rh, n)
    sd <- (2*qme/ri[1])^0.5# Teste 3
    
    # Dunnett
    
    teste1 <- abs(meanc - av$meantrat)/sd # Teste Dunnett
    
    # Decisao
    Decisao <- as.vector(teste1[k] > q1[posn, posr])
  return(list(Decisao = Decisao))
}  