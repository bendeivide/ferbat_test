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

qFerbat <- function(p, n, r, two.sided = TRUE){
  x <- rFerbat(n, r, ns = 50000, two.sided) # mid
  qut <- quantile(x, p) 
  return(qut)
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

# Distribuicao do teste Ferreira&Batista (Ferbat)
dn <- function(r, np = 64){
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


# Cenario de simulacao do teste Ferbat
Cen2Ferbat <- function(r, b){  
  n <- 5
  delta <- 2
  k <- 1
  alpha <- 0.05

    Y <- SimulateData(n = n, r = r, mu = 100, cvp = 0.10, delta = delta, cenario = 1)
    av <- anava(Y, controle = 1)
    qme <- av$QME
    ri <- rep(av$rh, n)
    
    # Teste Ferbat
    
    mataux <- matrix(Y$y, n, r, byrow = TRUE)
    range <- rowMaxs(mataux) - rowMins(mataux)
    meanrange <- mean(range)
    sdrange <- (meanrange * (2 / r)^0.5) / dn(r + r*b)
    teste2 <- (abs(av$meanc - av$meantrat) / sdrange) 
    
    # Decisao
    Decisao <- any(teste2> qFerbat(p = 1 - alpha, n, r))
  return(list(Decisao = Decisao))
}  

