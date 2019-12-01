### ================================================================== ###
### Modified versions of HW, AD, and DK tests from the homtest package ###
### Daniela Castro-Camilo                                              ###
### Email: daniela.castro.camilo@gmail.com                             ###
### ================================================================== ###
myHW.tests = function (x, cod, Nsim = 500){
  if (length(x) != length(cod)) {
    stop("x and cod must have the same length")
  }
  fac <- factor(cod)
  ni <- tapply(x, fac, length)
  k <- nlevels(fac)
  Lm <- sapply(split(x, fac), Lmoments)
  rLm <- regionalLmoments(x, fac)
  ti <- as.numeric(Lm[3, ])
  t3i <- as.numeric(Lm[4, ])
  t4i <- as.numeric(Lm[5, ])
  lambda1Reg <- as.numeric(rLm[1])
  lambda2Reg <- as.numeric(rLm[2])
  tauReg <- as.numeric(rLm[3])
  tau3Reg <- as.numeric(rLm[4])
  tau4Reg <- as.numeric(rLm[5])
  V1 <- (sum(ni * (ti - tauReg)^2)/sum(ni))^0.5
  V2 <- sum(ni * ((ti - tauReg)^2 + (t3i - tau3Reg)^2)^0.5)/sum(ni)
  parkappa <- par.kappa(lambda1Reg, lambda2Reg, tau3Reg, tau4Reg)
  xi <- parkappa$xi
  alfa <- parkappa$alfa
  kappa <- parkappa$k
  hacca <- parkappa$h
  V1s <- rep(NA, Nsim)
  V2s <- rep(NA, Nsim)
  for (i in 1:Nsim) {
    ti.sim <- rep(NA, k)
    t3i.sim <- rep(NA, k)
    for (j in 1:k) {
      campione <- rand.kappa(ni[j], xi, alfa, kappa, hacca)
      campione.ad <- campione/mean(campione)
      lmom <- Lmoments(campione.ad)
      ti.sim[j] <- lmom[3]
      t3i.sim[j] <- lmom[4]
    }
    tauReg.sim <- sum(ni * ti.sim)/sum(ni)
    tau3Reg.sim <- sum(ni * t3i.sim)/sum(ni)
    V1s[i] <- (sum(ni * (ti.sim - tauReg.sim)^2)/sum(ni))^0.5
    V2s[i] <- sum(ni * ((ti.sim - tauReg.sim)^2 + (t3i.sim - tau3Reg.sim)^2)^0.5)/sum(ni)
  }
  muV1 <- mean(V1s)
  stdV1 <- sd(V1s)
  muV2 <- mean(V2s)
  stdV2 <- sd(V2s)
  H1 <- (V1 - muV1)/stdV1
  H2 <- (V2 - muV2)/stdV2
  output <- c(H1, H2, tau3Reg.sim)
  names(output) <- c("H1", "H2", "t3R")
  return(output)
}

myADbootstrap.test = function (x, cod, Nsim = 500, index = 2, alpha = 0.05) {
  if (length(x) != length(cod)) {
    stop("x and cod must have the same length")
  }
  fac <- factor(cod)
  ni <- tapply(x, fac, length)
  k <- nlevels(fac)
  N <- sum(ni)
  if (index == 1) {
    indexflood <- function(x) {
      m <- mean(x)
      return(m)
    }
  }
  else if (index == 2) {
    indexflood <- function(x) {
      m <- median(x)
      return(m)
    }
  }
  med <- tapply(x, fac, indexflood)
  regione.adim <- x/unsplit(med, fac)
  A2kN <- ksampleA2(regione.adim, fac)
  A2kNs <- rep(NA, Nsim)
  for (i in 1:Nsim) {
    regione.simul <- nonparboot(regione.adim)
    med.simul <- tapply(regione.simul, fac, indexflood)
    regione.simul.adim <- regione.simul/unsplit(med.simul, fac)
    A2kNs[i] <- ksampleA2(regione.simul.adim, fac)
  }
  ecdfA2kNs <- ecdf(A2kNs)
  probabilita <- ecdfA2kNs(A2kN)
  output <- c(A2kN, probabilita)
  names(output) <- c("A2kN", "P")
  return(output)
}

myDK.test = function (x, cod, alpha = 0.05){
  if (length(x) != length(cod)) {
    stop("x and cod must have the same length")
  }
  fac <- factor(cod)
  Y <- x
  k <- nlevels(fac)
  nn <- tapply(Y, fac, length)
  N <- sum(nn)
  NN <- cumsum(nn)
  X <- matrix(0, nrow = N, ncol = k)
  X[1:nn[1], 1] <- 1
  for (j in 2:k) {
    X[(NN[j - 1] + 1):NN[j], j] <- 1
  }
  Z <- cbind(Y, X)
  U <- Z[(sort(Z[, 1], index.return = TRUE))$ix, ]
  U[, 1] <- (1:N)/N
  AA <- matrix(NA, nrow = max(nn), ncol = k)
  for (j in 1:k) {
    AA[1:nn[j], j] <- U[U[, j + 1] == 1, 1]
  }
  BB = cos(2 * pi * AA)
  CC = colSums(BB, na.rm = TRUE) * (2/nn)^0.5
  Ak = sum(CC^2)
  probabilita <- pchisq(Ak, df = k - 1, lower.tail = TRUE)
  if(probabilita < (1 - alpha)){print("DK accepts regional homogeneity")}else{print("DK rejects regional homogeneity")}
  output <- c(Ak, probabilita)
  names(output) <- c("Ak", "P")
  return(output)
}