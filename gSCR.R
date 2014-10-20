

########################################################################################
########################################################################################
# parm: vector of starting values
# y: data (doesn't do occasion specific anything just yet!)
# K: number of occasions (this might assume the same number of occasions for every trap
#    but easy to check or change)
# Xdata: a dataframe with trap XYs and a column for each covariate
#        (X Y cords MUST be labeled X and Y)
# Gdata: a dataframe with state space XYs and a column for each covariate
#        (X Y cords MUST be labeled X and Y)
# Pformula: a R like formula argument for detection (e.g. ~1  or  ~vgHt + temp etc)
# Dformula: same  (names will of course match names in the Xdata and Gdata)
# grdSize: cell area this is just used to give back the total area of the state
#          space at the end
# ssBuffer: if no state space is provided one is created using a buffer of this value
#           around the outer most trap locations
# delta: and has a resolution defined by delta
# model: this doesnt exist anymore but I will integrate this again, basically this
#        argument lets you fit either a binomial B or Poisson P encounter model.
#        So right now its just the binomial
# sex: the sex (1 or 0) of each individual (related to rows in Y) which can include
#      unknow sex (NA). Not tried this without sex yet so thinking about it, I wonder
#      whither an additional if #statement needs to go in here to say if sex is NULL
#      ignore any sex stuff.
# sex.on: you can specify sex specific parameters
#         'n'one, 'int'ercept, 'sig'ma or 'b'oth
# predict: False fits the model, but you can supply point estimates of the parameters
#          to make predicted density surfaces of observed and unobserved individuals
#          under the model (see code #at the end it will be clear)
# cov: a RASTER! That will be used to estimate cost parameters - this is
#       the ecological distance bit
# directions: this is the number of cells used to calculate least cost paths
########################################################################################
########################################################################################


gSCR <- function(parm, y = y, K = NULL, Xdata = trapdata, Gdata = NULL,
                 Pformula=~1, Dformula=~1, grdSize=0.4, ssbuffer = 2, delta = 0.2,
                 model="B", sex = NULL, sex.on = c("n","int","sig","b"), predict=FALSE,
                 cov=NULL, directions=8){

  sex <- sex + 1
  X <- as.matrix(Xdata[,c("X","Y")])
  G <- as.matrix(Gdata[,c("X","Y")])
  nG <- nrow(G)
  SSarea <- grdSize * nG
  if(is.null(cov)){
    D <- e2dist(X, G)
  }
  if(!is.null(cov)){
    alpha2 <- exp(parm[length(parm)])
    cost <- exp(alpha2 * cov)
    tr <- transition(cost, transitionFunction=function(x) (1/(mean(x))), direction = directions)
    trLayer <- geoCorrection(tr, scl = F)
    D <- costDistance(trLayer,as.matrix(X),as.matrix(G))
  }
  n0 <- exp(parm[1])

  if(sex.on=="n"){
    alpha1 <- exp(parm[c(2,2)])
    alpha0 <- parm[c(3,3)]
    psi.sex <- plogis(parm[4])
    p1 <- 4
  }
  if(sex.on=="int"){
    alpha1 <- exp(parm[c(2,2)])
    alpha0 <- parm[c(3,4)]
    psi.sex <- plogis(parm[5])
    p1 <- 5
  }
  if(sex.on=="sig"){
    alpha1 <- exp(parm[c(2,3)])
    alpha0 <- parm[c(4,4)]
    psi.sex <- plogis(parm[5])
    p1 <- 5
  }
  if(sex.on=="b"){
    alpha1 <- exp(parm[c(2,3)])
    alpha0 <- parm[c(4,5)]
    psi.sex <- plogis(parm[6])
    p1 <- 6
  }

  if(Pformula==~1){
    pK <- 0
    probcap1 <- plogis(alpha0[1]) * exp(-alpha1[1] * D * D)
    probcap2 <- plogis(alpha0[2]) * exp(-alpha1[2] * D * D)
  } else{
    pK <- length(all.vars(Pformula))
    pparm <- parm[(p1+1):(p1+pK)]
    P_cov <- data.frame(Xdata[,which(colnames(Xdata) %in% all.vars(Pformula))])
    colnames(P_cov) <- all.vars(Pformula)
    modelP <- model.matrix(Pformula, P_cov)
    P_beta1 <- c(alpha0[1],pparm)
    P_beta2 <- c(alpha0[2],pparm)
    probcap1 <- c(plogis(modelP %*% P_beta1)) * exp(-alpha1[1] * D * D)
    probcap2 <- c(plogis(modelP %*% P_beta2)) * exp(-alpha1[2] * D * D)
  }

  if(Dformula==~1){
    D_cov <- NULL
    modelD <- NULL
   probs <- (1/nG)
  }else{
    dK <- length(all.vars(Dformula))
    dparm <- parm[(p1+pK+1):(p1+pK+dK)]
    D_cov <- data.frame(Gdata[,which(colnames(Gdata) %in% all.vars(Dformula))])
    colnames(D_cov) <- all.vars(Dformula)
    modelD <- model.matrix(Dformula, D_cov)
    D_beta <- c(0,dparm)
    probs <- exp(modelD %*% D_beta)
    probs <- probs/sum(probs)
  }

  Pm1 <- Pm2 <- matrix(NA, nrow = nrow(D), ncol = ncol(D))
  ymat <- y
  ymat <- rbind(ymat,rep(0,nrow(X)))
  lik.marg <- lik.marg1 <- lik.marg2 <- rep(NA,nrow(ymat))

  for (i in 1:nrow(ymat)){

    Pm1[1:length(Pm1)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap1[1:length(Pm1)], log =
    TRUE))
    lik.cond1 <- exp(colSums(Pm1))
    lik.marg1[i] <- sum(lik.cond1 * probs)

    Pm2[1:length(Pm2)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap2[1:length(Pm2)], log =
    TRUE))
    lik.cond2 <- exp(colSums(Pm2))
    lik.marg2[i] <- sum(lik.cond2 * probs)

   if( sex[i]==1 & !is.na(sex[i]) )
    lik.marg[i]<- lik.marg1[i] * (1-psi.sex)

   if( sex[i]==2 & !is.na(sex[i]) )
    lik.marg[i]<- lik.marg2[i] * psi.sex

   if( is.na(sex[i]) )
    lik.marg[i]<- lik.marg1[i] * (1-psi.sex) + lik.marg2[i]*psi.sex
  }

  if(predict==FALSE){
    nv <- c(rep(1, length(lik.marg) - 1), n0)
    part1 <- lgamma(nrow(y) + n0 + 1) - lgamma(n0 + 1)
    part2 <- sum(nv * log(lik.marg))
    out <- -1 * (part1 + part2)
    attr(out, "SSarea") <- SSarea
    return(out)
  }

  if(predict==TRUE){
    posterior<-matrix(NA,nrow=nG,ncol=nrow(ymat))
   for(i in 1:nrow(ymat)){
     Pm1[1:length(Pm1)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap1[1:length(Pm1)], log =
     TRUE))
     lik.cond1 <- exp(colSums(Pm1))
     Pm2[1:length(Pm2)] <- (dbinom(rep(ymat[i, ], nG), rep(K, nG), probcap2[1:length(Pm2)], log =
     TRUE))
     lik.cond2 <- exp(colSums(Pm2))

    if( sex[i]==1 & !is.na(sex[i]) )
     lik.cond <- lik.cond1 * (1-psi.sex)

    if( sex[i]==2 & !is.na(sex[i]) )
     lik.cond <- lik.cond2 * psi.sex

   if( is.na(sex[i]) )
    lik.cond <- lik.cond1 * (1-psi.sex) + lik.cond2*psi.sex

     posterior[,i]<- (lik.cond*probs)/lik.marg[i]
   }
    return(cbind(G,posterior))
  }
}

