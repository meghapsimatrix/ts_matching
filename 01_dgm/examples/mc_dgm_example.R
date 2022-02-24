#####Set values#####

#Level 1: study/effect#
var.w <- .1                                #var w/in: constant
var.b <- 1/190                             #var btw: 1/190,1/40
rho <- var.b/(var.w+var.b)                 #ICC: .05,.20
muT <- 10                                  #Mean of treatment group: arbitrary
muC.C <- muT-.2*sqrt(var.w + var.b)        #Mean of control group in clustering: set to create a true mean ES of .2
muC.NC <- muT-.2*sqrt(var.w)               #Mean of control group in no clustering: set to create a true mean ES of .2
m <- 30                                    #Number of groups: constant (and same in T & C)
n <- 20                                    #Sample size w/in m: 20,100
Nt <- m*n
Nc <- m*n                                  #Assuming equal samples sizes in T & C.
N <- Nt + Nc

#Level 2: meta-analysis#
pcluster <- .10                            #prop of clustering: .10,.90
k <- 10                                    #Number of effects: 10,60
truevc <- .00                              #tau^2: 0,.05
truemu <- (muT-muC.C)/sqrt(var.w + var.b)  #True ES mean: constant at .2 [OR truemu = (muT-muC.NC)/sqrt(var.w) = .2]


#Sim#
reps <- 5000
output <- matrix(0,reps,24)
num.C <- k*pcluster
num.NC <- k - num.C


###Define models###

#RE model with the 'true' data#
likeRE.true <- function(pars) {
   vc <- pars[1]
   mu <- pars[2]
   thisone <- 1/2*log(var+vc) + 1/2*((d-mu)^2/(var+vc))
   return(sum(thisone)) }

#RE model with only naive data#
likeRE.naive <- function(pars) {
   vc <- pars[1]
   mu <- pars[2]
   thisone <- 1/2*log(vnaive+vc) + 1/2*((dnaive-mu)^2/(vnaive+vc))
   return(sum(thisone)) }

#RE model with only adjusted data#
likeRE.adj <- function(pars) {
   vc <- pars[1]
   mu <- pars[2]
   thisone <- 1/2*log(vT+vc) + 1/2*((dT-mu)^2/(vT+vc))
   return(sum(thisone)) }


###Start simulation###

for(z in 1:reps) {

#1) Generate clustered effects (d=dT) and calculate the appropriate variance (v=vT).
muC_grand_C <- rep(0,num.C); muT_grand_C <- rep(0,num.C)
var.pool_C <- rep(0,num.C); re.vc_C <- rep(0,num.C)

d_C <- rep(0,num.C); dnaive_C <- rep(0,num.C); dT_C <- rep(0,num.C)
v_C <- rep(0,num.C); vnaive_C <- rep(0,num.C); vT_C <- rep(0,num.C)

for(l in 1:num.C) {

#Generate two-group clustered raw data.
Mi_C <- rep(0,m); Yij_C <- matrix(0,m,n); muC_group <- rep(0,m)
Mi_T <- rep(0,m); Yij_T <- matrix(0,m,n); muT_group <- rep(0,m)

#Ys and means for clustered control group.
for(i in 1:m) { 
Mi_C[i] <- muC.C + sqrt(var.b)*rnorm(1,0,1)  #Qs: add vc?
for(j in 1:n) { Yij_C[i,j] <- Mi_C[i] + sqrt(var.w)*rnorm(1,0,1) }
muC_group[i] <- mean(Yij_C[i,])
}

#Ys and group means for clustered treatment group.
for(i in 1:m) { 
Mi_T[i] <- muT + sqrt(var.b)*rnorm(1,0,1)  #Qs: add vc?
for(j in 1:n) { Yij_T[i,j] <- Mi_T[i] + sqrt(var.w)*rnorm(1,0,1) }
muT_group[i] <- mean(Yij_T[i,])
}

#Calculate grand means and the total pooled within-groups variance.
muC_grand_C[l] <- mean(muC_group)
muT_grand_C[l] <- mean(muT_group)
var.pool_C[l] <- ( sum(rowSums((Yij_T - muT_grand_C[l])^2)) + sum(rowSums((Yij_C - muC_grand_C[l])^2)) ) / (N-2)

#Add random effects.
re.vc_C[l] <- rnorm(1,0,sqrt(truevc))

#Naive effects.
dnaive_C[l] <- ( (muT_grand_C[l] - muC_grand_C[l]) / sqrt(var.pool_C[l]) ) + re.vc_C[l]
vnaive_C[l] <- N/(Nt*Nc) + (dnaive_C[l]^2)/(2*(N-2))

#Adjusted effects.
dT_C[l] <- dnaive_C[l]*sqrt(1-((2*(n-1)*rho)/(N-2)))
h <- ((N-2)*((N-2)-2*(n-1)*rho)) / ((N-2)*((1-rho)^2)+n*(N-(2*n))*(rho^2)+2*(N-(2*n))*rho*(1-rho))
vT_C[l] <- (N/(Nt*Nc))*(1+(n-1)*rho) + (dT_C[l]^2)/(2*h)

#Truth = adjusted/clustered.
d_C[l] <- dT_C[l]
v_C[l] <- vT_C[l]
}


#2) Generate NOT clustered effects (d=dnaive) and calculate the appropriate variance (v=vnaive).
muC_grand_NC <- rep(0,num.NC); muT_grand_NC <- rep(0,num.NC)
var.pool_NC <- rep(0,num.NC); re.vc_NC <- rep(0,num.NC)

d_NC <- rep(0,num.NC); dnaive_NC <- rep(0,num.NC); dT_NC <- rep(0,num.NC)
v_NC <- rep(0,num.NC); vnaive_NC <- rep(0,num.NC); vT_NC <- rep(0,num.NC)

for(r in 1:num.NC) {

#Generate non-clustered raw data.
Yi_C <- rep(0,Nc); Yi_T <- rep(0,Nt)

#Ys for control group.
for(f in 1:Nc) { Yi_C[f] <- muC.NC + sqrt(var.w)*rnorm(1,0,1) }

#Ys for treatment group.
for(g in 1:Nt) { Yi_T[g] <- muT + sqrt(var.w)*rnorm(1,0,1) }

#Calculate grand means and the total pooled within-groups variance.
muC_grand_NC[r] <- mean(Yi_C)
muT_grand_NC[r] <- mean(Yi_T)
var.pool_NC[r] <- ( sum((Yi_T - muT_grand_NC[r])^2) + sum((Yi_C - muC_grand_NC[r])^2) ) / (N-2)

#Add random effects.
re.vc_NC[r] <- rnorm(1,0,sqrt(truevc))

#Naive effects.
dnaive_NC[r] <- ( (muT_grand_NC[r] - muC_grand_NC[r]) / sqrt(var.pool_NC[r]) ) + re.vc_NC[r]
vnaive_NC[r] <- N/(Nt*Nc) + (dnaive_NC[r]^2)/(2*(N-2))

#Adjusted effects.
dT_NC[r] <- dnaive_NC[r]*sqrt(1-((2*(n-1)*rho)/(N-2)))
h <- ((N-2)*((N-2)-2*(n-1)*rho)) / ((N-2)*((1-rho)^2)+n*(N-(2*n))*(rho^2)+2*(N-(2*n))*rho*(1-rho))
vT_NC[r] <- (N/(Nt*Nc))*(1+(n-1)*rho) + (dT_NC[r]^2)/(2*h)

#Truth = naive/not-clustered.
d_NC[r] <- dnaive_NC[r]
v_NC[r] <- vnaive_NC[r]
}


#3) Combine the C and NC data to create 3 data sets.
d <- c(d_C,d_NC); var <- c(v_C,v_NC)
dnaive <- c(dnaive_C,dnaive_NC); vnaive <- c(vnaive_C,vnaive_NC)
dT <- c(dT_C,dT_NC); vT <- c(vT_C,vT_NC)


#4) Estimate models.
pars <- c(truevc,truemu)

#TRUE#
REtrue <- nlminb(pars,likeRE.true,lower=c(0,-Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
#REtrue <- nlminb(pars,likeRE.true,lower=c(0,-Inf),upper=c(0,Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
trueest <- REtrue$par
output[z,1] <- trueest[1]
output[z,2] <- trueest[2]

SEREunadjustedest <- trueest
SEFEunadjustedest <- trueest[2]
yy <- d; y <- d; vv <- var; v <- var
if(trueest[1]==0) source("SEFE.unadjusted2.txt") else source("SERE.unadjusted2.txt")
if(trueest[1]==0) trueSEs <- FEunadjustedSEs else trueSEs <- REunadjustedSEs
output[z,3] <- trueSEs[1]
output[z,4] <- trueSEs[2]
output[z,5] <- REtrue$obj
output[z,6] <- REtrue$conv

#Naive#
REnaive <- nlminb(pars,likeRE.naive,lower=c(0,-Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
#REnaive <- nlminb(pars,likeRE.naive,lower=c(0,-Inf),upper=c(0,Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
naiveest <- REnaive$par
output[z,7] <- naiveest[1]
output[z,8] <- naiveest[2]

SEREunadjustedest <- naiveest
SEFEunadjustedest <- naiveest[2]
yy <- dnaive; y <- dnaive; vv <- vnaive; v <- vnaive
if(naiveest[1]==0) source("SEFE.unadjusted2.txt") else source("SERE.unadjusted2.txt")
if(naiveest[1]==0) naiveSEs <- FEunadjustedSEs else naiveSEs <- REunadjustedSEs
output[z,9] <- naiveSEs[1]
output[z,10] <- naiveSEs[2]
output[z,11] <- REnaive$obj
output[z,12] <- REnaive$conv

#Adjusted#
REadj <- nlminb(pars,likeRE.adj,lower=c(0,-Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
#REadj <- nlminb(pars,likeRE.adj,lower=c(0,-Inf),upper=c(0,Inf),control=list(eval.max=1000,iter.max=1000,abs.tol=10e-5,rel.tol=10e-5))
adjest <- REadj$par
output[z,13] <- adjest[1]
output[z,14] <- adjest[2]

SEREunadjustedest <- adjest
SEFEunadjustedest <- adjest[2]
yy <- dT; y <- dT; vv <- vT; v <- vT
if(adjest[1]==0) source("SEFE.unadjusted2.txt") else source("SERE.unadjusted2.txt")
if(adjest[1]==0) adjSEs <- FEunadjustedSEs else adjSEs <- REunadjustedSEs
output[z,15] <- adjSEs[1]
output[z,16] <- adjSEs[2]
output[z,17] <- REadj$obj
output[z,18] <- REadj$conv


#5) Evaluate heterogeneity: Q, I^2, & h.
w <- 1/var
FEmean <- sum(w*d)/sum(w)
Q <- sum(w*((d-FEmean)^2))
dfQ <- k-1
pQ <- 1-pchisq(Q,dfQ)
I2 <- max(0,100*((Q-dfQ)/Q))
H <- sqrt(Q/dfQ)
output[z,19] <- FEmean
output[z,20] <- Q
output[z,21] <- dfQ
output[z,22] <- pQ
output[z,23] <- I2
output[z,24] <- H


#6) Save raw data & output
write.table(rbind(d),file="d.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(var),file="var.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(dnaive),file="dnaive.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(vnaive),file="vnaive.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(dT),file="dT.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(vT),file="vT.txt",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(rbind(output[z,]),file="output.txt",append=TRUE,row.names=FALSE,col.names=FALSE)

if(z%%500==0) print(z) }
