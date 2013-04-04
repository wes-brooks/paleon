# 10/31
# Try a new way to calculate det(H1) and det(H0)

# setwd('C:/Users/JIN/Desktop/1028')

N=100                           # Number of loops
eps=0.05                        # TUNING 1 range of lambda
tau_Z=0.2                       # TUNING 2 SD of Z
I=2526
J=28

count.wi = read.csv("count.wi.csv", header=TRUE)  
index.wi = 1:4015 * (!is.na(count.wi[,4]))
cw = count.wi[index.wi,]

a1 = min(count.wi$Longitude)
a2 = max(count.wi$Longitude)
b1 = min(count.wi$Latitude)
b2 = max(count.wi$Latitude)

la = (a2-a1) / (length(unique(count.wi$Longitude)) - 1)
lb = (b2-b1) / (length(unique(count.wi$Latitude)) - 1)

#Neighbor structure
uv.wi = cbind(round((cw$Longitude-a1) / la + 1, 0), round((cw$Latitude-b1) / lb + 1, 0))
m_a = matrix(rep(uv.wi[,1], I), nrow=I)
m_b = matrix(rep(uv.wi[,2], I), nrow=I)
m.dis.abs = abs(m_a-t(m_a)) + abs(m_b-t(m_b))
C.abs = (m.dis.abs==1)                            
# sum(C.abs) = 9802

eigen.C = eigen(C.abs)$values
lambda_min = 1 / sort(eigen.C)[1]
lambda_max = 1 / sort(eigen.C)[I]

I_I = diag(1,I)
Y = as.matrix(cw[,4:31])

###########################
# Matrices to put results #
###########################

mu = sigma2 = lambda = array(0,c(J,N))
Z = array(0,c(I,J,N))

#Prior of mu
sigma2_mu = 10000                 

#Prior of sigma2
alpha0 = 1                             
beta0 = 1     

#Posterior of sigma2              
alpha1 = alpha0 + I/2                

#Prior of lambda
xi2 = 100                         
ratio.lambda = ratio.Z = 0

###################
# starting points #
###################

mu[,1]=apply(Y,2,sum)/sum(Y)
sigma2[,1]=1
lambda[,1]=0.1
Z[,,1]=log(1/J)

set.seed(1031)
date()

#########
# Loops #
#########

for(n in 2:N)
{
    cat(paste("n=", n, "\n", sep=''))

    for(j in 1:J)
    {
        cat(paste("j=", j, "\n", sep=''))

        ##########
        # 1.mu_j #
        ##########
     
        H0 = I_I - lambda[j,n-1] * C.abs
        sigma2.mu = (sum(H0)/sigma2[j,n-1] + 1/sigma2_mu)**(-1)
        mu.mu = sigma2.mu * sum(H0 %*% Z[,j,n-1])/sigma2[j,n-1]
        mu[j,n] = rnorm(1, mu.mu, sigma2.mu)
    
        ##############
        # 2.sigma2_j #
        ##############
    
        tmp = Z[,j,n-1] - mu[j,n]
        beta1 = beta0 + 0.5 * t(tmp) %*% H0 %*% tmp
        sigma2[j,n] = (rgamma(1, alpha1, rate=beta1))**(-1)
    
        ############
        # 3.lambda #
        ############
    
        lambda1 = runif(1, lambda[j,n-1]-eps, lambda[j,n-1]+eps)

        if((lambda1>lambda_min)&(lambda1<lambda_max))
        { 
            H1 = I_I-lambda1 * C.abs
            r = 0.5 * t(tmp) %*% (H1-H0) %*% tmp / sigma2[j,n] + (lambda1**2-lambda[j,n-1]**2)/xi2
            #ratio1 = det(H1)**0.5 / (det(H0)**0.5) * exp(-r)
            ratio1 = prod(1-lambda1*eigen.C)**0.5 / (prod(1-lambda[j,n-1]*eigen.C)**0.5) * exp(-r)

            if(runif(1)<ratio1)
            {
                lambda[j,n] = lambda1
                ratio1 = 1
            }
            else
            {
                lambda[j,n] = lambda[j,n-1]
                ratio1 = 0
            }
        } 
        
        if((lambda1<=lambda_min)|(lambda1>=lambda_max))
        {
            lambda[j,n] = lambda[j,n-1]
            ratio1=0
        }
    
        ratio.lambda = ratio.lambda + ratio1
    
    
    }
    
    ##########
    # 4.Z    #
    ##########
    
    Z[,,n] = Z[,,n-1]
    
    #Generate Z_ij one by one.
    #    for(i in 1:I)
    #    {
    #        for(j in 1:J)
    #        {
    #            Z[i,j,n] = rnorm(1, Z[i,j,n-1], sd=tau_Z)
    #            mu.z = mu[j,n]
    #
    #            for(k in which(C.abs[i,]!=0))
    #                mu.z = mu.z + lambda[j,n] * (Z[k,j,n]-mu[j,n])
    #
    #            r = 0.5 * ((Z[i,j,n-1]-mu.z)**2 - (Z[i,j,n]-mu.z)**2) / sigma2[j,n]
    #            pi1 = exp(Z[i,,n])  / sum(exp(Z[i,,n]))
    #            pi0 = exp(Z[i,,n-1]) / sum(exp(Z[i,,n-1]))
    #            ratio2 = prod(pi1**Y[i,]) / prod(pi0**Y[i,]) * exp(r)
    #
    #            if(runif(1)>ratio2)
    #            {
    #                Z[i,j,n] = Z[i,j,n-1]
    #                ratio2 = 0
    #            }
    #            else ratio2 = 1
    # 
    #            ratio.Z=ratio.Z+ratio2
    #        }
    #    }
    

    #Generate Z_i one by one.
    for(i in 1:I)
    {
        Z[i,,n] = rnorm(J, Z[i,,n-1], sd=tau_Z)
        mu.z = rep(mu[j,n], J)

        for(k in which(C.abs[i,]!=0))
            mu.z = mu.z + lambda[j,n] * (Z[k,,n] - mu[j,n])
      
        r = 0.5 * sum((Z[i,,n-1]-mu.z)**2 - (Z[i,,n]-mu.z)**2) / sigma2[j,n]
        pi1 = exp(Z[i,,n]) / sum(exp(Z[i,,n]))
        pi0 = exp(Z[i,,n-1]) / sum(exp(Z[i,,n-1]))
        ratio2 = prod(pi1**Y[i,]) / prod(pi0**Y[i,]) * exp(r)

        if(runif(1) > ratio2)
        {
            Z[i,,n] = Z[i,,n-1]
            ratio2 = 0
        }
        else ratio2 = 1
     
        ratio.Z = ratio.Z + ratio2
    }

}

date()
ratio.lambda/J/N
ratio.Z/I/N


##########################
# TRACE PLOT & HISTOGRAM # 
##########################

index = (N/5+1):N

jpeg(paste('1031.cw.mu1', N, eps, tau_Z,' .jpeg'), w=960, h=1280)
par(mfrow=c(7,8))
for(j in 1:J)
{
    plot(mu[j,], type='l', ylab=paste(colnames(Y)[j]), main=paste('mu_',j,sep=''), cex.main=1.5)
    abline(v=N/5, col=2)
    hist(mu[j,index], xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('mu_',j,sep=''), cex.main=1.5)
}
dev.off()

#Plot the difference between category j and 1.
jpeg(paste('1031.cw.mu2', N, eps, tau_Z, '.jpeg'), w=960, h=1280)
par(mfrow=c(7,8))
for(j in 2:J)
{
    plot(mu[j,]-mu[1,], type='l', ylab=paste(colnames(Y)[j]), main=paste('mu_',j,'-1',sep=''), cex.main=1.5)
    abline(v=N/5, col=2)
    hist((mu[j,index]-mu[1,index]), xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('mu_',j,'-1',sep=''), cex.main=1.5)
}
dev.off()


jpeg(paste('1031.cw.sigma2', N, eps, tau_Z, '.jpeg'), w=960, h=1280)
par(mfrow=c(7,8))
for(j in 1:J)
{
    plot(sigma2[j,], type='l', ylab=paste(colnames(Y)[j]), main=paste('sigma^2_',j,sep=''), cex.main=1.5)
    abline(v=N/5, col=2)
    hist(sigma2[j,index], xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('sigma^2_',j,sep=''), cex.main=1.5)
}
dev.off()


jpeg(paste('1031.cw.lambda', N, eps, tau_Z, '.jpeg'), w=960, h=1280)
par(mfrow=c(7,8))
for(j in 1:J)
{
    plot(lambda[j,], type='l', ylab=paste(colnames(Y)[j]), main=paste('lambda_',j,sep=''), cex.main=1.5)
    abline(v=N/5, col=2)
    hist(lambda[j,index], xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('lambda_',j,sep=''), cex.main=1.5)
}
dev.off()


I4=sample(1:I,4)
for(i in I4)
{
    jpeg(paste('1031.cw.Z1', N, eps, tau_Z, i, '.jpeg'), w=960, h=1280)
    par(mfrow=c(7,8))
    for(j in 1:J)
    {
        plot(Z[i,j,], type='l', ylab=paste(colnames(Y)[j]), main=paste('Z_',i,',',j,sep=''), cex.main=1.5)
        abline(v=N/5, col=2)
        hist(Z[i,j,index], xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('Z_',i,',',j, sep=''), cex.main=1.5)
    }
    dev.off()

}

I4=sample(1:I,4)
for(i in I4)
{
    #Plot the difference between category j and 1.
    jpeg(paste('1031.cw.Z2', N, eps, tau_Z, i, '.jpeg'), w=960, h=1280)
    par(mfrow=c(7,8))
    for(j in 2:J)
    {
        plot(Z[i,j,]-Z[i,1,],type='l',ylab=paste(colnames(Y)[j]),main=paste('Z_',i,',',j,'-1',sep=''),cex.main=1.5)
        abline(v=N/5, col=2)
        hist((Z[i,j,index]-Z[i,1,index]), xlab=paste(colnames(Y)[j]), freq=FALSE, main=paste('Z_',i,',',j,'-1',sep=''), cex.main=1.5)
    }
    dev.off()
}






















