library(tidyverse)

##########################################################################################
#   ______                                __        __                     
#  /      \                              /  |      /  |                    
# /$$$$$$  |  ______   ______    ______  $$ |____  $$/   _______   _______ 
# $$ | _$$/  /      \ /      \  /      \ $$      \ /  | /       | /       |
# $$ |/    |/$$$$$$  |$$$$$$  |/$$$$$$  |$$$$$$$  |$$ |/$$$$$$$/ /$$$$$$$/ 
# $$ |$$$$ |$$ |  $$/ /    $$ |$$ |  $$ |$$ |  $$ |$$ |$$ |      $$      \ 
# $$ \__$$ |$$ |     /$$$$$$$ |$$ |__$$ |$$ |  $$ |$$ |$$ \_____  $$$$$$  |
# $$    $$/ $$ |     $$    $$ |$$    $$/ $$ |  $$ |$$ |$$       |/     $$/ 
#  $$$$$$/  $$/       $$$$$$$/ $$$$$$$/  $$/   $$/ $$/  $$$$$$$/ $$$$$$$/  
#                             $$ |                                        
#                             $$ |                                        
#                             $$/  
# Plots the function and the area below it (or above), ideal for probability distributions
# col: color
graphics.areafunc <- function(func, x, col)
    {
        df = data.frame(x=x,y=func(x))
        df %>% ggplot(aes(x=x, y=y)) +
        geom_ribbon(data=df,aes(ymax=y),ymin=0,
                      fill=col,colour=NA,alpha=0.5)
    }

graphics.addplotdist <- function(dist,x,col1)
    {
        xp <- c(x[1],x,x[length(x)])
        yp <- c(0,dist,0)
    
        polygon(x = xp, y = yp, col = adjustcolor(col1,alpha.f=0.3), border = 0)
        lines(x,dist,col=col1,lwd=2)
    }

graphics.histandcurve <- function(dist, fun, col, xmin, xmax, breaks, title)
    {
        x <- seq(xmin,xmax, by = 1/(xmax - xmin))
        hist(dist,xlim=c(xmin -.75, xmax + .75),col=col, density = 30, main=title, 
             breaks = breaks, probability = TRUE)
        curve(fun(x), from=xmin, to=xmax, , xlab="x", ylab="y", col=col, lwd = 5, add = TRUE)
    }

graphics.2barplots <- function(dist1,dist2,x)
    {
        n1.col <- 'royalblue'
        n2.col <- 'steelblue4'
        ns <- rbind(dist1,dist2)
        barplot(ns,beside=T,names.arg=x,col=c(n1.col,n2.col),
                ylab='n observations', xlab='deaths',legend=TRUE)
    }

grafichs.printcredibilityinterval <- function(x,y,est,cred,xlim,col)
    {
        plot(x,y,pch='',xlim=xlim, main = 'Credibility interval')
    
        i <- 1
        j <- 1
        
        while(x[j] < cred[2])
            {
                if(x[i] < cred[1])
                    {
                        i <- i + 1
                    }
            
                j <- j + 1
            }
        k <- i
        yp <- c(0)
        while(k<=j)
            {
                yp <- c(yp,y[k])
                k <- k + 1
            }
        yp <- c(yp, 0)
        xp <- c(cred[1], seq(x[i], x[j], x[2]-x[1]), cred[2])
    
        polygon(x = xp, y = yp, col = adjustcolor(col,alpha.f=0.3), border = 0)
        lines(x,y,col=col,lwd=2)
        abline(v=est,lwd=2,col=col,lty=2)
        abline(v=cred[1],lwd=1,col=col,lty=2)
        abline(v=cred[2],lwd=1,col=col,lty=2)
    }

##########################################################################################
#  ______             ______                                                            
# /      |           /      \                                                           
# $$$$$$/  _______  /$$$$$$  |______    ______    ______   _______    _______   ______  
#   $$ |  /       \ $$ |_ $$//      \  /      \  /      \ /       \  /       | /      \ 
#   $$ |  $$$$$$$  |$$   |  /$$$$$$  |/$$$$$$  |/$$$$$$  |$$$$$$$  |/$$$$$$$/ /$$$$$$  |
#   $$ |  $$ |  $$ |$$$$/   $$    $$ |$$ |  $$/ $$    $$ |$$ |  $$ |$$ |      $$    $$ |
#  _$$ |_ $$ |  $$ |$$ |    $$$$$$$$/ $$ |      $$$$$$$$/ $$ |  $$ |$$ \_____ $$$$$$$$/ 
# / $$   |$$ |  $$ |$$ |    $$       |$$ |      $$       |$$ |  $$ |$$       |$$       |
# $$$$$$/ $$/   $$/ $$/      $$$$$$$/ $$/        $$$$$$$/ $$/   $$/  $$$$$$$/  $$$$$$$/ 

bayes.computeposterior <- function(likelihood, prior)
    {
        n.sample <- 2000
        delta.p  <- 1/n.sample
        p <- seq(from=1/(2*n.sample), by=1/n.sample, length.out=n.sample)
        
        likelihoodvec <- likelihood(p)
        priorvec <- prior(p)
    
        normalization <- delta.p * sum(likelihoodvec*priorvec)
        posterior <- likelihoodvec*priorvec/normalization
    
        return(posterior)
    }

bayes.getcredibilityinterval <- function(dist,conf,delta.p)
    {
        indx <- match(max(dist),dist)
        area <- max(dist)*delta.p
        
        i <- 1
        j <- 1
        while(area < conf)
            {
                if (indx - j > 0 )
                    {
                        area <- area + ( dist[indx+i] + dist[indx-j] ) * delta.p
                        i <- i + 1
                        j <- j + 1
                    }
                # This may solve some problems for distribution where the max is close to 0
                else
                    {
                        area <- area + ( dist[indx+i] ) * delta.p
                        i <- i + 1
                    }
            }
    
        return( c( p[indx-j] , p[indx+i] ) ) 
    }


##########################################################################################
#   ______                                     __  __                     
#  /      \                                   /  |/  |                    
# /$$$$$$  |  ______   _____  ____    ______  $$ |$$/  _______    ______  
# $$ \__$$/  /      \ /     \/    \  /      \ $$ |/  |/       \  /      \ 
# $$      \  $$$$$$  |$$$$$$ $$$$  |/$$$$$$  |$$ |$$ |$$$$$$$  |/$$$$$$  |
#  $$$$$$  | /    $$ |$$ | $$ | $$ |$$ |  $$ |$$ |$$ |$$ |  $$ |$$ |  $$ |
# /  \__$$ |/$$$$$$$ |$$ | $$ | $$ |$$ |__$$ |$$ |$$ |$$ |  $$ |$$ \__$$ |
# $$    $$/ $$    $$ |$$ | $$ | $$ |$$    $$/ $$ |$$ |$$ |  $$ |$$    $$ |
#  $$$$$$/   $$$$$$$/ $$/  $$/  $$/ $$$$$$$/  $$/ $$/ $$/   $$/  $$$$$$$ |
#                                  $$ |                        /  \__$$ |
#                                  $$ |                        $$    $$/ 
#                                  $$/                          $$$$$$/ 

# Parameters:
#            func: a function whose first argument is a real vector of parameters
#                  func returns a log10 of the likelihood function
#            theta.init: the initial value of the Markov Chain (and of func)
#            n.sample: number of required samples
#            sigma: standard deviation of the gaussian MCMC sampling pdf

metropolis.1dim <- function(func, theta.init, n.sample, sigma)
    {
        theta.cur <- theta.init        # current x value
        func.cur  <- func(theta.cur)   # current Log10(f(x)) value
    
        # We will return this matrix, it contains all the steps
        func.samp <- matrix(data=NA, nrow=n.sample, ncol=2)
        
        n.accept    <- 0
        rate.accept <- 0.0
    
        for(n in 1:n.sample)
            {
                # we pick a proposal sampling once from a gaussian distribution with mean theta.curr and
                # sigma as the sigma parameter
                theta.prop <- rnorm(n=1, mean = theta.cur, sigma)
                func.prop  <- func(theta.prop)
            
                # log10(Pprop/Pcurr) -> Log(Pprop) - Log(Pcurr)
                logMR <- func.prop - func.cur # Log10 of the metropolis ratio
                
                # We accept if
                # -> logMR > Log10(1) = 0
                # -> Log10(runif(1)) < logMR (accept reject)
                if( logMR>=0 || logMR>log10(runif(1)) )
                    {
                        # accept
                        theta.cur <- theta.prop
                        func.cur  <- func.prop
                        n.accept  <- n.accept + 1
                    }
                
                # Now we update our matrix with the n-th step:
                func.samp[n, 1] <- func.cur
                func.samp[n, 2] <- theta.cur
            }
    
        return(func.samp)
    }

#  This should work for every dim

# Parameters:
#            func: a function whose first argument is a real vector of parameters
#                  func returns a log10 of the likelihood function
#            theta.init: the initial value of the Markov Chain (and of func)
#            n.sample: number of required samples
#            sample.cov: standard deviations of the gaussian MCMC sampling pdf
metropolis.ndim <- function(func, theta.init, n.burnin, n.sample, sample.cov)
    {
        n.theta   <- length(theta.init)
    
        theta.cur <- theta.init
        func.cur  <- func(theta.init) # log10
        
        # We will return this matrix, it contains all the steps
        func.samp <- matrix(data=NA, nrow=n.sample, ncol=2+n.theta) 

        n.accept    <- 0
        rate.accept <- 0.0
  
        for(n in 1:(n.burnin+n.sample))
            {

                # Metropolis algorithm. No Hastings factor for symmetric proposal
                if(is.null(dim(sample.cov)))
                    { 
                        # theta and sampleCov are scalars
                        theta.prop <- rnorm(n=1, mean=thetaCur, sd=sqrt(sampleCov))
                    }
                else
                    {
                        theta.prop <- rmvnorm(n=1, mean=theta.cur, sigma=sample.cov, method="eigen")
                    }

                func.prop  <- func(theta.prop) 
            
                # log10(Pprop/Pcurr) -> Log(Pprop) - Log(Pcurr)
                logMR <- sum(func.prop) - sum(func.cur) # log10 of the Metropolis ratio
                
                # We accept if
                # -> logMR > Log10(1) = 0
                # -> Log10(runif(1)) < logMR (accept reject)
                if(logMR>=0 || logMR>log10(runif(1)))
                    {
                        theta.cur   <- theta.prop
                        func.cur    <- func.prop
                        n.accept    <- n.accept + 1
                        rate.accept <- n.accept/n
                    }
            
                if(n>n.burnin)
                    {
                        func.samp[n-n.burnin,1:2] <- func.cur
                        func.samp[n-n.burnin,3:(2+n.theta)] <- theta.cur
                    }

            }

        return(func.samp)
    }

metropolis.plotiter <- function(chain,dim)
    {
        bestparams <- c()
        # 10^(allSamp[,1]+allSamp[,2]) is the unnormalized posterior at each sample
        #thinning
        thinSel  <- seq(from=1, to=nrow(chain), by=200) # thin by factor 100
        postSamp <- chain[thinSel,]

        # Plot MCMC chains 
        par(mfrow=c(4,2), mar=c(3.0,3.5,0.5,0.5), oma=0.5*c(1,1,1,1), mgp=c(1.8,0.6,0), cex=0.9)

        parnames <- c(expression(b[0]), expression(paste(alpha, " / rad")), expression(b[2]), 
                      expression(paste(log, " ", sigma)))

        for(j in 3:dim+3) 
            { # columns of postSamp
                plot(1:nrow(postSamp), postSamp[,j], type="l", xlab="iteration", ylab=parnames[j-2])
                postDen <- density(postSamp[,j], n=2^10)
                bestparams <- c(bestparams,(postDen$x[match(max(postDen$y),postDen$y)]))
                plot(postDen$x, postDen$y, type="l", lwd=1.5, yaxs="i", ylim=1.05*c(0,max(postDen$y)),
                     xlab=parnames[j-2], ylab="density")
                abline(v=bestparams[j-2],col='red')

            }
        return(bestparams)
    }