library(dplyr)

K2 <- function(dataset, u)
    {
        n.nodes <- num.variables(dataset)# G
        curr.g <- matrix(0L,n.nodes,n.nodes) # integers!  G
        data <- as.data.frame(raw.data(dataset))

        for(i in 2:n.nodes)
            {

                # How many parents does the i-th node have? Assuming there is at most one arc between any two nodes and no cycles
                n.parents <- 0#sum(curr.g[,i])# At this stage, this is 0
                pred <- 1:(i-1) # Preceding nodes
                left <- pred
                parents <- c()
                lp.old <- f(data,i,parents)
                ok.to.proceed <- TRUE

                while(ok.to.proceed & (n.parents < u) & length(left) > 0)
                    {
                        g <- function(n.p){f(data,i,c(parents,n.p))}
                        proposal <- lapply(left,g)
                        j <- which.max(proposal)
                        lp.new <- proposal[[j]]
                        if(lp.new > lp.old){
                            curr.g[j,i] = 1L # The j-th node has been chosen as parent for the i-th one.
                            parents <- c(parents,left[j]) # Adding the chosen node as parent of the i-th node
                            n.parents <- n.parents + 1 #sum(curr.g[,i])
                            left <- left[-c(j)] # Removing from options the node that has been added as parent of the i-th node.
                            lp.old <- lp.new
                        }
                        else{ok.to.proceed <- FALSE}
                    }
                }
        return(curr.g)
    }


############ AUXILIARY FUNCTIONS ###############

prob.noparents <- function(data,namecol,lprod)
    {
        #lprod <- log(prod)
        col <- dplyr::pull(data, namecol)
        nunique <- length(unique(col))
        lprod <- lprod+lfactorial(nunique-1)
        den <- nunique - 1
    
        for(i in 1:nunique)
            {
                lprod <- lprod+lfactorial(length(col[col == unique(col)[i]]))
                den <- den + length(col[col == unique(col)[i]])
            }
        lprod <- lprod - lfactorial(den)

        #nprod <- exp(lprod)
    
        return(lprod)
    }

is.eq <- function(row1,row2){return(row1 == row2)}

prob.parents <- function(data,namecol,lprod,parents)
    {
        #lprod <- log(prod)
        col <- dplyr::pull(data, namecol)
        n.parents <- length(parents)
        col.parents <- data[parents]
        r <- length(unique(col))

        q <- 1
        combined <- list()
        for(j in 1:length(parents))
            {
                q <- q*length(unique(col.parents[,j]))
                combined[[j]] <- unique(col.parents[,j])    
            }
        combinations <- do.call(expand.grid, combined)
        # for j in 1:qi
        for(j in 1:q)
            {
                w  <- combinations[j,]
                # Compute Nijk!
                nij <- 0
                for(k in 1:r)
                    {
                        wij <- c(w,unique(col)[k])
                        nijk <- sum(apply(apply(cbind(col.parents, col),1,is.eq,row2=wij),2,all))
                        nij <- nij + nijk
                        lprod <- lprod+lfactorial(nijk)
                    }
                lprod <- lprod+lfactorial(r - 1) - lfactorial(nij + r - 1)
            }

        #nprod <- exp(lprod)
    
        return(lprod)

    }

#############################################

prob.model <- function(BN,D)
    {
        nvar <- length(nodes(BN))
        prod <- 1
        for(i in 1:nvar)
            {
                if(length(parents(BN, nodes(BN)[i])) == 0)
                    {prod <- prob.noparents(D,nodes(BN)[i],prod)}
                else
                    {prod <- prob.parents(BN,D,nodes(BN)[i],prod)}
                
            }
    
        return(prod)
    }

f <- function(data, i, parents)
            {
                #prod <- 1
                lprod <- 0
                colname <- colnames(data)[i]
                if(length(parents) == 0)
                    {lprod <- prob.noparents(data,colname,lprod)}
                else
                    {lprod <- prob.parents(data,colname,lprod,parents)}
               
                return(lprod)
            }