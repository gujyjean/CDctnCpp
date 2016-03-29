#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' @useDynLib CDctnCpp
#' @importFrom Rcpp sourceCpp
NULL

# Generate the distribution of in-degrees given the maximum in-degree and the total number of edges in a graph.
# maxdeg: maximum in-degree; node: the number of nodes; nedge: the total number of edges
degreeG <- function(maxdeg, node, nedge) {
  degree.dist <- rep(0,maxdeg+1)
  nx <- node; nex <- nedge
  if(maxdeg>=2) {
    for (i in maxdeg:2) {
      lower <- max(0, ceiling(nex-(nx-i/2)*(i-1)))
      upper <- min(nx-i, floor(nex/i))
      if(lower<upper)    {degree.dist[i+1] <- sample(lower:upper,1)}
      else if(lower==upper)    {degree.dist[i+1] <- lower}
      else  {next}
      nx <- nx-degree.dist[i+1]
      nex <- nex-i*degree.dist[i+1]
      if(nx<=0 || nex<0)    stop("Error")
    }
  }
  if (nx>=(nex+1))   degree.dist[2] <- nex
  else  stop("Error!")
  degree.dist[1] <- node-sum(degree.dist[2:(maxdeg+1)])
  return(degree.dist)
}



# Generate DAGs given the maximum in-degree and the total number of edges.

#' Generate DAGs given the maximum in-degree and the total number of edges
#' @param maxdeg maximum in-degree
#' @param node the number of nodes
#' @param nedge the total number of edges
#' @return A list of 4.
#' \item{ordex}{parent set matrix, of dimension maxdeg*node.
#'              Each nonzero entry in column j lists one parent of node j.
#'              The parents	of node j are sorted in ascending order.}
#' \item{trueG}{the adjacency matrix of the generated DAG}
#' \item{ts}{one topological sort of the generated DAG}
#' \item{degree.dist}{the distribution of in-degrees of all nodes.
#'                    Element j lists the number of nodes whose in-degree equals j-1}
#' @export
GGen <- function(maxdeg, node, nedge) {

  ordex <- matrix(0, maxdeg, node)		# matrix listing the parents of each node (similar to the adjacency list)
  ts <- sample(1:node)					# topological sort
  degree.dist <- degreeG(maxdeg, node, nedge)
  indegree <- rep(0,node)					# indegree is matched to ts
  ind <- rep(TRUE, node)
  for (i in (maxdeg+1):1) {
    if(degree.dist[i]>0) {
      candidate <- (i:node)[ind[i:node]]
      if(length(candidate)==1)   {
        indegree[candidate] <- (i-1)
        ind[candidate] <- FALSE
      }
      else  {
        pos <- sample(candidate,degree.dist[i])
        indegree[pos] <- (i-1)
        ind[pos] <- FALSE
      }
    }
  }
  for (i in 1:node) {
    np <- indegree[i]
    if(np==0)   next
    else if(i==2)    ordex[1,ts[i]] <- ts[1]
    else   {
      pa <- sample(ts[1:(i-1)], np)
      ordex[1:np,ts[i]] <- sort(pa)
    }
  }
  trueG <- matrix(0,node,node)
  for (i in 1:node) {
    trueG[ordex[ordex[,i]!=0,i],i] <- 1
  }

  return(list(ordex=ordex, trueG=trueG, ts=ts, degree.dist=degree.dist))

}



# Generate Gaussian data from DAGs using the outputs from the function GGen as the inputs.
# coeff: the beta coefficient; noi: the number of interventions/data blocks; nobs: the number of interventional samples in each data block

#' Generate Gaussian interventional data from DAGs
#' @param ordex parent set matrix, of dimension maxdeg*node. See GGen
#' @param ts one topological sort of the DAG. See GGen
#' @param coeff the beta value used to generate Gaussian samplesb$ob
#' @param noi the number of nodes that have interventional smaples
#' (i.e., the number of data blocks). Each data block contains
#' interventional samples where one particular node is experimentally manipulated
#' @param nobs the sample size of each data block
#' @return A list of 3
#' \item{data}{data matrix of dimension (nobs*noi)*node}
#' \item{ivn}{the vector of nodes that have interventional samples}
#' \item{obslist}{list of indices of observational samples for each node
#'                (excluding interventional samples).
#'                The list has length equal to node}
#' @export
DatGen <- function(ordex, ts, coeff, noi, nobs) {

  maxdeg <- nrow(ordex); node <- ncol(ordex)
  beta <- matrix(coeff, maxdeg, node)
  ivn <- sort(sample(1:node, noi))					# nodes that have interventional data
  fX <- matrix(rnorm(noi*nobs,0,1),ncol=noi)
  data <- datGenC(as.integer(node),as.integer(maxdeg),as.integer(ordex),as.integer(ts),
                as.integer(noi),as.integer(nobs),as.integer(ivn),as.double(fX),as.double(beta))
  data <- matrix(data, ncol=node, byrow=TRUE)
  obslist <- list()									# list of indices of observational samples for each node (excluding interventional samples)
  kdex <- 0
  for (k in 1:node) {
    if (any((k-ivn)==0)) {
      kdex <- kdex+1
      obslist[[k]] <- (1:nrow(data))[-((nobs*kdex-nobs+1):(nobs*kdex))]
    } else {
      obslist[[k]] <- 1:nrow(data)
    }
  }
  return(list(data=data, ivn=ivn, obslist=obslist))

}

# Compute solutions to the penalized likelihood problem on a grid of lambdas.

#' compute solutions to the penalized likelihood estimation
#' problem on a grid of lambdas
#' @param node the number of nodes of the DAG
#' @param data the data matrix of dimension (nobs*noi)*node. See DatGen
#' @param obslist list of observational samples. See DatGen
#' @param eorder matrix indicating the order to cycle through
#'        the pair of nodes in the CD algorithm. nrow of the matrix is
#'        0.5*node*(node-1), and ncol of the matrix is 2.
#'        Each row is a pair of nodes
#' @param fmlam the ratio of lambda_min/lambda_max
#' @param nlam the number of lambda values
#' @param eps threshold of convergence for the CD algorithm
#' @param beta the initial coefficient matrix B^0.
#'        If missing, will assume to be all zero matrix
#' @param lambda the sequence of lambda values to be used.
#'        If missing, it will be constructed from data using fmlam and nlam.
#' @param gamma the exponent in the weights for the adaptive lasso penalty
#' @param lowerbound the constant used to truncate the weights
#'        for the adaptive lasso penalty
#' @return A list of 2
#' \item{coef}{matrix of estimated coefficients, of dimension nlam*(node^2).
#'             Each row contains the estimate of the coefficient
#'             matrix (stored as a vector) at a single value of lambda}
#' \item{lambda}{the sequence of lambda values used}
#' @export
cdpathC <- function(node, data, obslist, eorder, fmlam=0.001, nlam=50, eps=1e-4, beta, lambda, gamma=0.15, lowerbound=1e-4)  {

  argind1 <- missing(lambda); argind2 <- missing(beta)

  obs.leng <- sapply(obslist, length)
  ncoef <- node*node
  sigma <- rep(0,node)
  xlist <- list()
  ylist <- list()
  Normx <- NULL
  inverseWeights <- matrix(0,node,node)
  scale.factor <- NULL
  for (i in 1:node) {
    xtemp <- scale(data[obslist[[i]], ], TRUE, FALSE)
    normx <- sqrt(colSums(xtemp^2))
    xtemp <- scale(xtemp, FALSE, normx)
    scale.factor <- c(scale.factor,normx[i]/normx)

    inverseWeights[-i,i] <- W <- pmax(abs(lm(xtemp[,i]~xtemp[,-i]-1)$coef),lowerbound)^gamma		# compute weights for the adaptive lasso
    xtemp[,-i] <- xtemp[,-i]*rep(W,rep(obs.leng[i],node-1))            								# rescale x for the adaptive lasso

    Normx <- c(Normx, colSums(xtemp^2))
    ylist[[i]] <- xtemp[,i]
    sigma[i] <- sum(xtemp[,i]^2)/obs.leng[i]
    xtemp[, i] <- 0
    xlist[[i]] <- xtemp
  }
  re.prod <- prodG(as.integer(node),xlist,ylist,as.integer(obs.leng))
  inprod <- re.prod$inprod                            # inner products of each feature with y initially
  feature.prod <- re.prod$fprod                     	# inner products between features

  if(argind1) {
    maxtmp1 <- max(abs(inprod)/rep(sigma,rep(node,node)))
    max.lam <- 7*maxtmp1
    min.lam <- fmlam*max.lam
    lambda <- exp(seq(log(max.lam), log(min.lam), log(fmlam)/(nlam-1)))
  }
  if(argind2)    beta <- matrix(0, node, node)


  norm <- penalty <- rep(0,node)
  for (i in 1:node) {
    residual <- ylist[[i]]-xlist[[i]] %*% beta[,i]
    norm[i] <- sum(residual^2)
    penalty[i] <- lambda[1]*sum(abs(beta[,i]))
  }
  rsq <- sum(obs.leng*log(sigma)/2)+sum(norm/(2*sigma))+sum(penalty)

  eorleng <- nrow(eorder)
  eorder <- t(eorder)
  result <- gridCD(as.integer(node), as.double(beta), as.double(lambda), as.double(inprod), as.double(feature.prod),
               as.integer(eorder), as.double(rsq), as.double(eps), coef=double(nlam*ncoef), RSQ=double(nlam), ITER=integer(nlam), as.integer(nlam-1),
               as.double(sigma), sig=double(nlam*node), as.double(norm), as.integer(obs.leng), as.integer(eorleng))
  coef <- matrix(result*rep(inverseWeights,nlam), ncol=ncoef, byrow=TRUE)




  # use new weights to recompute the solutions
  inverseWeights <- matrix(abs(coef[nlam,])^gamma,ncol=node)
  for (i in 1:node) {
    xtemp <- scale(data[obslist[[i]], ], TRUE, FALSE)
    normx <- sqrt(colSums(xtemp^2))
    xtemp <- scale(xtemp, FALSE, normx)

    xtemp[,-i] <- xtemp[,-i]*rep(inverseWeights[-i,i],rep(obs.leng[i],node-1))

    ylist[[i]] <- xtemp[,i]
    sigma[i] <- sum(xtemp[,i]^2)/obs.leng[i]
    xtemp[, i] <- 0
    xlist[[i]] <- xtemp
  }
  re.prod <- prodG(as.integer(node),xlist,ylist,as.integer(obs.leng))
  inprod <- re.prod$inprod
  feature.prod <- re.prod$fprod

  if(argind1) {
    maxtmp1 <- max(abs(inprod)/rep(sigma,rep(node,node)))
    max.lam <- 7*maxtmp1
    min.lam <- fmlam*max.lam
    lambda <- exp(seq(log(max.lam), log(min.lam), log(fmlam)/(nlam-1)))
  }
  if(argind2)    beta <- matrix(0, node, node)

  norm <- penalty <- rep(0,node)
  for (i in 1:node) {
    residual <- ylist[[i]]-xlist[[i]] %*% beta[,i]
    norm[i] <- sum(residual^2)
    penalty[i] <- lambda[1]*sum(abs(beta[,i]))
  }
  rsq <- sum(obs.leng*log(sigma)/2)+sum(norm/(2*sigma))+sum(penalty)

  result <- gridCD(as.integer(node), as.double(beta), as.double(lambda), as.double(inprod), as.double(feature.prod),
               as.integer(eorder), as.double(rsq), as.double(eps), coef=double(nlam*ncoef), RSQ=double(nlam), ITER=integer(nlam), as.integer(nlam-1),
               as.double(sigma), sig=double(nlam*node), as.double(norm), as.integer(obs.leng), as.integer(eorleng))
  coef <- matrix(result*rep(inverseWeights,nlam), ncol=ncoef, byrow=TRUE)




  coef <- coef*rep(scale.factor,rep(nlam,ncoef))
  return(list(coef=coef, lambda=lambda))

}



# Compute the negative log-likelihood using unpenalized coefficients.

#' compute the negative log-likelihood using unpenalized coefficients
#' @param node the number of nodes of the DAG
#' @param G the adjacency matrix of the DAG
#' @param data the data matrix of dimension (nobs*noi)*node. See DatGen
#' @param obslist list of observational samples. See DatGen
#' @return pLog, the negative log-likelihood
#' @export
plog <- function(node, G, data, obslist)  {

  obs.leng <- sapply(obslist, length)
  RSS <- rep(0,node)
  for (i in 1:node) {
    pa <- which(G[,i]!=0)
    xtemp <- scale(data[obslist[[i]], ], TRUE, FALSE)
    if(length(pa)>0) {
      RSS[i] <- sum(lm(xtemp[,i]~xtemp[,pa]-1)$res^2)
    } else {
      RSS[i] <- sum(xtemp[,i]^2)
    }
  }
  pLog <- sum(obs.leng*log(RSS/obs.leng)/2)+sum(obs.leng/2)
  return(pLog)

}
