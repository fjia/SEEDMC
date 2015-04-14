#' Generate Mplus syntax for hypothesized model
#'
#' This is the second step of SEEDMC. 
#' It creates the Mplus syntax for the hypothesized model, 
#' based on the user-specified latent variance-covariance
#' matrix, latent means, and residual variances. 
#' Currently, it only accommodates latent linear or quadratic 
#' growth curve models.
#'
#' @param fl The hypothesized factor loadings matrix 
#' for a latent growth curve model. 
#' Rows indicate indicators/repeated measures, 
#' and columns represent latent intercept and slope(s).
#' @param latcov The hypothesized variance-covariance matrix of  
#' the latent intercept and slope(s).
#' @param latmean The hypothesized mean vector of 
#' the latent intercept and slope(s).
#' @param resvar The residual variances of the indicators. 
#' If \code{resvar} is a scalar, then the 
#' residual variances are the same across time. 
#' It could also be a vector that contains the 
#' residual variances for each of the time points.
#' @param ... Other arguments, currently unused.
#'
#' @export
#'
#' @return Returns an object of class "seedmod", 
#' which contains the following information.
#' \enumerate{
#'   \item Population values of model parameters.
#'   \item Syntax that can be used to define the 
#' 			population model in Mplus.
#' }
#'
#' @references Wu. W., Jia, F., Rhemtulla, M., & Little, T. D. 
#'			(revise and resubmit). Search for efficient complete and 
#'			planned missing data designs for analysis of change. 
#'			Behavioral Research Methods.
#'
#' @examples
#' ### Linear
#' fl <- matrix(NA, 5, 2)  
#' fl[,1] <- 1
#' fl[,2] <- 0:4
#' latcov <- matrix(NA, 2, 2)
#' diag(latcov) <- c(28.776, 8.201)
#' latcov[1, 2] <- 1.56
#' latcov[2, 1] <- 1.56
#' latmean <- c(39.457, 8.063)
#' resvar <- 30
#'
#' modL <- modelMplus (fl = fl, latcov = latcov, 
#'			latmean = latmean, resvar = resvar) 
#' modL 
#'
#' #### Quadratic
#' fl <- matrix(NA, 5, 3)
#' fl[,1] <- 1
#' fl[,2] <- 0:4
#' fl[,3] <- fl[,2]^2
#' latcov <- matrix(NA, 3, 3)
#' diag(latcov) <- c(107.08, 24.60, 1.22)
#' latcov[1, 2] <- -3.69
#' latcov[1, 3] <- -1.36
#' latcov[2, 3] <- -4.96
#' latmean <- c(13.97, -1.15, 0.2)
#' resvar <- c(rep(41.98,4),40)
#'
#' modQ <- modelMplus(fl = fl, latcov = latcov, 
#'		latmean = latmean, resvar = resvar) 
#' modQ

modelMplus <- function (fl = NULL, latcov = NULL, 
						latmean = NULL, resvar = NULL, ...)
{ 

	if (!isSymmetric(latcov)) {
		if (nrow(latcov) == ncol(latcov)) {
			warning("The latcov matrix is not symmetric. Only the upper triangle (including diagonal) is used.")
		} else {
		    stop ("latcov needs to be a square matrix.")
		}
	}
	if (ncol(fl) != nrow(latcov))  stop ("The column number of fl must be equal to the row/column number of latcov.")
	if (ncol(fl)!= length(latmean)) stop ("The column number of fl must be equal to the length of latmean.")
	if (length(resvar)!= 1 && length(resvar)!= nrow(fl)) stop ("The length of resvar must be 1 or equal to the row number of fl.")

	
	iname <- paste0("y", 1: nrow(fl))

	if (ncol(fl)==2){
		fname <- c("i", "s")
	} else {
		fname <- c("i", paste0("s", 1:(ncol(fl) - 1)))
	} 

	rownames(fl) <- iname
	colnames(fl) <- fname
	flscpt <- rep(NA, ncol(fl))
	for (x in (1:ncol(fl))){
		flscpt[x] <-  paste(paste(fname[x], "by"), paste0(iname, "@", fl[, x], collapse = " "), ";")
	}


	rownames(latcov) <- fname
	colnames(latcov) <- fname
	latcovscptM <- matrix(NA, nrow = nrow(latcov), ncol = ncol(latcov))
	for (x in (1:nrow(latcov))) {
		diag(latcovscptM)[x] <- paste(paste0(fname[x], "*", diag(latcov)[x]), ";", collapse = " ")
		if (x == nrow(latcov)) break
		for (y in ((x+1):nrow(latcov))){
			latcovscptM[x, y] <- paste(paste(fname[x],"with",fname[y]), paste0("*", latcov[fname[x],fname[y]]),";")
		}
	}
	latcovscpt <- c(diag(latcovscptM), latcovscptM[upper.tri(latcovscptM)])
	names(latmean) <- fname  
	latscpt <- c(latcovscpt, paste0("[", paste0(fname, "*", latmean, collapse = " "), "] ;"))


	indscpt_1 <- NULL
	indscpt_2 <- NULL
	if (length(resvar)==1){
		resvarV <- rep(resvar, nrow(fl))
		indscpt_1 <- paste0(iname[1], "-", iname[length(iname)], "*", resvar, ";") 
	} else {
		resvarV <- resvar
		indscpt_1 <- paste0(paste0(iname, "*", resvarV), ";") 
	} 

	names(resvarV) <- iname

	resmean <- rep(0, nrow(fl))
	names(resmean) <- iname

	indscpt_2 <- paste0("[", iname[1], "-", iname[length(iname)],"@", resmean[1], "];") 
	indscpt <- c(indscpt_1,indscpt_2)

	out <- list(fl = fl, latcov = latcov, latmean = latmean, 
				resvar = resvarV, resmean = resmean,
				popMplusFL  = flscpt, popMplusLAT = latscpt, popMplusIND = indscpt)
				### call = match.call())
	class(out) <- "seedmod"
	return(out)
}
NULL

#' @export
print.seedmod <-  function(x, ...) 
{   
	out <- list(fl = x$fl, latcov = x$latcov, latmean = x$latmean, 
				resvar = x$resvar, resmean = x$resmean)
	cat("Model Population: \n")
	print(out)
	cat("\n")
	cat("Model Population in Mplus: \n")
	cat("MODEL POPULATION: \n")
	cat(c(x$"popMplusFL", x$"popMplusLAT", x$"popMplusIND"), sep="\n")
}
NULL