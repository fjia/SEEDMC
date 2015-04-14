#' Find the most efficient designs
#'
#' This is the fourth step of SEEDMC. 
#' It summarizes the "seedfoud" object obtained from 
#' the \code{seedmcMplus()} function, and provides 
#' relative efficiency measures of each design for
#' the latent means, variances and covariances.
#'
#' @param x The "seedfound" object obtained from the 
#' \code{seedmcMplus()} function.
#' @param param A vector of the parameters of interests. 
#' Here "I" indicates the intercept, 
#' "S" represents the slope in a linear model,
#' and "S1" and "S2 represent the linear and qudratic sloples
#' in a qudratic model, respectively. 
#' Valid prameter names include:
#' \itemize{
#'   \item "Means.I"
#'   \item "Means.S" for linear models only,
#'			  or "Means.S1" and "Means.S2" 
#'				for quadratic models.
#'   \item "Variances.I"
#'   \item "Variances.S" for linear models only, 
#'			 or "Variances.S1" and "Variances.S2"
#'			 for quadratic models.
#'   \item "I.WITH.S" for linear models only, 
#'		or "I.WITH.S1", "I.WITH.S2" and 
#'		"S2.WITH.S1"for quadratic models.
#' }
#' @param out A character string indicating how the results 
#' should be organized. If "com", the function will compute 
#' relative efficiency only based on the complete data designs. 
#' If "miss", it will only compute relative efficiency for 
#' missing data designs. If "mix", it will compute relative efficiency 
#' based on the entire design pool.
#' @param maxerr A number between 0 and 1. 
#' Only the designs with the proportion of error/warning messages 
#' smaller than \code{maxerr} will 
#' be included when computing the relative efficiency. 
#' The default is 0.05, meaning that the designs, 
#' for which more than 5\% replications have error/warning messages,
#' will be omitted from the computation of relative efficiency.
#' @param recut A number between 0 and 1. 
#' Only the designs whose relative efficiencies 
#' are greater than \code{recut} will be included 
#' in the final result. The default is 0, 
#' meaning that all designs in the pool will be included.
#' @param acrosscut A number between 0 and 1, which 
#' must be equal to or greater than \code{recut}.
#' When a number is specified, an additional table 
#' will be created showing the designs whose relative efficiencies 
#' are greater than \code{recut} across all parameters of interest.
#' @param ... Other arguments, currently unused.
#'
#' @export
#'
#' @return For each parameter of interests, 
#' return the relative efficiency (RE) of 
#' each design, and the extra information.
#'
#' @references Wu. W., Jia, F., Rhemtulla, M., & Little, T. D. 
#'			(revise and resubmit). Search for efficient complete and 
#'			planned missing data designs for analysis of change. 
#'			Behavioral Research Methods.
#'
#' @examples
#' patQ <- designPool(type ="both",  time = 5, 
#'			traj = "quadratic", budget = 100000, unitcost = 20)
#'
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
#' modQ <- modelMplus (fl = fl, latcov = latcov, 
#'		latmean = latmean, resvar = resvar) 
#'
#' seedmcoutQ <- seedmcMplus(pattern = patQ, model = modQ, nreps = 2)
#' #seedmcoutQ$results$Complete
#' seedmcoutQ$results$Complete[[1]]
#' # For presentation purpose, we set nreps = 2. 
#' # In a real study, nrep should be a large number (e.g. 5000).
#'
#' finalQ <- topDesigns(x = seedmcoutQ, 
#'			param = c("Means.I", "Means.S1", "Means.S2"), 
#'			out = "mix", recut = 0, maxerr = 0.4)

topDesigns <- function(x = NULL, param = NULL, out = NULL, 
					maxerr = 0.05, recut = 0, acrosscut = NULL,...) 
{
    ## extrEmpSE() function to extract the empirical 
	## standard errors from the seedmc object
	extrEmpSE <- function(results, label){  
		resChunk <- results 
		##__________________________________________
		numChunk <- rep(NA, length(resChunk))
		for (j in 1:length(resChunk)){
			numChunkLong <- strsplit(gsub("\\D", " ", names(resChunk)[j]), "\\s")
			numChunk[j] <- as.numeric(numChunkLong[[1]][numChunkLong[[1]] != ""][2]) 
		}
		##__________________________________________
		desLabel <- label
		p <- list()
		ps <- list()
		for (i in 1:length(param)){
			p[[i]] <- rep(NA, length(resChunk))
			###for (j in numChunk) {
			for (j in 1:length(resChunk)){
			
				##__________________________________________
				jNew <- which(numChunk == j)
				resDes <- resChunk[[jNew]]$unstandardized
				##__________________________________________
				
				###resDes <- resChunk[[j]]$unstandardized
				rownames(resDes) <- paste(resDes[,1], resDes[,2], sep=".")
				p[[i]][j] <- resDes[param[i], "population_sd"]  
			}
			ps[[i]] <- data.frame(paste0(desLabel, 1:length(resChunk)), param[i], p[[i]]) 
			colnames(ps[[i]]) <- c("Design", "Parameter", "Emp_Std_Error")
		}
		names(ps) <- param
		return(ps)
	}
	
	## effInfo() function compute relative efficiency and 
	## add other information to each element in an extrEmpSE() result
	effInfo <- function(empseout, outindex){  
		extr <- list()
		patN <- list()
		patCode <- list()
		for (i in 1:length(outindex)) {
		
			##__________________________________________
			extrChuck <- x$extraInfo[[outindex[i]]]
			extrChuck <- as.data.frame(extrChuck)
			
			numChunk <- rep(NA, nrow(extrChuck))
			for (j in 1:nrow(extrChuck)){
				numChunkLong <- strsplit(gsub("\\D", " ", rownames(extrChuck)[j]), "\\s")
				numChunk[j] <- as.numeric(numChunkLong[[1]][numChunkLong[[1]] != ""][2]) 
			}
			extrChuck$jNew <- numChunk
			extr[[i]] <- extrChuck[order(extrChuck$jNew),]
			##__________________________________________
			
			###extr[[i]] <- x$extraInfo[[outindex[i]]]
			patN[[i]] <- do.call("rbind", x$patterns[[outindex[i]]]$N)
			patCode[[i]] <- do.call("rbind",x$patterns[[outindex[i]]]$patterns_notation)
		}
		extrM <- do.call("rbind", extr)
		empseout$RE <- 0 
		empseout$Done_Reps <- extrM[,"comprep"]
		empseout$Err_Msgs <- extrM[,"errmsg"]
		empseout$N <- do.call("rbind",patN)
		empseout$Patterns <- do.call("rbind", patCode)
			
		## Compute RE for the designs whose number of 
		## error messages are smaller than rqstrep*maxerr
		nErrTol <- round(extrM[1 ,"rqstrep"]*maxerr, 0)
		good <- empseout[empseout$Err_Msgs <= nErrTol, ] 
		
		good$RE <- min((good$Emp_Std_Error)^2)/(good$Emp_Std_Error)^2  
		
		## Sort based on Emp_Std_Error
		goodSorted <- good[order(good$Emp_Std_Error),]
		
		## Add back the designs whose number of 
		## error messages are greater than rqstrep*maxerr
		failed <- empseout[empseout$Err_Msgs > nErrTol, ] 
		
		allSorted  <- rbind(goodSorted, failed)
		
		## Keep designs with RE > recut
		allSortedKeep <- allSorted[allSorted$RE >= recut,]
		
		return(allSortedKeep)
	}
	
	comSE <- NULL
	missSE <- NULL

	if (out == "com") {
		if (x$patterns$Type == "miss") { 
			stop ("There is no complete data design results in the seedmc object.")
		} else {
		comSE <- extrEmpSE(results = x$results$Complete, label = "C")
		Effs <- lapply(comSE, effInfo, outindex = "Complete")
		}
	} else if (out == "miss") {
		if (x$patterns$Type == "com") { 
			stop ("There is no missing data design results in the seedmc object.")
		} else {
		missSE <- extrEmpSE(results = x$results$Missing, label = "M")
		Effs <- lapply(missSE, effInfo, outindex = "Missing")
		}
	} else if (out == "mix"){
	    comSE <- extrEmpSE(results = x$results$Complete, label = "C")
	    missSE <- extrEmpSE(results = x$results$Missing, label = "M")
		mixSE <- mapply(rbind, comSE, missSE, SIMPLIFY = FALSE, USE.NAMES = TRUE)
		Effs <- lapply(mixSE, effInfo, outindex = c("Complete","Missing"))
    }

	
	if (length(acrosscut) != 0) {
		if (acrosscut < 0L && acrosscut > 1L) stop ("acrosscut must be a value between 0 and 1.")
		if (acrosscut < recut) stop ("acrosscut must be greater than or equal to recut.")
		numEffs <- list()
		someEffs <- list()
		for (i in 1:length(Effs)) {
			temp <- Effs[[i]]
			numEffs[[i]] <- which(temp$RE >= acrosscut)
			someEffs[[i]] <- temp$Design[numEffs[[i]]]
		}
		pos <- unlist(someEffs)
		postab <- table(pos)
		intsec <- names(postab)[postab==length(someEffs)]
		temp2 <- Effs[[1]]
		mat <- match(intsec, temp2$Design)
		if (length(mat) == 0){
			acEffs <- data.frame(Design = NA, Parameters = paste0(param, collapse = "_"), REcut = acrosscut)	
		} else {
			acEffs <- cbind(Design = temp2[mat, "Design"], 
						Parameters = paste0(param, collapse = "_"), REcut = acrosscut, 
						temp2[mat, c("Done_Reps", "Err_Msgs", "N", "Patterns")])
		}
	
		out <- list(Single_Parameter = Effs, Across_Parameters = acEffs)
	} else {
		out <- Effs
	}


	return(out)  
}
NULL
