#' Generate design pool
#'
#' This is the first step of SEEDMC. It generates a design pool. 
#' This step returns the response patterns of the possible designs, 
#' proportion of each pattern, and sample size for each design.
#'
#' @param type A character string indicating the design type. 
#' Options include "com" (for complete data designs only), 
#' "miss" (for planned missing data designs only), and "both".
#' @param traj A character string defining the shape of trajectory. 
#' Options include "linear" and "quadratic". 
#' @param time The number of maximum possible measurement occasions. 
#' If \code{type} = "miss" or "both", 
#' the maximum number of \code{time} is 6; 
#' if \code{type} = "com", 
#' the maximum number of \code{time} is 10.
#' If \code{traj} = "linear", 
#' the minimum number of \code{time} is 3; 
#' if \code{traj} = "quadratic", 
#' the minimum number of \code{time} is 4.
#' @param budget The total amount of budget for the longitudinal study. 
#' It is used together with \code{unitcost} to determine the 
#' sample size (N) for each design.
#' @param unitcost The cost of a single measurement for each participant.
#' It is used together with \code{budget} to determine the 
#' sample size (N) for each design. 
#' @param costratio The cost ratio, indicating how many times more 
#' it costs by recruiting a new participants than 
#' a repeated measure (\code{unitcost}). The default is 1.
#' Given \code{budget}, \code{unitcost} and \code{costratio}, 
#' the sample size (N) can be computed as 
#' \code{budget/(costratio * unitcost + (time - 1) * unitcost)}.
#' @param nfix The fixed sample size. 
#' Used only when the sample size is fixed across different designs,
#' DO NOT use \code{nfix} together with
#' \code{budget}, \code{unitcost} and \code{costratio}, 
#' the sample size determined by \code{budget}, \code{unitcost}
#' and \code{costratio}, will overwrite \code{nfix}. The default is NULL. 
#' @param attrition A vector of attrition rate at all possible
#' measurement occasions. The default is NULL. 
#' @param ... Other arguments, currently unused.
#'
#' @export
#'
#' @return Returns an object of class "seedpool", 
#'	which describes the design pool. 
#'	It contains the following information for 
#'	each design: 
#' \enumerate{
#'   \item Response patterns (0 = complete, 1 = missing);
#'   \item Attrition rates;
#'   \item Proportions of each pattern;
#'   \item Sample sizes.
#' }
#'
#' @references Wu. W., Jia, F., Rhemtulla, M., & Little, T. D. 
#'			(revise and resubmit). Search for efficient complete and 
#'			planned missing data designs for analysis of change. 
#'			Behavioral Research Methods.
#'
#' @examples
#' patC <- designPool(type ="com", time = 5, 
#'			traj = "linear", budget = 100000, 
#'			unitcost = 20, nfix = NULL) 
#' summary(patC)
#' patC$Complete$patterns
#' patC$Complete$patterns_notation
#'
#' patM <- designPool(type ="miss", time = 6, 
#'			traj = "linear", budget = 100000, 
#'			unitcost = 20, attrition = rep(0.5, 6))
#' summary(patM)
#' patM$Missing$patterns
#' patM$Missing$patterns_notation
#'
#' patQ <- designPool(type ="both",  time = 5, 
#'			traj = "quadratic", budget = 100000, unitcost = 20)
#' summary(patQ) 

designPool <- function (type = NULL, traj = NULL, time = NULL, budget = NULL, 
						unitcost = NULL, costratio = 1L, nfix = NULL, attrition = NULL,...)
{ 
	
	if ((time %% 1) != 0) stop ("time needs to an integer.")
	
	if (traj == "linear" && time < 3) stop ("time needs to equal to or greater than 3 for a linear trajectory.")
	if (traj == "quadratic" && time < 4) stop ("time needs to be equal to or greater than 4 for a quadratic trajectory.")
	
	if (type == "com" && time > 10) stop ("time needs to be smaller or equal to 10 with type = 'com'.")
	if (type == "miss" && time > 6) stop ("time needs to be smaller or equal to 6 with type = 'miss'.")
	if (type == "both" && time > 6) stop ("time needs to be smaller or equal to 6 with type = 'both'.")
	

	if (length(attrition) != 0 && length(attrition) != time) {
		stop ("attrition must be NULL or a numeric vector with T elements (T = the number of measurement occassions, i.e., time).")
	} else if (length(attrition) == time && !is.numeric(attrition)) {
		stop ("attrition must be NULL or a numeric vector with all elements between 0 and 1.")
	} else if (length(attrition) == time && (sum(attrition >= 1) >= 1 || sum(attrition < 0) >= 1)) {
		stop ("attrition must be NULL or a numeric vector with all elements between 0 and 1.")
	}
	
	#######
	
	T <- as.integer(time)
	Budget <- as.numeric(budget)
	C2 <- as.numeric(unitcost) 
	CR <- as.numeric(costratio)
	outCom <- NULL
	outMiss <- NULL	

	if (type == "com" || type == "both" ) {	

		## Tc: number of waves that have complete data
		## Tm: number of waves that have missing data
		## Design: number of possible designs
		if (traj == "linear") {
			leastT <- 3
		} else if (traj == "quadratic") {
			leastT <- 4
		}

		Tc <- leastT:T
		Tm <- T - Tc
		Design <- list()    
		
		temp1 <- list()
		temp2 <- list()
		for (i in 1:length(Tc)){
			if (Tc[i] < T) {
				Design[[i]] <- ncol(combn(2:(T - 1), (Tc[i]-2)))
			} else { 
				Design[[i]] <- 1 
			}
			temp1[[i]] <- expand.grid(Tc[i], Tm[[i]])
			temp2[[i]] <- cbind(T, temp1[[i]], Design[[i]])
		}
		cond1 <- do.call("rbind",temp2)
		colnames(cond1) <- c("Time", "Tc", "Tm", "Design")
		
		
		##Pattern: number of missing patterns for each design
		Pattern <- rep(NA, nrow(cond1)) 
		for (i in 1:nrow(cond1)){
			Pattern[i] <- 1
		}
		
		## j: the lable of a design in each condition
		## jj: the lable of a design in all conditions
		if (length(Budget) == 0) Budget <- NA 
		if (length(C2) == 0) C2 <- NA	

		cond2 <- cbind(cond1, Pattern, Budget, C2, CR)
		temp3 <- list()
		for (i in 1:nrow(cond2)){
			temp3[[i]] <- cond2[rep(i, each=cond2[i,"Design"]),]
			j <- 1:nrow(temp3[[i]])
			temp3[[i]] <- cbind(temp3[[i]], j)
		}
		cond3 <- do.call("rbind",temp3)
		## create a string vector: "1", "2",..., "7"
		rownames(cond3) <- sprintf("%.0f", 1:nrow(cond3))  
		jj <- 1:nrow(cond3)
		allDesigns <- cbind(cond3, jj)
		## split data frame to list by row 
		allDesignList <- split(allDesigns, 1:NROW(allDesigns)) 
		################################################################################################
		
		
		## patternCom() function
		patternCom <- function(xdesign)
		{
			aTime <- as.integer(xdesign["Time"])
			aTc <- as.integer(xdesign["Tc"])
			aBudget <- as.numeric(xdesign["Budget"])
			aC2 <- as.integer(xdesign["C2"])
			aCR <- as.integer(xdesign["CR"])
			aj <- as.integer(xdesign["j"])
			
			## 1. PMD patterns for each design: 0 = complete, 1 = missing
			PMDPat <- list()
			if (aTime == aTc) {
				PMDPat <- rep(0, aTime)  
			} else {
				PMDPat <- c(0, rep(1, aTime-2), 0)
				comWs <- combn(2:(aTime - 1), (aTc-2))
				comW <- comWs[, aj]
				PMDPat[comW]  <- 0	
			}
			
			
			## 2. PMD patterns for each design using short notation
			whichcom <- toString(which(PMDPat == 0))
			patCode <- paste0("{", whichcom, "}")
			
			
			## 3. Probs of each pattern
			patProbs <- 1
			
			
			## 4. Sample size for each design
			if (!is.na(aBudget) && !is.na(aC2)) {
				aC1 <- aC2*aCR
				N <- floor(aBudget/(aC1 + (aTc - 1)*aC2)) 
				if (N < 1) {
				stop ("make sure that budget/(costratio * unitcost + (time - 1) * unitcost) >= 1.")
				}
				
			} else if (is.na(aBudget) && is.na(aC2) && length(nfix) != 0)  {
				N <- nfix
				if (N < 1) {
					stop ("make sure that nfix >= 1.")
				}
				
			} else {
				stop ("budget and unitcost need to be both specified (to determine varing sample sizes) 
						or both unspecified (to use a fixed sample size).")
			} 
			return(list(patterns = PMDPat, patterns_notation = patCode, probs = patProbs, N = N))
		}
		
		outCom <- lapply(allDesignList,  patternCom)
		
	} 
	
	if (type =="miss" || type == "both") {	

		## Tc: number of waves at which all participants have complete data
		## Tm: number of waves at which each participant has missing data 
		## Design: number of possible designs
		Tc <- 0:(T - 2)  
		Tm <- list()   
		Design <- list()  
			
		temp1 <- list()
		temp2 <- list()
		for (i in 1:length(Tc)){
			Tm[[i]] <- 1:(T - Tc[i] -1)  
			if (Tc[i]==0) {Design[[i]] <- 1
			} else {
				Design[[i]] <- ncol(combn(T, Tc[i]))
			}
			temp1[[i]] <- expand.grid(Tc[i], Tm[[i]])
			temp2[[i]] <- cbind(T, temp1[[i]], Design[[i]])
			}
		cond1 <- do.call("rbind",temp2)
		colnames(cond1) <- c("Time", "Tc", "Tm", "Design")
		## get ride of the design which has no overlap one anny time point
		cond1 <- cond1[(cond1$Time - cond1$Tm) > 1, ]  
		
		
		## Pattern: number of missing patterns for each design
		Pattern <- rep(NA, nrow(cond1))   
		for (i in 1:nrow(cond1)){
			Pattern[i] <- ncol(combn((cond1[i,"Time"] - cond1[i,"Tc"]),  cond1[i,"Tm"]))
		}
		
		## j: the lable of a design in each condition
		## jj: the lable of a design in all conditions
		if (length(Budget) == 0) Budget <- NA 
		if (length(C2) == 0) C2 <- NA	
		
		cond2 <- cbind(cond1, Pattern, Budget, C2, CR)
		temp3 <- list()
		for (i in 1:nrow(cond2)){
			temp3[[i]] <- cond2[rep(i, each=cond2[i,"Design"]),]
			j <- 1:nrow(temp3[[i]])
			temp3[[i]] <- cbind(temp3[[i]], j)
		}
		cond3 <- do.call("rbind",temp3)
		## create a string vector: "1", "2",..., "49"
		rownames(cond3) <- sprintf("%.0f", 1:nrow(cond3))  
		jj <- 1:nrow(cond3)
		allDesigns <- cbind(cond3, jj)
		## split data frame to list by row 
		allDesignList <- split(allDesigns, 1:NROW(allDesigns)) 
		################################################################################################

		
		### patternMiss() function
		patternMiss <- function(xdesign)
		{		
			aTime <- as.integer(xdesign["Time"])
			aTc <- as.integer(xdesign["Tc"])
			aTm <- as.integer(xdesign["Tm"])
			aPattern <- as.integer(xdesign["Pattern"])
			aBudget <- as.numeric(xdesign["Budget"])
			aC2 <- as.integer(xdesign["C2"])
			aCR <- as.integer(xdesign["CR"])
			aj <- as.integer(xdesign["j"])	
			
			## 1. PMD patterns for each design: 0 = complete, 1 = missing
			## Starting from all missing
			mat <- matrix(1, aPattern, aTime) 
			
			## Determining the location of complete wave for everyone
			if (aTc==0) {
				matC <- mat
			} else {
				matC <- mat
				comWs <- combn(aTime, aTc)
				comW <- comWs[, aj]
				matC[,comW]  <- 0	 ## 0 = complete
			}
			
			## Allocating missing cells
			matP <- matC 
			allTime <- 1:aTime
			if (aTc==0) {
				comPs <- combn(allTime, (aTime - aTc - aTm))
			} else {
				comPs <- combn(allTime[-comW], (aTime - aTc - aTm))
			}
			for (i in 1:nrow(matP)){
				comP <- comPs[, i]
				matP[i, comP] <- 0 ## 0 = complete
			}
			
			## 2. PMD patterns for each design using short notation
			comWs <- combn(aTime, aTc)
			comW <- comWs[, aj]
			whichcom <- toString(comW)

			whichcomSpl <- unlist(strsplit(whichcom, ","))
			comOc <- sub(" ","", whichcomSpl)
			wholeOc <- sprintf("%.0f", 1:aTime)
			missOc <- wholeOc[!(wholeOc %in% comOc)]
			
			aT_Tm <- aTime - aTm 			
			missComb <- combn(missOc, (aT_Tm - aTc))
			pc <- rep(NA, ncol(missComb))
			for (j in 1:ncol(missComb)){
				pc0 <- c(comOc, missComb[, j])
				pc0 <- sort(as.numeric(pc0[pc0!="NA"]))
				pc[j] <- toString(pc0)
			}
			patCode <- paste0("{", paste0(pc, collapse = "| "), "}")
			
			
			## 3. Probs of each pattern
			patProbs <- c(rep(round(1/aPattern,3), (aPattern - 1)), 
						1 - round(1/aPattern,3)* (aPattern - 1))
						
			## 4. Sample size for each design
			if (!is.na(aBudget) && !is.na(aC2)) {
				aC1 <- aC2*aCR
				N <- floor(aBudget/(aC1 + (aTime - aTm - 1)*aC2))
				if (N < 1) {
					stop ("make sure that budget/(costratio * unitcost + (time - 1) * unitcost) >= 1.")
				}
			} else if (is.na(aBudget) && is.na(aC2) && length(nfix) != 0)  {
				N <- nfix
				if (N < 1) {
					stop ("make sure that nfix >= 1.")
				}
			
			} else {
				stop ("budget and unitcost need to be both specified (to determine varing sample sizes) 
						or both unspecified (to use a fixed sample size).")
			} 
			
			return(list(patterns = matP, patterns_notation = patCode, probs = patProbs, N = N))
		}	
		
		outMiss <- lapply(allDesignList,  patternMiss)
	}	

	
	## reOrgList() function
	reOrgList <- function(multilist, element)
	{
		out <- multilist[[element]]
	}
	
	Complete_patterns <- lapply(outCom, reOrgList, element = "patterns")
	Complete_patterns_notation <- lapply(outCom, reOrgList, element = "patterns_notation")
	Complete_probs <- lapply(outCom, reOrgList, element = "probs")
	Complete_N <- lapply(outCom, reOrgList, element = "N")
	
	
	Missing_patterns <- lapply(outMiss, reOrgList, element = "patterns")
	Missing_patterns_notation <- lapply(outMiss, reOrgList, element = "patterns_notation")
	Missing_probs <- lapply(outMiss, reOrgList, element = "probs")
	Missing_N <- lapply(outMiss, reOrgList, element = "N")
	
	out <- list(Type = type, Time = time, Traj = traj, Attrition = attrition, 
				Complete = list(patterns = Complete_patterns, patterns_notation = Complete_patterns_notation,
								probs = Complete_probs, N = Complete_N, Attrition = attrition),
				Missing = list(patterns = Missing_patterns, patterns_notation = Missing_patterns_notation,
								probs = Missing_probs, N = Missing_N, Attrition = attrition),
				call = match.call())
	class(out) <- "seedpool"	
	return(out)
}
NULL

#' @export
print.seedpool <-  function(x, ...) 
{   
	out <- list(Type = x$Type, Traj = x$Traj, Time = x$Time,
				Attrition = x$Attrition)
				
	cat("Design Pool: \n")
	print(out)
	cat("\n")
	if (x$Type=="com" || x$Type=="both" ){
		cat("To check patterns of complete data designs, call element $Complete$patterns or $Complete$patterns_notation. \n")
		cat("To check sample sizes of complete data designs, call element $Complete$N. \n")
	}
	
	cat("\n")
	if (x$Type=="miss" || x$Type=="both" ){
		cat("To check patterns of missing data designs, call element $Missing$patterns or $Missing$patterns_notation. \n")
		cat("To check sample sizes of missing data designs, call element $Missing$N. \n")
	}
	
}
NULL

#' @export
summary.seedpool <- function(object,...) 
{
	ndesignCom <- NULL
	ndesignMiss <- NULL
	nrepeatCom <- NULL
	nrepeatMiss <- NULL
	
	if (object$Type=="com" || object$Type=="both" ){
		comPat <- object$Complete$patterns
		ndesignCom <- length(comPat)
		for (i in 1:ndesignCom) {
			nrepeatCom[i] <- sum(comPat[[i]] == 0) 
		}
	}
	if (object$Type=="miss" || object$Type=="both" ){
		missPat <- object$Missing$patterns
		ndesignMiss <- length(missPat)
		for (i in 1:ndesignMiss) {
			nrepeatMiss[i] <- rowSums(missPat[[i]] == 0)[1]
		}
	}
	res <- list(call = object$call, type = object$Type, 
				ndesignCom = ndesignCom, ndesignMiss = ndesignMiss,
				nrepeatCom = nrepeatCom, nrepeatMiss = nrepeatMiss)
	class(res) <- "summary.seedpool"
	return(res)

}
NULL

#' @export
print.summary.seedpool <- function(x,...) 
{
	tblCom <- NULL
	tblMiss <- NULL
	if (length(x$nrepeatCom)!=0){
		tblCom <- table(x$nrepeatCom)
		tblComSum <- as.data.frame(addmargins(tblCom))
		colnames(tblComSum) <- c("repeated_measures", "num_of_designs")
	}
	
	if (length(x$nrepeatMiss)!=0){
		tblMiss <- table(x$nrepeatMiss)
		tblMissSum <- as.data.frame(addmargins(tblMiss))
		colnames(tblMissSum) <- c("repeated_measures", "num_of_designs")
	}
	
	cat("Call: \n")
	print(x$call)
	cat("\n")
	if (x$type=="com" ){
		cat("Number of complete data designs =", x$ndesignCom, "\n")
		print(tblComSum, row.names = FALSE)
		
	} else if (x$type=="miss"){
		cat("Number of missing data designs =", x$ndesignMiss, "\n")
		print(tblMissSum, row.names = FALSE)
		
	} else {
		cat("Number of complete data designs =", x$ndesignCom, "\n")
		print(tblComSum, row.names = FALSE)
		cat("\n \n")
		
		cat("Number of missing data designs =", x$ndesignMiss, "\n")
		print(tblMissSum, row.names = FALSE)
		cat("\n")
	}
}
NULL
