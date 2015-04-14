#' Search for efficient designs based on design pool and hypothesized model
#'
#' This is the third step of SEEDMC. 
#' Based on the "seedpool" and "seedmod" objects 
#' generated from the first and second steps. 
#' This function automates the Mplus Monte Carlo approach 
#' to evaluate the efficiency of each design. 
#' Besides the resulting R object, 
#' a directory named "seedmc_mplus_output_yyyy-mm-dd_hh-mm-ss" 
#' will be created under the working directory. 
#' Here "yyyy-mm-dd_hh-mm-ss" reflects system time.
#' All Mplus output files are saved in the directory.
#'
#' @param pattern The "seedpool" object created by 
#' the \code{designPool()} function.
#' @param model The "seedmod" object created by 
#' the \code{modleMplus()} function.
#' @param nreps The number of replications used in the 
#' Monte carlo approach for each design. The default is 5000.
#' @param seed A number used as the seed in Mplus. 
#' If not specified, a random number will be assigned.
#' @param ... Other arguments, currently unused.
#'
#' @export
#'
#' @return Returns an object of class "seedfoud", 
#' which contains the following information.
#' \enumerate{
#'   \item Information in the given "seedpool" object. 
#'   \item Information in the given "seedmod" object.
#'   \item Information from the Mplus output 
#'		for all designs in the design pool, 
#' 		including the model results and messages 
#'		about completion error for each replication.
#' }
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
#' modQ <- modelMplus(fl = fl, latcov = latcov, 
#'		latmean = latmean, resvar = resvar) 
#'
#' seedmcoutQ <- seedmcMplus(pattern = patQ, model = modQ, nreps = 2)
#' #seedmcoutQ$results$Complete
#' seedmcoutQ$results$Complete[[1]]
#' # For presentation purpose, we set nreps = 2. 
#' # In a real study, nrep should be a large number (e.g. 5000).

seedmcMplus <- function (pattern = NULL, model = NULL, nreps = 5000, seed = NULL, ...)
{ 

	if (pattern$Time != nrow(model$fl)) {
		stop ("pattern and model do not mathch. hint: check total number of measurement occasions.")
	}
	
	if (pattern$Traj == "linear" && ncol(model$fl) != 2) {
		stop ("pattern and model do not mathch. hint: check trajectory.")
	}
	
	if (pattern$Traj == "quadratic" && ncol(model$fl) != 3) {
		stop ("pattern and model do not mathch. hint: check trajectory.")
	}
	

	parent.directory <- getwd()
    T <- as.integer(pattern$Time)
	
	if (length(seed) == 0) {
		runif(1)
		seed <- sample(.Random.seed[-1:-2], 1)
	}
	
	## seedmcMplusCom () function, nd = the number of com design
	seedmcMplusCom <- function(nd){ 
		type <- "com"
		pat <- pattern$Complete
		
		## File name 
		FNAME <- paste0("seedmc_", type , "_t", T, "_dsn", nd) 
	
	    TIMECODE <- 1:T
		
		## Varaible names	
		VNAMES <- paste0("y", TIMECODE)
	
		## PATMISS and PATPROBS
		if (length(pattern$Attrition)==0 || all(pattern$Attrition == 0)){
			patmx <- rep(0, T)
		} else {
		    patmx <- pattern$Attrition
		}
		PATMISS <- paste0(paste0("y", 1:T, "(", patmx, ")", collapse = " "), ";")  
		PATPROBS <- paste(pat$probs[nd], ";")
		
		## pop FL
        FL <- model$popMplusFL
		
		## pop LAT
		LAT <- model$popMplusLAT
		
		## pop IND
		IND <- model$popMplusIND
		
		## funtion to omit time points
		otpFUN <- function(scpt){
			sp <- list()
			short <- list()
			for (i in 1: length(scpt)) {
				if (length(TIMECODE[pat$patterns[[nd]] == 1])== 0){
					short[[i]] <- scpt[i] 
				} else {
					sp[[i]] <- strsplit(scpt[i], " ")[[1]]
					om <- paste0("y", TIMECODE[pat$patterns[[nd]] == 1])
					oml <- rep(NA, length(om))
					for (j in 1:length(om)) {
						oml[j] <- grep(om[j], sp[[i]] )
					}
					short[[i]] <- paste0(sp[[i]][-oml], collapse = " ")
				}
			}
			return(unlist(short))
		}
		anlsFL <- otpFUN(FL)
		

		indscpt_orgV <- paste(paste0("y", 1: T, "*", model$resvar, collapse = " ")) 
		indscpt_orgM <- paste0(paste0("y", 1: T, "@", model$resmean, collapse = " "))  	

		anlsIND_orgV <- otpFUN(indscpt_orgV)
		anlsIND_orgM <- otpFUN(indscpt_orgM)
		
		if (all(model$resvar[1] == model$resvar)){
			anlsIND_1 <- paste0(strsplit(anlsIND_orgV, split=" ")[[1]], "(1);")
		} else {
			anlsIND_1 <- paste0(strsplit(anlsIND_orgV, split=" ")[[1]], ";")
		}
		anlsIND_2 <- paste0("[", strsplit(anlsIND_orgM, split=" ")[[1]], "];")
		anlsIND <- c(anlsIND_1, anlsIND_2)

		
		## Mplus Syntax 
		scriptMplus <- c(
			paste0("TITLE: ", FNAME, ";"),
			"MONTECARLO: ",
			paste0("NAMES ARE ", paste(VNAMES, collapse = " "), ";"),
			paste0("NOBSERVATIONS = ", pat$N[[nd]], ";"),
			paste0("NREPS = ", nreps, ";"),
			paste0("SEED = ", seed, ";"),
			"PATMISS =",    
			PATMISS,          
			"PATPROBS =", 
			PATPROBS,
			"MODEL POPULATION:",
			FL, 
			LAT,
			IND,
			"Analysis:",
			"COVERAGE = 0;",
			"MODEL:",	
			anlsFL, 
			LAT,
			anlsIND,
			"OUTPUT: TECH9;")
		
		
		fileConn <- file(paste0(FNAME, ".inp"))
		writeLines(scriptMplus, fileConn)
		close(fileConn)

		## Run Mplus using functions from MplusAutomation package
		MplusAutomation::runModels(recursive=FALSE) # run Mplus models
				
		## Remove input and output files
		file.remove(paste0(FNAME, ".inp"))
		
	}

	## seedmcMplusMiss () function, nd = the number of miss designs
	seedmcMplusMiss <- function(nd){  
		type <- "miss"
		pat <- pattern$Missing
		
		## File name 
		FNAME <- paste0("seedmc_", type , "_t", T, "_dsn", nd) 
	
		TIMECODE <- 1:T
		
		## Varaible names	
		VNAMES <- paste0("y", TIMECODE)
	
		## PATMISS
		if (length(pattern$Attrition)==0 || all(pattern$Attrition == 0)){
			patmx <- pat$patterns[[nd]]
		} else {
		    oldpatmx <- pat$patterns[[nd]]
			patmx <- pat$patterns[[nd]]
			for (i in 1:nrow(oldpatmx)) {
				wwm <- which(oldpatmx[i,]==0)
				wwmc <- wwm - wwm[1] + 1
				patmx[i,wwm] <- oldpatmx[i, wwm] + pattern$Attrition[wwmc]
			}
		}
		PATMISS <- rep(NA, nrow(patmx))
		for (i in 1:(nrow(patmx) - 1)){
			onePat <- paste0("y", 1:T, "(", patmx[i,], ")", collapse = " ")
			PATMISS[i] <- paste(onePat, "|", sep=" ")
		}
		onePatLine <- paste0("y", 1:T, "(", patmx[nrow(patmx),], ")", collapse = " ")
		PATMISS[nrow(patmx)] <- paste(onePatLine , ";")

		
		## PATPROBS (with length controlled)
		longProbLine <- floor((nrow(patmx) - 1) / 6)
		shortProbLine <- (nrow(patmx)-1) %% 6
		if (shortProbLine == 0 ) {  
			fProbLine <- paste(rep(round(1/nrow(patmx),3), 6), "|", collapse = " ") 
			lProbLine <- paste((1 - round(1/nrow(patmx),3)* (nrow(patmx) - 1)), ";")
			PATPROBS <- c(rep(fProbLine, longProbLine), lProbLine)
		} else {
			fProbLine <- paste(rep(round(1/nrow(patmx),3), 6), "|", collapse = " ") 
			mProbLine <- paste(rep(round(1/nrow(patmx),3), nrow(patmx) - 6*longProbLine - 1), "|", collapse = " ") 
			lProbLine <- paste((1 - round(1/nrow(patmx),3)* (nrow(patmx) - 1)), ";")
			PATPROBS <- c(rep(fProbLine, longProbLine), mProbLine , lProbLine)
		} 

		
		## pop FL
		FL <- model$popMplusFL
		
		## pop LAT
		LAT <- model$popMplusLAT
		
		## pop IND
		IND <- model$popMplusIND

		## anlsIND
		indscpt_orgV <- paste(paste0("y", 1: T, "*", model$resvar, collapse = " ")) 
	
		if (all(model$resvar[1] == model$resvar)){
			anlsIND_1 <- paste0(strsplit(indscpt_orgV, split=" ")[[1]], "(1);")
		} else {
			anlsIND_1 <- paste0(strsplit(indscpt_orgV, split=" ")[[1]], ";")
		}
		anlsIND_2 <- IND[2]
		anlsIND <- c(anlsIND_1, anlsIND_2)

		
		## Mplus Syntax 
		scriptMplus <- c(
			paste0("TITLE: ", FNAME, ";"),
			"MONTECARLO: ",
			paste0("NAMES ARE ", paste(VNAMES, collapse = " "), ";"),
			paste0("NOBSERVATIONS = ", pat$N[[nd]], ";"),
			paste0("NREPS = ", nreps, ";"),
			paste0("SEED = ", seed, ";"),
			"PATMISS =",    
			PATMISS,          
			"PATPROBS =", 
			PATPROBS,
			"MODEL POPULATION:",
			FL, 
			LAT,
			IND,
			"Analysis:",
			"COVERAGE = 0;",
			"MODEL:",	
			FL, 
			LAT,
			anlsIND,
			"OUTPUT: TECH9;")
		
		
		fileConn <- file(paste0(FNAME, ".inp"))
		writeLines(scriptMplus, fileConn)
		close(fileConn)

		## Run Mplus using functions from MplusAutomation package
		MplusAutomation::runModels(recursive=FALSE) 
				
		## Remove input and output files
		file.remove(paste0(FNAME, ".inp"))
	}	
	
	readOut <- function(FNAME){
		## Number of Requested and completed replications from Mplus output 
		mplusOut <- readLines(FNAME)

		rqstRepLine <- strsplit(mplusOut[grep("Number of replications", mplusOut) + 1], "\\s")
		rqstRep <- as.numeric(rqstRepLine[[1]][rqstRepLine[[1]] != ""][2])
		
		compRepLine <- strsplit(mplusOut[grep("Number of replications", mplusOut) + 2], "\\s")
		compRep <- as.numeric(compRepLine[[1]][compRepLine[[1]] != ""][2])

		
		## Number of replications with error/warning from Mplus output 
		errMsgLine <- mplusOut[grep("REPLICATION", mplusOut)]
		errMsg <- length(errMsgLine) - 4   
		
		
		## Number of different types of error/warning messages 
		messLine1 <- mplusOut[grep("THE STANDARD ERRORS OF THE MODEL PARAMETER ESTIMATES COULD NOT", mplusOut)]
		mess1 <- length(messLine1)

		messLine2 <- mplusOut[grep("THE CHI-SQUARE STATISTIC IS NEGATIVE", mplusOut)]
		mess2 <- length(messLine2)

		messLine3 <- mplusOut[grep("THE CHI-SQUARE COULD NOT BE COMPUTED", mplusOut)]
		mess3 <- length(messLine3)
		
		messLine4 <- mplusOut[grep("NO CONVERGENCE", mplusOut)]
		mess4 <- length(messLine4)
		
		messLine5 <- mplusOut[grep("THE LATENT VARIABLE COVARIANCE MATRIX", mplusOut)]
		mess5 <- length(messLine5)
		
		messLine6 <- mplusOut[grep("THE RESIDUAL COVARIANCE MATRIX", mplusOut)]
		mess6 <- length(messLine6)
		
		extraInfo <- c(rqstrep =rqstRep, comprep = compRep, errmsg = errMsg,
						mess1 = mess1 , mess2 = mess2, mess3 = mess3,
						mess4 = mess4 , mess5 = mess5, mess6 = mess6)
	} 
	

	DRNAME_TIME0 <- Sys.time()
	DRNAME_TIME1 <- gsub(" ", "_", DRNAME_TIME0)
	DRNAME_TIME2 <- gsub(":", "-", DRNAME_TIME1)
	DRNAME <- paste0("seedmc_mplus_output", "_", DRNAME_TIME2)

	dir.create(paste0(parent.directory, "/", DRNAME), showWarnings = TRUE) 
	if (pattern$Type == "com") { 
		dir.create(paste0(parent.directory, "/", DRNAME, "/com")) 
		setwd(paste0(parent.directory, "/", DRNAME, "/com"))
		comRun <- sapply(1:length(pattern$Complete$patterns), seedmcMplusCom)
		comEst <- MplusAutomation::extractModelParameters(recursive=FALSE) 
		comExtra <- t(sapply(names(comEst), readOut))
		setwd(parent.directory)
		
		missEst <- NULL
		missExtra <- NULL
		
	} else if (pattern$Type == "miss") {  
		comEst <- NULL
		comExtra <- NULL
		
		dir.create(paste0(parent.directory, "/", DRNAME, "/miss"))  
		setwd(paste0(parent.directory, "/", DRNAME, "/miss"))
		missRun <- sapply(1:length(pattern$Missing$patterns), seedmcMplusMiss)
		missEst <- MplusAutomation::extractModelParameters(recursive=FALSE) 
		missExtra <- t(sapply(names(missEst), readOut))
		setwd(parent.directory)	
		
	} else if (pattern$Type == "both")  {
		dir.create(paste0(parent.directory, "/", DRNAME, "/com")) 
		setwd(paste0(parent.directory, "/", DRNAME, "/com"))
		comRun <- sapply(1:length(pattern$Complete$patterns), seedmcMplusCom)
		comEst <- MplusAutomation::extractModelParameters(recursive=FALSE)
		comExtra <- t(sapply(names(comEst), readOut))
		
		dir.create(paste0(parent.directory, "/", DRNAME, "/miss")) 
		setwd(paste0(parent.directory, "/", DRNAME, "/miss"))
		missRun <- sapply(1:length(pattern$Missing$patterns), seedmcMplusMiss)
		missEst <- MplusAutomation::extractModelParameters(recursive=FALSE)
		missExtra <- t(sapply(names(missEst), readOut))
	}
	
	out <- list(patterns = pattern, model = model, 
				results = list(Complete = comEst, Missing = missEst), 
				extraInfo = list(Complete = comExtra, Missing = missExtra)) 
				
				
	class(out) <- "seedfound"
	
	setwd(paste0(parent.directory, "/", DRNAME)) 
	dput(out, file = "seedmc_mplus_res.RData") 
	setwd(parent.directory) 
	
	return(out)
}
NULL

#' @export
print.seedfound <-  function(x, ...) 
{  
	outPool <- list(Type = x$patterns$Type, Traj = x$patterns$Traj, Time = x$patterns$Time,
				Attrition = x$patterns$Attrition) ### 2015/04/04
				
	cat("Design Pool: \n")
	print(outPool)
	cat("\n")

	outMod <- list(fl = x$model$fl, latcov = x$model$latcov, latmean = x$model$latmean, 
				resvar = x$model$resvar, resmean = x$model$resmean)
	cat("Model Population: \n")
	print(outMod)
	cat("\n")
	
	if (x$patterns$Type=="com" || x$patterns$Type=="both" ){
		cat("To check simulation results of complete data designs, call element $results$Complete. \n")
	}
	
	if (x$patterns$Type=="miss" || x$patterns$Type=="both" ){
		cat("To check simulation results of each missing data design, call element $results$Missing. \n")
	}
	
	cat("Use function topDesignsMplus() for further information. \n")
 }
NULL