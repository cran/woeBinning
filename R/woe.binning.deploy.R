##### This is the actual binning deployment function for numeric variables and factors. #####

woe.binning.deploy.2 <- function(df, pred.var, look.up.table, add.woe.or.dum.var) {


	### Binning in case a numerical variable was selected
	if ( is.numeric(df[, c(pred.var)]) ) {
	
		# Add variable with binned intervals
		dfrm.binned <- cut(df[, c(pred.var)], look.up.table[,2], labels = NULL,
			include.lowest = FALSE, right = TRUE, dig.lab = 10,
			ordered_result = FALSE)
		dfrm.binned <- as.data.frame(dfrm.binned)
		colnames(dfrm.binned)[1] <- paste(pred.var,"binned",sep=".")
		levels(dfrm.binned[,1]) <- c(levels(dfrm.binned[,1]), "Missing")   # add factor level 'Missing'
		dfrm.binned[,1][is.na(dfrm.binned[,1])] <- 'Missing'   # replace NA with string 'Missing'
		
		# Add variable with corresponding WOE values
		if ( add.woe.or.dum.var=="woe" ) {
			dfrm.binned[,2] <- look.up.table[,1][match(dfrm.binned[,1] , row.names(look.up.table))]
			colnames(dfrm.binned)[2] <- paste("woe",pred.var,"binned",sep=".")
			dfrm.binned[,2][is.na(dfrm.binned[,2])] <- look.up.table[length(look.up.table[,1]),1]   # replace NA in original numeric variable with corresponding WOE value
		}

		# Add dummy variables for binned classes
		if ( add.woe.or.dum.var=="dum" ) {
			binned.var <- dfrm.binned[,1]
			for ( level in unique(binned.var) ){
				# Remove special characters from binned intervals
				level <- gsub("(","",level, fixed = TRUE)
				level <- gsub("]","",level, fixed = TRUE)
				level <- gsub(",",".",level, fixed = TRUE)
				dfrm.binned[paste("dum",pred.var,gsub(" ","",level),"binned",sep=".")] <- ifelse(binned.var==level,1,0)
			}
		}
	
	}
	
	
	### Binning in case a factor was selected
	if ( is.factor(df[, c(pred.var)]) ) {

		# Add variable with binned levels		
		dfrm.binned <- df[, c(pred.var)]   # Copy original predictor variable
		dfrm.binned <- as.data.frame(dfrm.binned)
		levels(dfrm.binned[,1]) <- c(levels(dfrm.binned[,1]), "Missing")   # add factor level 'Missing'
		dfrm.binned[,1][is.na(dfrm.binned[,1])] <- 'Missing'   # replace NA with string 'Missing'
		dfrm.binned[,1] <- look.up.table[,1][match(dfrm.binned[,1], look.up.table[,2])]   # replace original factor level with aggregated level from look-up table
		colnames(dfrm.binned)[1] <- paste(pred.var,"binned",sep=".")
		levels(dfrm.binned[,1]) <- c(levels(dfrm.binned[,1]), "unknown")   # add factor level 'unknown'
		dfrm.binned[,1][is.na(dfrm.binned[,1])] <- "unknown"   # in case original factor level is unknown replace with "unknown"
		
		# Add variable with corresponding WOE values
		if ( add.woe.or.dum.var=="woe" ) {
			dfrm.binned[,2] <- look.up.table[,3][match(dfrm.binned[,1], look.up.table[,1])]
			colnames(dfrm.binned)[2] <- paste("woe",pred.var,"binned",sep=".")
		}

		# Add dummy variables for binned levels
		if ( add.woe.or.dum.var=="dum" ) {
			for ( level in unique(dfrm.binned[,1]) ){
				dfrm.binned[paste("dum",pred.var,gsub("[^[:alnum:]]","",level),"binned",sep=".")] <- ifelse(dfrm.binned[,1]==level,1,0)   # only alphanumeric characters are allowed
			}
		}
		
	}


	### Pass dataframe with binned variables
	dfrm.binned


}


#' @title Deployment of Binning
#'
#' @description
#' \code{woe.binning.deploy} applies the binning solution generated and saved via the \code{\link{woe.binning}}
#' or \code{\link{woe.tree.binning}} function to (new) data.
#'
#' @section General Procedure:
#' \code{woe.binning.deploy} applies the binning information that was generated from the \code{woe.binning}
#' or \code{woe.tree.binning} function to a data frame. In this data frame the names of the variables
#' to be binned need to be identical to the ones used with the \code{woe.binning}
#' or \code{woe.tree.binning} function. For each variable a binned version will be added.
#' Optionally a variable with associated weight of evidence (WOE) values or corresponding
#' dummy variables (one dummy variable for each final bin) are provided. 
#' @section Handling of Missing Data:
#' In case NAs already occurred during the \code{woe.binning} or \code{woe.tree.binning} binning process the code
#' \sQuote{Missing} is displayed and a corresponding WOE value can be computed.
#' In case NAs only occur in the deployment scenario \sQuote{Missing} is
#' displayed for numeric variables and \sQuote{unknown} for factors; and
#' the corresponding WOE values will be NAs then, as well.
#' @section Handling of Unknown Factor Levels:
#' For factor levels that have not been provided in generating the
#' binning solution via the \code{woe.binning} or \code{woe.tree.binning} function a new factor
#' level \sQuote{unknown} is displayed and the corresponding WOE value will be NA.
#'
#' @usage
#' woe.binning.deploy(df, binning, min.iv.total, add.woe.or.dum.var)
#'
#' @param df
#' Name of the data frame the binning solution - that was generated via the function \code{woe.binning}
#' or \code{woe.tree.binning} - should be applied to. The variable names and types (numerical or factor)
#' need to be identical to the ones used during the generation of the binning solution.
#' @param binning
#' Binning information generated from the \code{woe.binning} or \code{woe.tree.binning} function.
#' Contains names of the input predictor variables and the
#' corresponding binning, WOE and IV information, which is used to
#' add a binned variable to a copy of the input data.
#' @param min.iv.total
#' If the IV total value of a binned variable falls below this limit (e.g. 0.1)
#' it will not be added to the data. Just omit this parameter in case you would
#' like to add all binned variables (default).
#' @param add.woe.or.dum.var
#' \emph{add.woe.or.dum.var=\dQuote{woe}} adds an additional variable with WOE scores
#' and \emph{=\dQuote{dum}} additional dummy variables for each (aggregated) level
#' of the binned variable. In case of dummy variables make sure that you have set
#' an appropriate \emph{abbrev.fact.levels} parameter in the \code{woe.binning} or \code{woe.tree.binning} function
#' to avoid too long variable names. In principle, only alphanumeric characters
#' and dots (.) will be used for variable names. Just omit this parameter in case you
#' don't need additional variables.
#' 
#' @examples
#' # Load German credit data and create a subset
#' data(germancredit)
#' df <- germancredit[, c('creditability', 'credit.amount', 'duration.in.month',
#'                   'savings.account.and.bonds', 'purpose')]
#'
#' # Bin all variables of the data frame (apart from the target variable)
#' # with default parameter settings
#' binning <- woe.binning(df, 'creditability', df)
#'						  
#' # Deploy the binning solution to the data frame
#' # (add all binned variables and corresponding WOE variables)
#' df.with.binned.vars.added <- woe.binning.deploy(df, binning,
#'                                                add.woe.or.dum.var='woe')		
#'						  
#' # Deploy the binning solution to the data frame
#' # (add binned variables with IV>=0.1 and corresponding dummy variables)
#' df.with.binned.vars.added <- woe.binning.deploy(df, binning,
#'                                                min.iv.total=0.1,
#'                                                add.woe.or.dum.var='dum')		
#'
#' @export

##### This function calls the actual binning deployment function above for every specified predictor variable that needs to be binned. #####

woe.binning.deploy <- function(df, binning, min.iv.total, add.woe.or.dum.var) {

	#### Warning message and default in case iv.limits parameter is not (correctly) specified
	if ( !missing(min.iv.total)==TRUE ) {
		if ( min.iv.total<=0 || !is.numeric(min.iv.total) ) { warning("Incorrect parameter specification; accepted min.iv.total parameter needs to be > 0.") }
		if ( min.iv.total>max(unlist(binning[,3])) ) { warning("Incorrect parameter specification; min.iv.total parameter is > all observed IVs. Try smaller parameter value or omit parameter.") }
	} else {
		 min.iv.total=0
	}

	#### Default in case add.woe.or.dum.var parameter is not specified
	if ( missing(add.woe.or.dum.var)==TRUE ) { add.woe.or.dum.var="none" }
	
	#### Subset of binning solution with binned variables with IV total > min.iv.total
	if ( min.iv.total>0 ) {
		binning <- binning[binning[,3]>=min.iv.total,]
	}
	
	if ( (length(binning)/3) == 1 ) {
		dfrm.binned <- woe.binning.deploy.2(df, binning[1][[1]], binning[2][[1]], add.woe.or.dum.var)
		df <- cbind(df, dfrm.binned)   # add binned variables to input data frame
	} else {	
		for (i in 1:(length(binning)/3)) {
			dfrm.binned <- woe.binning.deploy.2(df, binning[i,1][[1]], binning[i,2][[1]], add.woe.or.dum.var)
			df <- cbind(df, dfrm.binned)   # add binned variables to input data frame
		}
	}

	#### Pass dataframe with binned variables
	df

}