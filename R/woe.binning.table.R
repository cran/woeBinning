##### This is the actual function for tabulation of the binned numeric variables and factors. #####

woe.binning.table.2 <- function(pred.var, look.up.table) {

	# In case of a look-up table for a factor remove duplicates and save binned categories as row names
	if ( colnames(look.up.table)[1]=='Group.2' ) {

		look.up.table <- subset(look.up.table, !duplicated(look.up.table$Group.2))   # remove duplicates (i.e. keep only rows with merged levels, not with original levels)

		woe.table <- as.data.frame(look.up.table$Group.2)
		colnames(woe.table)[1] <- paste("Final.Bin")

		woe.table$Total.Count <- look.up.table[,5] + look.up.table[,6]
		woe.table$Total.Distr. <- woe.table$Total.Count/sum(woe.table$Total.Count)
		woe.table[,4] <- look.up.table[,5]
		colnames(woe.table)[4] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[5]),"Count",sep=".")
		woe.table[,5] <- look.up.table[,6]
		colnames(woe.table)[5] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Count",sep=".")
		woe.table[,6] <- look.up.table[,7]
		colnames(woe.table)[6] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[5]),"Distr.",sep=".")
		woe.table[,7] <- look.up.table[,8]
		colnames(woe.table)[7] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Distr.",sep=".")

		woe.table[nrow(woe.table)+1,c(2:7)] <- colSums(woe.table[c(2:7)], na.rm=TRUE)
		levels(woe.table$Final.Bin) <- c(levels(woe.table$Final.Bin), "Total")
		woe.table$Final.Bin[nrow(woe.table)] <- "Total"

		woe.table[,3] <- paste(format(round(100*woe.table[,3], 1), nsmall = 1),"%",sep="")
		woe.table[,6] <- paste(format(round(100*woe.table[,6], 1), nsmall = 1),"%",sep="")
		woe.table[,7] <- paste(format(round(100*woe.table[,7], 1), nsmall = 1),"%",sep="")
		woe.table[,8] <- woe.table[,5]/woe.table[,2]
		woe.table[,8] <- paste(format(round(100*woe.table[,8], 1), nsmall = 1),"%",sep="")
		colnames(woe.table)[8] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Rate",sep=".")
		WOE <- as.data.frame(c(look.up.table[,3],NA))
		woe.table[,9] <- format(round(WOE, 1), nsmall = 1)
		colnames(woe.table)[9] <- "WOE"
		IV <- as.data.frame(c(look.up.table[,9], look.up.table[1,4]))
		woe.table[,10] <- format(round(IV, 3), nsmall = 3)
		colnames(woe.table)[10] <- "IV"

		woe.table		

	} else {

		woe.table <- as.data.frame(row.names(look.up.table))
		colnames(woe.table)[1] <- paste("Final.Bin")
		# Make bin names more readable
		woe.table[,1] <- gsub("]", "", woe.table[,1], fixed = TRUE)
		woe.table[,1] <- gsub(" ", "", woe.table[,1], fixed = TRUE)
		woe.table[,1] <- gsub("^.*\\,","<= ",woe.table[,1])

		woe.table$Total.Count <- look.up.table[,5] + look.up.table[,6]
		woe.table$Total.Distr. <- woe.table$Total.Count/sum(woe.table$Total.Count)
		woe.table[,4] <- look.up.table[,5]
		colnames(woe.table)[4] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[5]),"Count",sep=".")
		woe.table[,5] <- look.up.table[,6]
		colnames(woe.table)[5] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Count",sep=".")
		woe.table[,6] <- look.up.table[,7]
		colnames(woe.table)[6] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[5]),"Distr.",sep=".")
		woe.table[,7] <- look.up.table[,8]
		colnames(woe.table)[7] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Distr.",sep=".")

		woe.table[nrow(woe.table)+1,c(2:7)] <- colSums(woe.table[c(2:7)], na.rm=TRUE)
		levels(woe.table$Final.Bin) <- c(levels(woe.table$Final.Bin), "Total")
		woe.table$Final.Bin[nrow(woe.table)] <- "Total"

		woe.table[,3] <- paste(format(round(100*woe.table[,3], 1), nsmall = 1),"%",sep="")
		woe.table[,6] <- paste(format(round(100*woe.table[,6], 1), nsmall = 1),"%",sep="")
		woe.table[,7] <- paste(format(round(100*woe.table[,7], 1), nsmall = 1),"%",sep="")
		woe.table[,8] <- woe.table[,5]/woe.table[,2]
		woe.table[,8] <- paste(format(round(100*woe.table[,8], 1), nsmall = 1),"%",sep="")
		colnames(woe.table)[8] <- paste(gsub("[^[:alnum:]]","",colnames(look.up.table)[6]),"Rate",sep=".")
		WOE <- as.data.frame(c(look.up.table[,1],NA))
		woe.table[,9] <- format(round(WOE, 1), nsmall = 1)
		colnames(woe.table)[9] <- "WOE"
		IV <- as.data.frame(c(look.up.table[,9], look.up.table[1,4]))
		woe.table[,10] <- format(round(IV, 3), nsmall = 3)
		colnames(woe.table)[10] <- "IV"
		
		# Remove Missing Data row in case of no NAs
		if ( woe.table[nrow(woe.table)-1,2]==0 ) { woe.table <- woe.table[-c(nrow(woe.table)-1),] }
		
		woe.table		

	}

}



#' @title Tabulation of Binning
#'
#' @description
#' \code{woe.binning.table} tabulates the binning solution generated and saved via the \code{\link{woe.binning}} function.
#'
#' @details
#' For each predictor variable \code{woe.binning.table} generates a table (data frame).
#' This table contains the final bin labels, total counts, total distribution (column percentages),
#' counts for the first and the second target class, distribution of the first and the second target
#' class (column percentages), rate (row percentages) of the target event specified via the
#' \emph{event.class} parameter in the \code{woe.binning} function, as well as weight of evidence
#' (WOE) and information values (IV).
#'
#' @usage
#' woe.binning.table(binning)
#'
#' @param binning
#' Binning information generated from the \code{\link{woe.binning}} function.
#' Contains names of the input predictor variables and the
#' corresponding binning, counts, WOE and IV information, which is used to
#' generate the tables.
#' 
#' @family woe binning functions
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
#' # Tabulate the binned variables
#' tabulate.binning <- woe.binning.table(binning)
#' tabulate.binning
#'
#' \dontrun{
#'
#' # Plot a layouted table (using the gridExtra library) for a specific
#' # variable (in this example for the first binned variable
#' # with the highest IV value)
#' library(gridExtra)
#' grid.table(tabulate.binning[[1]],
#'           theme = ttheme_default(core=list(bg_params=
#'                   list(fill=c(rep(c('grey95','grey90'),
#'	               length.out=nrow(tabulate.binning[[1]])-1),
#'                   '#BCC7BD')),fg_params=list(cex=0.8)),
#'	               colhead=list(fg_params=list(cex=0.8))),
#'           rows=NULL)
#' }
#'
#' @export

##### This function calls the actual tabulation function above for every specified predictor variable that has been binned. #####

woe.binning.table <- function(binning) {
	
	# Declare list
	woe.tables <- list()

	# Add WOE tables to list and use name of predictor variable as name of the respective list element
	if ( (length(binning)/3)==1 ) {
		woe.tables[[paste("WOE Table for",binning[1][[1]])]] <- woe.binning.table.2(binning[1][[1]], binning[2][[1]])
	} else {
		for (i in 1:(length(binning)/3)) {
				woe.tables[[paste("WOE Table for",binning[i,1][[1]])]] <- woe.binning.table.2(binning[i,1][[1]], binning[i,2][[1]])
		}
	}
	
	# Pass WOE tables
	woe.tables

}