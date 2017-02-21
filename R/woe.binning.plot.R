##### This is the actual ploting function for binned numeric variables and factors. #####

woe.binning.plot.2 <- function(pred.var, look.up.table, multiple.plots) {


	# In case of a look-up table for a factor remove duplicates and save binned categories as row names
	# In case of a look-up table for a numeric variable make bin names more readable
	if ( colnames(look.up.table)[1]=='Group.2' ) {
		look.up.table <- subset(look.up.table, !duplicated(look.up.table$Group.2))
		row.names(look.up.table) <- look.up.table$Group.2
	} else {
		row.names(look.up.table) <- gsub("]", "", row.names(look.up.table), fixed = TRUE)
		row.names(look.up.table) <- gsub(" ", "", row.names(look.up.table), fixed = TRUE)
		row.names(look.up.table) <- gsub("^.*\\,","<= ",row.names(look.up.table))
	}

	# Remove rows with NAs (for WOE values)
	look.up.table <- na.omit(look.up.table)

	# Format plot: wrapping of too long (x axis) labels
	wrap.it <- function(x, len)
	{ 
	  sapply(x, function(y) paste(strwrap(y, len), 
	                              collapse = "\n"), 
	         USE.NAMES = FALSE)
	}
	wrap.labels <- function(x, len)
	{
	  if (is.list(x))
	  {
	    lapply(x, wrap.it, len)
	  } else {
	    wrap.it(x, len)
	  }
	}
	wr.lap <- wrap.labels(row.names(look.up.table), 25)   # wrapping after 30 characters
	# Display plot
	barplot(look.up.table$woe,
		main=pred.var,
		ylab=paste("WOE"),
		names.arg=wr.lap,
		ylim=c(-100*ceiling(max(abs(look.up.table$woe),na.rm=TRUE)/100),100*ceiling(max(abs(look.up.table$woe),na.rm=TRUE)/100)),
		las=2,
		cex.main=1,
		cex.names=0.9,
		cex.axis=0.8,
		cex.lab=0.8,
		col=gray.colors(nrow(look.up.table)))
	mtext(paste("IV =", format(round(look.up.table$iv.total.final[1], 3), nsmall=3, scientific=FALSE)), line=-0.5, cex=0.8)
	
	
}



#' @title Visualization of Binning
#'
#' @description
#' \code{woe.binning.plot} visualizes the binning solution generated and saved via the \code{\link{woe.binning}} function.
#'
#' @details
#' For each predictor variable \code{woe.binning.plot} generates a weight of evidence
#' (WOE) plot. In case of multiple predictors an additional plot with variables ranked
#' via the information value (IV) will be displayed.
#'
#' @usage
#' woe.binning.plot(binning, multiple.plots, plot.range)
#'
#' @param binning
#' Binning information generated from the \code{\link{woe.binning}} function.
#' Contains names of the input predictor variables and the
#' corresponding binning, WOE and IV information, which is used to
#' generate the WOE and IV plots.
#' @param multiple.plots
#' In case the binning solution contains several predictor variables they will
#' be visualized via multiple plots (max. four WOE plots per graph window).
#' Use \emph{multiple.plots=FALSE} to avoid this and to display single plots in
#' separate windows.
#'
#' @param plot.range
#' Range of variables that should be plotted in quotes. For example \dQuote{1:10}
#' will generate WOE plots and one IV plot for the ten variables with the
#' highest IV values, \dQuote{11:20} for the next ten variables and so on.
#' Just omit this parameter to visualize all binned variables (default).
#' 
#' @family woe binning functions
#' 
#' @examples
#' # Load German credit data
#' data(germancredit)
#' df <- germancredit
#'
#' # Bin all variables of the data frame (apart from the target variable)
#' # with default parameter settings
#' binning <- woe.binning(df, 'creditability', df)
#'
#' # Plot all binned variables as multiple plots
#' woe.binning.plot(binning)
#'
#' # Plot only the first four binned variables with the highest IV value
#' # as multiple plots
#' woe.binning.plot(binning, plot.range='1:4')
#'
#' # Plot the binned variables in single plots
#' woe.binning.plot(binning, multiple.plots=FALSE)
#'
#' @importFrom stats na.omit
#' @importFrom grDevices dev.new
#' @importFrom grDevices gray.colors
#' @importFrom graphics barplot
#' @importFrom graphics par
#' @importFrom graphics mtext
#'
#' @export

##### This function calls the actual plotting function above for every specified predictor variable that has been binned. #####

woe.binning.plot <- function(binning, multiple.plots, plot.range) {

	# Set default in case multiple plots parameter is not specified
	if ( missing(multiple.plots)==TRUE ) { multiple.plots<-TRUE } 

	# Specify list of variables that should be plotted
	if ( missing(plot.range)==FALSE ) {
		binning <- binning[as.numeric(strsplit(plot.range, ":")[[1]][1]):as.numeric(strsplit(plot.range, ":")[[1]][2]),]
	}

	for (i in 1:(length(binning)/3)) {

		# Check if multiple plots (max. 4 per window) should be used
		if (  (multiple.plots==TRUE) && ((length(binning)/3)>1) ) {
			if ( (i==1) && (i<=4) ) {
				if ( (length(binning)/3)==2 ) {
					par(mfrow=c(1,2), mai=c(2.4,1,0.8,0.8))
				} else {
					par(mfrow=c(2,2), mai=c(1.6,0.8,0.4,0.4))
				}
			}
		}

		# Open a new window for each plot and set margins (in case of one variable only or multiple plot option is disabled)
		if (  (multiple.plots==FALSE) || ((length(binning)/3)==1) ) {
			dev.new()
			par(oma = c(6, 0, 0, 0))
		}
		
		# Generate WOE plot for each predictor variable
		if ( (length(binning)/3)==1 ) {
			df.input <- woe.binning.plot.2(binning[1][[1]], binning[2][[1]], multiple.plots)
		} else {
			df.input <- woe.binning.plot.2(binning[i,1][[1]], binning[i,2][[1]], multiple.plots)
		}
		
		# Check if further multiple plots (max. 4 per window) are needed		
		if ( (multiple.plots==TRUE) && ((i/4)%%1==0) && (i<(length(binning)/3)) ) {
			# Open new window for the plot and set margins
			dev.new()
			par(mfrow=c(2,2), mai=c(1.6,0.8,0.4,0.4))
		}
		
	}


	# In case of more than one predictor variable display an additional plot with variables ranked via IV
	if ( (length(binning)/3)>1 ) {

		predictor.names <- as.data.frame(unlist(binning[,1]))
		predictor.names[,1] <- strtrim(predictor.names[,1], 40)   # truncate predictor names in plot
		iv.totals.with.duplicates <- lapply(binning[,2], `[[`, 'iv.total.final')
		iv.totals <- unlist(lapply(iv.totals.with.duplicates, `[[`, 1))
		iv.table <- cbind(predictor.names,iv.totals)
		iv.table <- iv.table[order(iv.totals),]
		iv.table$combined <- paste(iv.table[,1], " IV=", format(round(iv.table$iv.totals, 3), nsmall=3, scientific=FALSE), sep="")

		# Open new window for the plot
		dev.new()   
		# Format plot: wrapping of too long axis labels
		wrap.it <- function(x, len)
		{ 
		  sapply(x, function(y) paste(strwrap(y, len), 
		                              collapse = "\n"), 
		         USE.NAMES = FALSE)
		}
		wrap.labels <- function(x, len)
		{
		  if (is.list(x))
		  {
		    lapply(x, wrap.it, len)
		  } else {
		    wrap.it(x, len)
		  }
		}

		wr.lap <- wrap.labels(iv.table$combined, 1)   # wrapping in case of blanks
		
		par(mar=c(0, 15, 3, 1), xpd=TRUE)   # defining margins
		barplot(iv.table$iv.totals,
			horiz=TRUE,
			main="Variables Ranked by Information Value",
			names.arg=wr.lap,
			xlim=c(0,ceiling(10*max(iv.table$iv.totals,na.rm=TRUE))/10),
			las=2,
			cex.main=1,
			cex.names=0.8,
			cex.axis=0.8,
			col=gray.colors(nrow(iv.table), start=0.9, end=0.3))

	}


}