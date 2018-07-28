##### This is the actual tree-like binning function for numeric variables and factors. #####

woe.tree.binning.2 <- function(df, target.var, pred.var, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, bad, good) {


#### Build subsets with target and predictor variable
df <- df[, c(target.var, pred.var)]  # used for final binning
dfrm <- df[, c(target.var, pred.var)]   # used for iterative merging and splitting of bins
colnames(dfrm)[1] <- paste("target.var")
colnames(dfrm)[2] <- paste("predictor.var")


#### Check if numerical variable or factor was provided as predictor and apply appropriate binning technique

### Binning in case a numerical variable was selected
if ( length(unique(dfrm[,1]))==2 && is.numeric(dfrm[,2]) ) {

	## Derive number of initial bins from min.perc.total parameter
	max.bins <- trunc(1/min.perc.total)

	## Derive cutpoints for bins (with similar frequency)
	cutpoints <- quantile(dfrm$predictor.var,(0:max.bins)/max.bins, na.rm=TRUE)
	innercutpoints <- cutpoints[2:(length(cutpoints)-1)]   # remove outer (observed) boudaries
	cutpoints <- c(-Inf, innercutpoints, +Inf)   # add -Inf, +Inf to cutpoints
	cutpoints <- unique(cutpoints)   # remove multiple cutpoints with same value

	## Calculate initial crosstab from binned variable and target variable
	## to identify and merge sparse bins
	
	# Compute binned variable from cutpoints and add it to the subset data frame
	dfrm$predictor.var.binned <- cut(dfrm$predictor.var, cutpoints, labels = NULL,
    		include.lowest = FALSE, right = TRUE, dig.lab = 10,
		ordered_result = TRUE)
	
	# Compute crosstab from binned variable and target variable and covert it to a data frame
	freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var, useNA="always")
	row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
	woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
	woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
	# Compute columns percents for target classes from crosstab frequencies
	woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
	woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
	# Correct column percents in case of 0 frequencies (in case of no NA skip last row)
	if ( !anyNA(df[,2]) ) {
		if ( min(woe.dfrm[-nrow(woe.dfrm),1],na.rm=TRUE)==0 | min(woe.dfrm[-nrow(woe.dfrm),2],na.rm=TRUE)==0 ) {
			woe.dfrm$col.perc.a[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)
			woe.dfrm$col.perc.b[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)	
		}
	} else {
		if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
			woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
			woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
		}
	}

	# Check for bins (without last regular and without NA bin) if frequencies < percentage limit specified above
	# (in reverse order to remain correct reference to cutpoints)
	for (i in (nrow(woe.dfrm)-2):1) {
		if (woe.dfrm$col.perc.a[i]<min.perc.class | woe.dfrm$col.perc.b[i]<min.perc.class | ((woe.dfrm[i,1]+woe.dfrm[i,2])/(sum(woe.dfrm[,1],na.rm=TRUE)+sum(woe.dfrm[,2],na.rm=TRUE)))<min.perc.total) {
			# Remove cutpoint			
			cutpoints <- cutpoints[-c((i+1))]
			# Compute binned variable from cutpoints and add it to the subset data frame
			dfrm$predictor.var.binned <- cut(dfrm$predictor.var, cutpoints, labels = NULL,
					include.lowest = FALSE, right = TRUE, dig.lab = 10,
					ordered_result = TRUE)
			# Compute crosstab from binned variable and target variable and covert it to a data frame
			freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var, useNA="always")
			row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
			woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
			woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
			# Compute columns percents for target classes from crosstab frequencies
			woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
			woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
			# Correct column percents in case of 0 frequencies (in case of no NA skip last row)
			if ( !anyNA(df[,2]) ) {
				if ( min(woe.dfrm[-nrow(woe.dfrm),1],na.rm=TRUE)==0 | min(woe.dfrm[-nrow(woe.dfrm),2],na.rm=TRUE)==0 ) {
					woe.dfrm$col.perc.a[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)
					woe.dfrm$col.perc.b[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)	
				}
			} else {
				if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
					woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
					woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
				}
			}
		}
		# Stop in case 3 cutpoints (-Inf, x, +Inf) are reached
		if ( length(cutpoints)==3 ) { break } 
	}
	
	# Check for last regular bin if frequencies < percentage limit specified above (only in case number of cutpoints > 3
	if ( length(cutpoints)>3 ) {
		if (woe.dfrm$col.perc.a[(nrow(woe.dfrm)-1)]<min.perc.class | woe.dfrm$col.perc.b[(nrow(woe.dfrm)-1)]<min.perc.class | ((woe.dfrm[nrow(woe.dfrm)-1,1]+woe.dfrm[nrow(woe.dfrm)-1,2])/(sum(woe.dfrm[,1],na.rm=TRUE)+sum(woe.dfrm[,2],na.rm=TRUE)))<min.perc.total) {
			# Remove cutpoint
			cutpoints <- cutpoints[-c(nrow(woe.dfrm)-1)]
			# Compute binned variable from cutpoints and add it to the subset data frame
			dfrm$predictor.var.binned <- cut(dfrm$predictor.var, cutpoints, labels = NULL,
					include.lowest = FALSE, right = TRUE, dig.lab = 10,
					ordered_result = TRUE)
			# Compute crosstab from binned variable and target variable and covert it to a data frame
			freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var, useNA="always")
			row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
			woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
			woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
			# Compute columns percents for target classes from crosstab frequencies
			woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
			woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
			# Correct column percents in case of 0 frequencies (in case of no NA skip last row)
			if ( !anyNA(df[,2]) ) {
				if ( min(woe.dfrm[-nrow(woe.dfrm),1],na.rm=TRUE)==0 | min(woe.dfrm[-nrow(woe.dfrm),2],na.rm=TRUE)==0 ) {
					woe.dfrm$col.perc.a[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)
					woe.dfrm$col.perc.b[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)	
				}
			} else {
				if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
					woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
					woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
				}
			}
		}
	}


	## After sparse levels are merged:
	## Tree-based iterative partitioning of bins until IV based stop criteria is reached
	## or 2 aggregated bins are left (i.e. 3 cutpoints: -Inf, middle cutpoint, +Inf).

	innercutpoints <- cutpoints[2:(length(cutpoints)-1)]
	
	if ( length(cutpoints)>2 ) {
	
		for (i in 1:(length(innercutpoints)-1)) {

			for (i in 1:length(innercutpoints)) {
			
				if ( exists('selected.cuts', inherits=FALSE) ) {
					pred.var.cut <- cut(dfrm$predictor.var, c(-Inf, selected.cuts, innercutpoints[i], Inf), labels=NULL, include.lowest=FALSE, right=TRUE, dig.lab=10, ordered_result=TRUE)
				} else {
					pred.var.cut <- cut(dfrm$predictor.var, c(-Inf, innercutpoints[i], Inf), labels=NULL, include.lowest=FALSE, right=TRUE, dig.lab=10, ordered_result=TRUE)
				}
			
				freq.table <- table(pred.var.cut, dfrm$target.var, useNA="always")
				row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
				woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
				woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
				woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
				woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
				# Correct column percents in case of 0 frequencies (in case of no NA skip last row)
				if ( !anyNA(df[,2]) ) {
					if ( min(woe.dfrm[-nrow(woe.dfrm),1],na.rm=TRUE)==0 | min(woe.dfrm[-nrow(woe.dfrm),2],na.rm=TRUE)==0 ) {
						woe.dfrm$col.perc.a[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.a[-nrow(woe.dfrm)]+0.0001)
						woe.dfrm$col.perc.b[-nrow(woe.dfrm)] <- (woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)/sum(woe.dfrm$col.perc.b[-nrow(woe.dfrm)]+0.0001)	
					}
				} else {
					if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
						woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
						woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
					}
				}
				woe.dfrm$woe <- 100*log(woe.dfrm$col.perc.a/woe.dfrm$col.perc.b)
				woe.dfrm$woe[is.finite(woe.dfrm$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
				woe.dfrm$iv.bins <- (woe.dfrm$col.perc.a-woe.dfrm$col.perc.b)*woe.dfrm$woe/100
				iv.total <- sum(woe.dfrm$iv.bins, na.rm=TRUE)
				ifelse (exists('iv.total.collect', inherits=FALSE), iv.total.collect <- cbind(iv.total.collect, iv.total), iv.total.collect <- iv.total)
				
			}

			# Restore former solution in case stop criteria is reached and exit loop
			if ( exists('max.iv.total.collect.backup', inherits=FALSE) ) {
				if ( (max.iv.total.collect.backup+max.iv.total.collect.backup*stop.limit)>max(iv.total.collect) ) {
					innercutpoints <- innercutpoints.backup
					break
				}
			}

			# Backups to be able to restore former solution in case stop criteria is reached
			max.iv.total.collect.backup <- max(iv.total.collect)
			innercutpoints.backup <- innercutpoints
			
			# Get index of cutpoint with highest IV and reset iv.total.collect
			index.optimal.cut <- which(iv.total.collect==max(iv.total.collect))[1]
			iv.total.collect <- NULL

			# collect and sort selected cuts
			ifelse (exists('selected.cuts', inherits=FALSE), selected.cuts <- cbind(selected.cuts, innercutpoints[index.optimal.cut[sort.list(index.optimal.cut)]]), selected.cuts <- innercutpoints[index.optimal.cut[sort.list(index.optimal.cut)]])
			selected.cuts <- sort(selected.cuts)
			selected.cuts <- unique(selected.cuts)

			# Remove selected cutpoint from cutpoint list
			innercutpoints <- innercutpoints[-index.optimal.cut]

		}

			#print(selected.cuts)
			pred.var.cut <- cut(dfrm$predictor.var, c(-Inf, selected.cuts, Inf), labels = NULL, include.lowest = FALSE, right = TRUE, dig.lab = 10, ordered_result=TRUE)
			freq.table <- table(pred.var.cut, dfrm$target.var, useNA="always")
			row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
			woe.dfrm.final <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
			woe.dfrm.final <- woe.dfrm.final[, c(good, bad)]   # Select columns with raw frequencies only			
			woe.dfrm.final$col.perc.a <- woe.dfrm.final[,1]/sum(woe.dfrm.final[,1])
			woe.dfrm.final$col.perc.b <- woe.dfrm.final[,2]/sum(woe.dfrm.final[,2])
			# Correct column percents in case of 0 frequencies (in case of no NA skip last row)
			if ( !anyNA(df[,2]) ) {
				if ( min(woe.dfrm.final[-nrow(woe.dfrm.final),1],na.rm=TRUE)==0 | min(woe.dfrm.final[-nrow(woe.dfrm.final),2],na.rm=TRUE)==0 ) {
					woe.dfrm.final$col.perc.a[-nrow(woe.dfrm.final)] <- (woe.dfrm.final$col.perc.a[-nrow(woe.dfrm.final)]+0.0001)/sum(woe.dfrm.final$col.perc.a[-nrow(woe.dfrm.final)]+0.0001)
					woe.dfrm.final$col.perc.b[-nrow(woe.dfrm.final)] <- (woe.dfrm.final$col.perc.b[-nrow(woe.dfrm.final)]+0.0001)/sum(woe.dfrm.final$col.perc.b[-nrow(woe.dfrm.final)]+0.0001)	
				}
			} else {
				if ( min(woe.dfrm.final[,1],na.rm=TRUE)==0 | min(woe.dfrm.final[,2],na.rm=TRUE)==0 ) {
					woe.dfrm.final$col.perc.a <- (woe.dfrm.final$col.perc.a+0.0001)/sum(woe.dfrm.final$col.perc.a+0.0001)
					woe.dfrm.final$col.perc.b <- (woe.dfrm.final$col.perc.b+0.0001)/sum(woe.dfrm.final$col.perc.b+0.0001)	
				}
			}
		
	}

	woe.dfrm.final$woe <- 100*log(woe.dfrm.final$col.perc.a/woe.dfrm.final$col.perc.b)
	woe.dfrm.final$woe[is.finite(woe.dfrm.final$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
	woe.dfrm.final$iv.bins <- (woe.dfrm.final$col.perc.a-woe.dfrm.final$col.perc.b)*woe.dfrm.final$woe/100
	# Add cutpoints needed for deployment
	cutpoints.final <- c(-Inf, selected.cuts, Inf)
	woe.dfrm.final$cutpoints.final <- cutpoints.final
	upper.cutpoints.final.dfrm <- rbind(as.data.frame(cutpoints.final[-1]),'Missing')
	woe.dfrm.final <- cbind(woe.dfrm.final, upper.cutpoints.final.dfrm)
	# Compute final IV
	iv.total.final <- sum(woe.dfrm.final$iv.bins, na.rm=TRUE)
	woe.dfrm.final$iv.total.final <- iv.total.final
	## Save final binning solution via look-up-table for deployment
	look.up.table <- woe.dfrm.final[,c(5,7:9,1:4,6)]

}


### Binning in case a factor was selected
if ( length(unique(dfrm[,1]))==2 && is.factor(dfrm[,2]) ) {
	
	## Copy predictor variable to prepare binning/recoding
	dfrm$predictor.var.binned <- dfrm$predictor.var

	## Handling of NAs
	if ( anyNA(dfrm$predictor.var.binned)==TRUE ) {
		levels(dfrm$predictor.var.binned) <- c(levels(dfrm$predictor.var.binned), "Missing")   # add factor level 'Missing'
		dfrm$predictor.var.binned[is.na(dfrm$predictor.var.binned)] <- "Missing"   # replace NA with string 'Missing'
	}
	
	## Prepare binned factor in INPUT data (levels may be merged in subsequent steps)
	df[,ncol(df)+1] <- df[, c(pred.var)]
	colnames(df)[ncol(df)] <- paste(pred.var,".binned",sep="")
	# Handling of NAs
	if ( anyNA(df[,ncol(df)])==TRUE ) {
		levels(df[,ncol(df)]) <- c(levels(df[,ncol(df)]), "Missing")   # add factor level 'Missing'
		df[,ncol(df)][is.na(df[,ncol(df)])] <- "Missing"   # replace NA with string 'Missing'
	}

	
	## Calculate initial crosstab from binned variable and target variable
	## to identify and merge sparse bins
	
	# Compute crosstab from binned variable and target variable and covert it to a data frame
	freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var)
	woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
	woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only

	# Compute WOE and information value (IV) from crosstab frequencies
	woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
	woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
	# Correct column percents in case of 0 frequencies
	if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
		woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
		woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
	}

	# Merge factor levels with frequencies < percentage limit specified above to "misc. level" (associated with pos. and neg. WOE values)
	woe.dfrm$sparse.merge[woe.dfrm$col.perc.a<min.perc.class | woe.dfrm$col.perc.b<min.perc.class | ((woe.dfrm[,1]+woe.dfrm[,2])/(sum(woe.dfrm[,1],na.rm=TRUE)+sum(woe.dfrm[,2],na.rm=TRUE)))<min.perc.total] <- 1
	woe.dfrm.sparse.subset <- na.omit(woe.dfrm)
	woe.dfrm.sparse.subset$sparse.merge[woe.dfrm.sparse.subset$col.perc.a <= woe.dfrm.sparse.subset$col.perc.b] <- -1
	woe.dfrm.sparse.subset.pos <- woe.dfrm.sparse.subset[woe.dfrm.sparse.subset$sparse.merge==1, ]
	woe.dfrm.sparse.subset.neg <- woe.dfrm.sparse.subset[woe.dfrm.sparse.subset$sparse.merge==-1, ]
	levels(dfrm$predictor.var.binned)[levels(dfrm$predictor.var.binned)%in%(row.names(woe.dfrm.sparse.subset.pos))] <- "misc. level pos."
	levels(dfrm$predictor.var.binned)[levels(dfrm$predictor.var.binned)%in%(row.names(woe.dfrm.sparse.subset.neg))] <- "misc. level neg."


	## After sparse levels are merged:
	## Tree-based partitioning of bins sorted by WOE vales

	# Compute crosstab from binned variable and target variable and covert it to a data frame
	freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var)
	#row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
	woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
	woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
	# Compute WOE and information value (IV) from crosstab frequencies
	woe.dfrm$col.perc.a <- woe.dfrm[,1]/sum(woe.dfrm[,1])
	woe.dfrm$col.perc.b <- woe.dfrm[,2]/sum(woe.dfrm[,2])
	# Correct column percents in case of 0 frequencies
	if ( min(woe.dfrm[,1],na.rm=TRUE)==0 | min(woe.dfrm[,2],na.rm=TRUE)==0 ) {
		woe.dfrm$col.perc.a <- (woe.dfrm$col.perc.a+0.0001)/sum(woe.dfrm$col.perc.a+0.0001)
		woe.dfrm$col.perc.b <- (woe.dfrm$col.perc.b+0.0001)/sum(woe.dfrm$col.perc.b+0.0001)	
	}
	woe.dfrm$woe <- 100*log(woe.dfrm$col.perc.a/woe.dfrm$col.perc.b)
	woe.dfrm$woe[is.finite(woe.dfrm$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
	woe.dfrm <- woe.dfrm[order(woe.dfrm$woe),]   # sort data via WOE values
	woe.dfrm$iv.bins <- (woe.dfrm$col.perc.a-woe.dfrm$col.perc.b)*woe.dfrm$woe/100

	# In case there are more than 2 regulare bins (+ Missing bin) left:
	# iterative split bins into binary subsets (tree-like, i.e. 1. split
	# -> 2 aggregated bins, 2. split -> 3 aggregated bins, etc.) and realize
	# solution with total IV value that fullfills the stop crieria
	if ( (anyNA(df[,2]) && nrow(woe.dfrm)>3) || (!anyNA(df[,2]) && nrow(woe.dfrm)>2) )  {
		for ( i in 1:1:(nrow(woe.dfrm-2)) ) {
			for ( i in 1:(nrow(woe.dfrm)-1) ) {
				woe.dfrm$trycut[1:i] <- 'a'   # 1. node
				woe.dfrm$trycut[(i+1):nrow(woe.dfrm)] <- 'b'   # 2. node
				if ( !'cut' %in% names(woe.dfrm) ) {
					woe.dfrm.try <- aggregate(woe.dfrm[,3:4], by=list(woe.dfrm$trycut), 'sum')
				} else {
					woe.dfrm.try <- aggregate(woe.dfrm[,3:4], by=list(woe.dfrm$trycut, woe.dfrm$cut), 'sum')
				}
				woe.dfrm.try$woe <- 100*log(woe.dfrm.try$col.perc.a/woe.dfrm.try$col.perc.b)
				woe.dfrm.try$woe[is.finite(woe.dfrm.try$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
				woe.dfrm.try$iv.bins <- (woe.dfrm.try$col.perc.a-woe.dfrm.try$col.perc.b)*woe.dfrm.try$woe/100
				iv.total <- sum(woe.dfrm.try$iv.bins, na.rm=TRUE)
				ifelse (exists('iv.total.collect', inherits=FALSE), iv.total.collect <- cbind(iv.total.collect, iv.total), iv.total.collect <- iv.total)
			}
			
			index.optimal.cut <- which(iv.total.collect==max(iv.total.collect))[1]
	
			# Restore former solution in case stop criteria is reached and exit loop
			if ( exists('max.iv.total.collect.backup', inherits=FALSE) ) {
				if ( (max.iv.total.collect.backup+max.iv.total.collect.backup*stop.limit)>max(iv.total.collect) ) {
					break
				}
			}	
		
			# Backup to be able to restore former solution in case stop criteria is reached
			max.iv.total.collect.backup <- max(iv.total.collect)
		
			iv.total.collect <- NULL
			
			woe.dfrm$cuttemp <- 'm'   # all incl. Missing
			woe.dfrm$cuttemp[1:index.optimal.cut] <- 'a'   # 1. node
			woe.dfrm$cuttemp[(index.optimal.cut+1):nrow(woe.dfrm)] <- 'b'   # 2. node
			
			if ( !'cut' %in% names(woe.dfrm) ) {
				woe.dfrm$cut <- woe.dfrm$cuttemp
			} else {
				woe.dfrm$cut <- paste(woe.dfrm$cut, woe.dfrm$cuttemp, sep="")
			}
			
		}
	}

	woe.dfrm$Group.1 <- row.names(woe.dfrm)
	woe.dfrm$Group.2 <- row.names(woe.dfrm)
	
	# Merge names of factor levels to be joined in a new variable
	if ( (anyNA(df[,2]) && nrow(woe.dfrm)>3) || (!anyNA(df[,2]) && nrow(woe.dfrm)>2) )  {
		for ( i in (nrow(woe.dfrm)-1):1 ) {
			if ( woe.dfrm$cut[i]==woe.dfrm$cut[i+1] ) {
				woe.dfrm$Group.2[i] <- paste(row.names(woe.dfrm)[i], "+", woe.dfrm$Group.2[i+1])
			}
		}
		for ( i in 2:nrow(woe.dfrm) ) {
			if ( woe.dfrm$cut[i]==woe.dfrm$cut[i-1] ) {
				woe.dfrm$Group.2[i] <- woe.dfrm$Group.2[i-1]
			}
		}
	} else {   # In case of only 2 regular bins (+ 1 missing data bin) build the data frame structure that is expected by the final procedure
		woe.dfrm$trycut <- NA
		woe.dfrm$cuttemp <- NA
		woe.dfrm$cut <- "a"
		woe.dfrm$cut[2] <- "b"
		if ( nrow(woe.dfrm)>2 ) { woe.dfrm$cut[3] <- "c" }
		woe.dfrm <- woe.dfrm[,c(1:6,9:11,7,8)]
	}

	# Restore original factor level names and original counts via outer join (because they may have be lost by former aggregating to misc. levels)
	woe.dfrm.sparse.subset$misc[woe.dfrm.sparse.subset$sparse.merge==1] <- "misc. level pos."
	woe.dfrm.sparse.subset$misc[woe.dfrm.sparse.subset$sparse.merge==-1] <- "misc. level neg."
	woe.dfrm.sparse.subset$original.names <- row.names(woe.dfrm.sparse.subset)		
	# Rename variables with aggregated count vor misc. bins to avoid name conflicts in merging
	colnames(woe.dfrm)[1] <- paste(colnames(woe.dfrm)[1], "aggr", sep=".")
	colnames(woe.dfrm)[2] <- paste(colnames(woe.dfrm)[2], "aggr", sep=".")
	# Merge
	woe.dfrm <- merge( woe.dfrm.sparse.subset[,c(6:7,1:2)], woe.dfrm, by.x=1, by.y=10, all=TRUE)
	# Restore original factor level names
	woe.dfrm$Group.1 <- woe.dfrm$misc
	woe.dfrm$Group.1[!is.na(woe.dfrm$original.names)] <- woe.dfrm$original.names[!is.na(woe.dfrm$original.names)]
	# Restore original counts
	woe.dfrm[,3][is.na(woe.dfrm[,3])] <- woe.dfrm[,5][is.na(woe.dfrm[,3])]
	woe.dfrm[,4][is.na(woe.dfrm[,4])] <- woe.dfrm[,6][is.na(woe.dfrm[,4])]
	# Remove unnecessary count variables
	woe.dfrm <- woe.dfrm[, -c(5,6)]

	# Realize final bin aggregation resulting from the tree-like procedure above
	# and compute corresponding WOE and IV values
	woe.dfrm.aggr <- aggregate(woe.dfrm[,3:4], by=list(woe.dfrm$cut), 'sum')
	colnames(woe.dfrm.aggr)[1] <- 'cut'
	woe.dfrm.aggr$col.perc.a <- woe.dfrm.aggr[,2]/sum(woe.dfrm.aggr[,2])
	woe.dfrm.aggr$col.perc.b <- woe.dfrm.aggr[,3]/sum(woe.dfrm.aggr[,3])

	# Correct column percents in case of 0 frequencies
	if ( min(woe.dfrm.aggr[,2],na.rm=TRUE)==0 | min(woe.dfrm.aggr[,3],na.rm=TRUE)==0 ) {
		woe.dfrm.aggr$col.perc.a <- (woe.dfrm.aggr$col.perc.a+0.0001)/sum(woe.dfrm.aggr$col.perc.a+0.0001)
		woe.dfrm.aggr$col.perc.b <- (woe.dfrm.aggr$col.perc.b+0.0001)/sum(woe.dfrm.aggr$col.perc.b+0.0001)	
	}	
	woe.dfrm.aggr$woe <- 100*log(woe.dfrm.aggr$col.perc.a/woe.dfrm.aggr$col.perc.b)
	woe.dfrm.aggr$woe[is.finite(woe.dfrm.aggr$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
	woe.dfrm.aggr <- woe.dfrm.aggr[order(woe.dfrm.aggr$woe),]   # sort data via WOE values
	woe.dfrm.aggr$iv.bins <- (woe.dfrm.aggr$col.perc.a-woe.dfrm.aggr$col.perc.b)*woe.dfrm.aggr$woe/100
	woe.dfrm.aggr$iv.total.final <- sum(woe.dfrm.aggr$iv.bins, na.rm=TRUE)

	# Merge the table with the final WOE and IV values with the table containing the original and aggregated bin names
	look.up.table <- merge(woe.dfrm.aggr, woe.dfrm[11:13], by.x=1, by.y=1)
	look.up.table <- look.up.table[,c(9,10,6,8,2:5,7)]
	look.up.table <- look.up.table[order(look.up.table$woe, look.up.table$Group.2),]   # sort by woe value and merged bin name

	# Convert variables with original and aggregated factor levels into factors
	look.up.table$Group.1 <- factor(look.up.table$Group.1)
	look.up.table$Group.2 <- factor(look.up.table$Group.2)
	
	# In case the misc. level consists only of only NA rename it 'Missing'
	if ( length(which(look.up.table[,2]=='Missing'))==1 && length(which(look.up.table[,1]=="misc. level neg."))==1 ) {
		if ( (which(look.up.table[,2]=='Missing') == which(look.up.table[,1]=='misc. level neg.')) ) {
			levels(look.up.table[,1]) <- c(levels(look.up.table[,2]), 'Missing')   # add factor level 'Missing'
			look.up.table[,1][look.up.table[,2]=='Missing'] <- 'Missing'
		}
	}
	if ( length(which(look.up.table[,2]=='Missing'))==1 && length(which(look.up.table[,1]=="misc. level pos."))==1 ) {
		if ( (which(look.up.table[,2]=='Missing') == which(look.up.table[,1]=='misc. level pos.')) ) {
			levels(look.up.table[,1]) <- c(levels(look.up.table[,2]), 'Missing')   # add factor level 'Missing'
			look.up.table[,1][look.up.table[,2]=='Missing'] <- 'Missing'
		}
	}

	# Abbreviate long factor levels (in case they are longer than specified or longer than 1000 characters)
	if ( abbrev.fact.levels==0 && 1000<max(nchar(as.character(look.up.table$Group.2))) ) {
		abbrev.fact.levels <- 1000
	}	
	if ( abbrev.fact.levels>0 && abbrev.fact.levels<max(nchar(as.character(look.up.table$Group.2))) ) {
		look.up.table$Group.2 <- as.factor(abbreviate(look.up.table$Group.2, abbrev.fact.levels))   # actual abbrevation
		look.up.table$Group.2 <- as.factor(gsub("[*+*]", " ", look.up.table$Group.2))   # remove + signs
		look.up.table$Group.2 <- as.factor(gsub("  +", " ", look.up.table$Group.2))   # remove double blanks
	}

}


#### Check for correct variable specification and
#### generate requested output, in case specification is correct

### Display warning message in case of incorrect predictor variable specification

if ( (is.numeric(dfrm[,2])==FALSE) && (is.factor(dfrm[,2])==FALSE)  ) {
	warning("Incorrect variable specification.\nPredictor variable needs to be a numeric variable or a factor.")
}

### Generate requested output, in case specification is correct

else {

	## Function passes the final binning solution as look-up table
	look.up.table
	
}


}



#' @title Binning via Tree-Like Segmentation
#'
#' @description
#' \code{woe.tree.binning} generates a supervised tree-like segmentation of numeric variables
#' and factors with respect to a dichotomous target variable. Its parameters provide
#' flexibility in finding a binning that fits specific data characteristics and practical
#' needs.
#'
#' @section Binning of Numeric Variables:
#' Numeric variables (continuous and ordinal) are binned beginning with initial classes with
#' similar frequencies. The number of initial bins results from the \emph{min.perc.total}
#' parameter: min.perc.total will result in trunc(1/min.perc.total) initial bins,
#' whereby \emph{trunc} is needed to guarantee bins with similar frequencies.
#' For example \emph{min.perc.total=0.07} will cause trunc(14.3)=14 initial classes.
#' Next, if \emph{min.perc.class}>0, bins with sparse target classes will be merged with
#' the next upper bin, and in case of the last bin with the next lower one. NAs have
#' their own bin and will not be merged with others. Finally the actual tree-like procedure
#' starts: binary splits iteratively assign nearby classes with similar weight of evidence
#' (WOE) values to segments in a way that maximizes the resulting information value (IV).
#' The procedure stops when the IV increases less then specified by a percentage value
#' (\emph{stop.limit} parameter).
#' @section Binning of Factors:
#' Factors (categorical variables) are binned via factor levels. As a start sparse levels
#' (defined via the \emph{min.perc.total} and \emph{min.perc.class} parameters) are merged
#' to a \sQuote{miscellaneous} level: if possible, respective levels (including sparse NAs)
#' are bundled as \sQuote{misc. level pos.} (associated with positive WOE values), respectively
#' as \sQuote{misc. level neg.} (associated with negative WOE values). In case a misc. level
#' contains only NAs it will be named \sQuote{Missing}. Afterwards the actual tree-like
#' procedure starts: binary splits iteratively assign levels with similar WOE values to
#' segments in a way that maximizes the resulting information value (IV). The procedure stops
#' when the IV increases less then specified by a percentage value (\emph{stop.limit} parameter).
#' @section Adjustment of 0 Frequencies:
#' In case the crosstab of the bins with the target classes contains frequencies = 0
#' the column percentages are adjusted to be able to compute the WOE and IV values:
#' the offset 0.0001 (=0.01\%) is added to each column percentage cell and the column
#' percentages are recomputed then. This allows considering bins associated with one target
#' class only, but may cause extreme WOE values for these bins. If a correction is not
#' appropriate choose \emph{min.perc.class}>0; bins with sparse target classes will be
#' merged then before computing any WOE or IV value.
#' @section Handling of Missing Data:
#' Cases with NAs in the target variable will be ignored. For predictor variables the following
#' applies: in case NAs already occurred when generating the binning solution
#' the code \sQuote{Missing} is displayed and a corresponding WOE value can be computed.
#' (Note that factor NAs may be joined with other sparse levels to a \sQuote{miscellaneous}
#' level - see above; only this \sQuote{miscellaneous} level will be displayed then.)
#' In case NAs occur in the deployment scenario only \sQuote{Missing} is
#' displayed for numeric variables and \sQuote{unknown} for factors; and
#' the corresponding WOE values will be NA then, as well.
#'
#' @usage
#' woe.tree.binning(df, target.var, pred.var, min.perc.total,
#'                 min.perc.class, stop.limit, abbrev.fact.levels, event.class)
#'
#' @return
#' \code{woe.tree.binning} generates an object with the information necessary
#' for studying and applying the realized binning solution. When saved
#' it can be used with the functions \code{\link{woe.binning.plot}}, \code{\link{woe.binning.table}}
#' and \code{\link{woe.binning.deploy}}.
#'
#' @param df
#' Name of data frame with input data.
#' @param target.var
#' Name of dichotomous target variable in quotes. Only target variables with
#' two distinct values (e.g. 0, 1 or \dQuote{Y}, \dQuote{N}) are accepted;
#' cases with NAs in the target variable will be ignored.
#' @param pred.var
#' Name of predictor variable(s) to be binned in quotes.
#' A single variable name can be provided, e.g. \dQuote{varname1}, or a list of
#' variable names, e.g. c(\dQuote{varname1}, \dQuote{varname2}). Alternatively one
#' can repeat the name of the input data frame; the function will be applied
#' to all its variables apart from the target variable then.
#' Numeric variables and factors are supported and may contain NAs.
#' @param min.perc.total
#' For numeric variables this parameter defines the number of initial
#' classes before any merging or tree-like splitting is applied. For example
#' \emph{min.perc.total=0.05} (5\%) will result in 20 initial classes. For factors
#' the original levels with a percentage below this limit are collected in a
#' \sQuote{miscellaneous} level before the merging based on the \emph{min.perc.class}
#' and the tree-like splitting based on the WOE values starts. Increasing the
#' \emph{min.perc.total} parameter will avoid sparse bins. Accepted range: 0.0001-0.2;
#' default: 0.01.
#' @param min.perc.class
#' If a column percentage of one of the target classes within a bin is
#' below this limit (e.g. below 0.01=1\%) then the respective bin will be
#' joined with others. In case of numeric variables adjacent predictor classes
#' are merged. For factors respective levels (including sparse NAs) are
#' assigned to a \sQuote{miscellaneous} level. Setting \emph{min.perc.class}>0
#' may provide more reliable WOE values. Accepted range: 0-0.2;
#' default: 0, i.e. no merging with respect to sparse target classes
#' is applied.
#' @param stop.limit
#' Stops WOE based segmentation of the predictor's classes/levels in case the
#' resulting information value (IV) increases less than \emph{x}\% (e.g. 0.05 = 5\%)
#' compared to the preceding binning step. Increasing the \emph{stop.limit} will
#' simplify the binning solution and may avoid overfitting. Accepted range: 0-0.5;
#' default: 0.1.
#' @param abbrev.fact.levels
#' Abbreviates the names of new (merged) factor levels via the base R
#' \code{\link{abbreviate}} function in case the specified number of
#' characters is exceeded. Accepted range: 0-1000; default: 200.
#' 0 will prevent applying any abbreviation, i.e. only factor levels with
#' more than 1000 characters will be truncated then.
#' This option is particularly relevant in case one wants to generate dummy
#' variables via the \code{\link{woe.binning.deploy}} function, because the
#' factor levels will be part of the dummy variable names then.
#' @param event.class
#' Optional parameter for specifying the class of the target event. This
#' class typically indicates a negative event like a loan default or a
#' disease. Use integers (e.g. 1) or characters in quotes (e.g. \dQuote{bad}).
#' This class will be represented by negative WOE values then.
#' 
#' @family binning functions
#' 
#' @examples
#' # Load German credit data and create subset
#' data(germancredit)
#' df <- germancredit[, c('creditability', 'credit.amount', 'duration.in.month',
#'                   'savings.account.and.bonds', 'purpose')]
#'
#' # Bin a single numeric variable
#' binning <- woe.tree.binning(df, 'creditability', 'duration.in.month',
#'                            min.perc.total=0.01, min.perc.class=0.01,
#'                            stop.limit=0.1, event.class='bad')
#'
#' # Bin a single factor
#' binning <- woe.tree.binning(df, 'creditability', 'purpose',
#'                            min.perc.total=0.05, min.perc.class=0, stop.limit=0.1,
#'                            abbrev.fact.levels=50, event.class='bad')
#'
#' # Bin two variables (one numeric and one factor)
#' # with default parameter settings
#' binning <- woe.tree.binning(df, 'creditability', c('credit.amount','purpose'))
#'
#' # Bin all variables of the data frame (apart from the target variable)
#' # with default parameter settings
#' binning <- woe.tree.binning(df, 'creditability', df)
#'
#' @importFrom stats aggregate
#' @importFrom stats embed
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @export

##### This function calls the actual tree-like binning function above for every specified predictor variable that needs to be binned. #####

woe.tree.binning <- function(df, target.var, pred.var, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, event.class) {


	#### Warning message and defaults in case parameters are not specified
	if ( missing(df)==TRUE || missing(target.var)==TRUE || missing(pred.var)==TRUE ) { warning("Incorrect specification of data frame and/or variables.") }	
	if ( missing(min.perc.total)==TRUE ) { min.perc.total=0.01 }
	if ( min.perc.total<0.0001 || min.perc.total>0.2 || !is.numeric(min.perc.total) ) {
		warning("Incorrect parameter specification; accepted min.perc.total parameter range is 0.0001-0.2. Parameter was set to default (0.01).")
		min.perc.total=0.01
	}
	if ( missing(min.perc.class)==TRUE ) { min.perc.class=0 }
	if ( min.perc.class<0 || min.perc.class>0.2 || !is.numeric(min.perc.class) ) {
		warning("Incorrect parameter specification; accepted min.perc.class parameter range is 0-0.2. Parameter was set to default (0).")
		min.perc.class=0
	}
	if ( missing(stop.limit)==TRUE ) { stop.limit=0.1 }
	if ( stop.limit<0 || stop.limit>0.5 || !is.numeric(stop.limit) ) {
		warning("Incorrect parameter specification; accepted stop.limit parameter range is 0-0.05. Parameter was set to default (0.1).")
		stop.limit=0.1
	}
	if ( missing(abbrev.fact.levels)==TRUE ) { abbrev.fact.levels=200 }
	if ( abbrev.fact.levels<0 || abbrev.fact.levels>1000 ) {
		warning("Incorrect parameter specification; accepted abbrev.fact.levels parameter range is 0-10000. Parameter was set to default (200).")
		abbrev.fact.levels=200
	}

	#### Display warning message in case of incorrect target variable specification
	if ( !(length(unique(df[,target.var][!is.na(df[,target.var])]))==2) ) {
		warning("Incorrect variable specification.\nTarget variable must have two distinct values (NAs are accepted).")
	}

	#### Display warning message in case none of the target classes matches the specified event.class parameter
	if ( !missing(event.class) ) {
		if ( (unique(df[,target.var])[1]==event.class || unique(df[,target.var])[2]==event.class)==FALSE ) {
			warning("None of the target classes matches the specified event.class parameter.")
		}
	}
	
	#### In case bad class was specified assign 'good' and 'bad' codes (the latter will be associated with negative WOE values then)
	if ( !missing(event.class) ) { 
		if ( unique(df[,target.var])[1]==event.class ) {
			bad <- unique(df[,target.var])[1]
			good <- unique(df[,target.var])[2]
		} else {
			bad <- unique(df[,target.var])[2]
			good <- unique(df[,target.var])[1]
		}
	} else {
		bad <- unique(df[,target.var])[1]
		good <- unique(df[,target.var])[2]	
	}
	bad <- toString(bad)
	good <- toString(good)

	#### Gather names and look-up tables (with binned classes and WOE values) for each predictor variable in a list
	if ( is.data.frame(pred.var)==TRUE ) {
		pred.var <- as.list(colnames(subset(df, select=-c(which( colnames(df)==target.var )))))   # convert variable names of data frame into a list (without target variable)
	} else {
		as.list(pred.var)   # provide variable name(s) as a list
	}

	#### Subset: consider only cases without NA in target variable
	df <- df[!is.na(df[,target.var]),]
		
	#### Call actual binning function and put binning solutions together with respective variable names into a list
	binning <- lapply(pred.var, function(x) woe.tree.binning.2(df, target.var, x, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, bad, good))

	#### Read names and IV total values in the list and put them together with the binning tables
	names.of.pred.var <- lapply(pred.var, function(x) x)
	iv.total.list <- lapply(binning, function(x) colMeans(x[4]))
	binning <- matrix(c(names.of.pred.var, binning, iv.total.list),ncol=3)

	#### Sort via IV total
	binning <- binning[rev(sort.list(as.numeric(binning[,3]))),]
		
	binning

		
}