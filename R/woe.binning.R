##### This is the actual binning function for numeric variables and factors. #####

woe.binning.2 <- function(df, target.var, pred.var, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, bad, good) {


#### Build subsets with target and predictor variable
df <- df[, c(target.var, pred.var)]  # used for final binning
dfrm <- df[, c(target.var, pred.var)]   # used for iterative merging of bins
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


	## After sparse bins are merged:
	## Merge bins with similar WOE values and calculate corresponding WOE table and IV step by step
	## until 2 bins are left (i.e. 3 cutpoints: -Inf, middle cutpoint, +Inf)
	
	while ( length(cutpoints)>2 ) {
	
		# Compute binned variable from cutpoints and add it to the subset data frame
		dfrm$predictor.var.binned <- cut(dfrm$predictor.var, cutpoints, labels = NULL,
	    		include.lowest = FALSE, right = TRUE, dig.lab = 10,
			ordered_result = TRUE)
		
		# Compute crosstab from binned variable and target variable and covert it to a data frame
		freq.table <- table(dfrm$predictor.var.binned, dfrm$target.var, useNA="always")
		row.names(freq.table)[is.na(row.names(freq.table))] <- 'Missing'   # Replace NA in row.names with string 'Missing'
		woe.dfrm <- as.data.frame.matrix(freq.table)   # Convert frequency table to data frame
		woe.dfrm <- woe.dfrm[, c(good, bad)]   # Select columns with raw frequencies only
		
		# Compute WOE and information value (IV) from crosstab frequencies
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
		woe.dfrm$woe.lag <- c(NA, embed(woe.dfrm$woe,2)[,2])
		woe.dfrm$woe.diff <- abs(woe.dfrm$woe-woe.dfrm$woe.lag)
		woe.dfrm$iv.bins <- (woe.dfrm$col.perc.a-woe.dfrm$col.perc.b)*woe.dfrm$woe/100
					
		# Calculate total IV for current binning
		iv.total <- sum(woe.dfrm$iv.bins, na.rm=TRUE)

		# Collect total IVs for different binning solutions
		ifelse (exists('iv.total.collect', inherits=FALSE), iv.total.collect <- cbind(iv.total.collect, iv.total), iv.total.collect <- iv.total)
		
		# In case IV decreases by more than percentage specified by stop.limit parameter above
		# restore former binning solution (cutpoints) and leave loop
		if ( length(iv.total.collect)>1 ) {
			actual.iv.decrease <- ((iv.total.collect[length(iv.total.collect)-1]-iv.total.collect[length(iv.total.collect)])/(iv.total.collect[length(iv.total.collect)-1]))

			if ( (actual.iv.decrease>stop.limit) && (exists('stop.limit.exceeded', inherits=FALSE)==FALSE) ) {
				cutpoints.final <- cutpoints.backup
				woe.dfrm.final <- woe.dfrm.backup
				stop.limit.exceeded <- TRUE   # indicates that stop limit is exceeded to prevent overriding the final solution
			}
		}
		
		# Save first cutpoint solution and corresponding WOE values as final solution (is used in case no WOE merging will be applied)
		if ( exists('cutpoints.backup', inherits=FALSE)==FALSE ) {
			cutpoints.final <- cutpoints
			woe.dfrm.final <- woe.dfrm
		}

		# Saves binning solution after last merging step in case the IV stop limit was not exceeded
		if ( (exists('stop.limit.exceeded', inherits=FALSE)==FALSE) && (length(cutpoints)==3) ) {
			cutpoints.final <- cutpoints
			woe.dfrm.final <- woe.dfrm
		}
		
		# Save backups of current cutpoints and corresponding WOE values before merging to be able to retrieve solution in case IV decrease is too strong
		cutpoints.backup <- cutpoints
		woe.dfrm.backup <- woe.dfrm

		# Determine the index of the minimum WOE difference between adjacent bins and
		# merge bins with minimum WOE difference (apart from the last 'Missing' bin)	
		min.woe.diff <- which(woe.dfrm$woe.diff[-nrow(woe.dfrm)]==min(woe.dfrm$woe.diff[-nrow(woe.dfrm)], na.rm=TRUE))
		cutpoints <- cutpoints[-c(min.woe.diff)]
		
	}


	## Compute final IV
	iv.total.final <- sum(woe.dfrm.final$iv.bins, na.rm=TRUE)
		
	## Save final binning solution via look-up-table for deployment
	lower.cutpoints.final.dfrm <- as.data.frame(cutpoints.final)
	upper.cutpoints.final.dfrm <- rbind(as.data.frame(cutpoints.final[-1]),'Missing')
	look.up.table <- cbind(woe.dfrm.final[, 5,drop=FALSE], lower.cutpoints.final.dfrm, upper.cutpoints.final.dfrm)
#	if ( look.up.table[nrow(look.up.table),1]==0 ) { look.up.table[nrow(look.up.table),1] <- NA }   # replace WOE=0 in Missing row with NA (because this only occurs in case Missing Data does not occur during binning)
	look.up.table <- cbind.data.frame(look.up.table, iv.total.final, woe.dfrm.final[, c(1,2,3,4,8),drop=FALSE])   # add column with final total Information Value


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
	## Merge levels with similar WOE values and calculate corresponding WOE table and IV step by step until
	## 2 regular bins (+ Missing or 'misc. level') are left

	while ( length(levels(dfrm$predictor.var.binned))>3 ) {
	
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
		woe.dfrm$woe.lag <- c(NA, embed(woe.dfrm$woe,2)[,2])
		woe.dfrm$woe.diff <- abs(woe.dfrm$woe-woe.dfrm$woe.lag)
		woe.dfrm$iv.bins <- (woe.dfrm$col.perc.a-woe.dfrm$col.perc.b)*woe.dfrm$woe/100
		
		# Calculate total IV for current binning
		iv.total <- sum(woe.dfrm$iv.bins, na.rm=TRUE)
		
		# Collect total IVs for different binning solutions
		ifelse (exists('iv.total.collect', inherits=FALSE), iv.total.collect <- cbind(iv.total.collect, iv.total), iv.total.collect <- iv.total)
		
		# In case IV decreases by more than percentage specified by stop.limit parameter above
		# restore former binning solution (cutpoints) and leave loop
		if ( length(iv.total.collect)>1 ) {
			actual.iv.decrease <- ((iv.total.collect[length(iv.total.collect)-1]-iv.total.collect[length(iv.total.collect)])/(iv.total.collect[length(iv.total.collect)-1]))
			if ( (actual.iv.decrease>stop.limit) && (exists('stop.limit.exceeded', inherits=FALSE)==FALSE) ) {
				stop.limit.exceeded <- TRUE   # indicates that stop limit is exceeded to prevent overriding the final solution
			}
		}
		
		# Merge until 2 regular bins remain
		
		if ( length(levels(dfrm$predictor.var.binned))>3 ) {
			
			# Merge levels with most similar WOE values
			min.woe.diff <- which(woe.dfrm$woe.diff==min(woe.dfrm$woe.diff, na.rm=TRUE))
			levels(dfrm$predictor.var.binned)[levels(dfrm$predictor.var.binned)%in%c(row.names(woe.dfrm)[min.woe.diff][[1]][1],row.names(woe.dfrm)[min.woe.diff-1][[1]][1])] <- paste(row.names(woe.dfrm)[min.woe.diff][[1]][1], "+", row.names(woe.dfrm)[min.woe.diff-1][[1]][1])
			
			# Save names of the factor levels that are merged
			list.level.a <- as.list(row.names(woe.dfrm)[min.woe.diff][[1]][1])
			list.level.b <- as.list(row.names(woe.dfrm)[min.woe.diff-1][[1]][1])
			
			# Collect names of the factor levels that are merged in lists (until stop criteria is reached)
			if ( exists('list.level.a.collected', inherits=FALSE)==FALSE ) {
				list.level.a.collected <- list.level.a
				list.level.b.collected <- list.level.b
			}
			else {
				if ( exists('stop.limit.exceeded', inherits=FALSE)==FALSE ) {
					list.level.a.collected <- c(list.level.a.collected, list.level.a)
					list.level.b.collected <- c(list.level.b.collected, list.level.b)
				}
				else {
					list.level.a.collected <- list.level.a.collected[1:length(list.level.a.collected)-1]
					list.level.b.collected <- list.level.b.collected[1:length(list.level.b.collected)-1]
				}
			}
	
		}

	
	}
	
	### Apply FINAL binning to INPUT data
	
	## Merge factor levels
	# Merge sparse levels
	levels(df[,ncol(df)])[levels(df[,ncol(df)])%in%(row.names(woe.dfrm.sparse.subset.pos))] <- "misc. level pos."
	levels(df[,ncol(df)])[levels(df[,ncol(df)])%in%(row.names(woe.dfrm.sparse.subset.neg))] <- "misc. level neg."
	# Merge levels with similar WOE values
	if ( exists('list.level.a.collected', inherits=FALSE)==TRUE ) {
		for ( i in 1:length(list.level.a.collected) ) {
			levels(df[,ncol(df)])[levels(df[,ncol(df)])%in%c(list.level.a.collected[i],list.level.b.collected[i])] <- paste(list.level.a.collected[i], "+", list.level.b.collected[i])
		}
	}

	## Repeat generating WOE table for selected binning solution
	
	# Compute crosstab from binned variable and target variable and covert it to a data frame
	freq.table.final <- table(df[,ncol(df)], dfrm$target.var)
	#row.names.final(freq.table.final)[is.na(row.names(freq.tabl.finale))] <- 'Missing'   # replace NA in row.names with string 'Missing'
	woe.dfrm.final <- as.data.frame.matrix(freq.table.final)   # convert frequency table to data frame
	woe.dfrm.final <- woe.dfrm.final[, c(good, bad)]   # Select columns with raw frequencies only

	# Compute WOE and information value (IV) from crosstab frequencies
	woe.dfrm.final$col.perc.a <- woe.dfrm.final[,1]/sum(woe.dfrm.final[,1])
	woe.dfrm.final$col.perc.b <- woe.dfrm.final[,2]/sum(woe.dfrm.final[,2])
	# Correct column percents in case of 0 frequencies
	if ( min(woe.dfrm.final[,1],na.rm=TRUE)==0 | min(woe.dfrm.final[,2],na.rm=TRUE)==0 ) {
		woe.dfrm.final$col.perc.a <- (woe.dfrm.final$col.perc.a+0.0001)/sum(woe.dfrm.final$col.perc.a+0.0001)
		woe.dfrm.final$col.perc.b <- (woe.dfrm.final$col.perc.b+0.0001)/sum(woe.dfrm.final$col.perc.b+0.0001)	
	}	
	woe.dfrm.final$woe <- 100*log(woe.dfrm.final$col.perc.a/woe.dfrm.final$col.perc.b)
	woe.dfrm.final$woe[is.finite(woe.dfrm.final$woe)==FALSE] <- NA   # convert Inf, -Inf and NaN to NA
	woe.dfrm.final <- woe.dfrm.final[order(woe.dfrm.final$woe),]   # sort data via WOE values
	woe.dfrm.final$iv.bins <- (woe.dfrm.final$col.perc.a-woe.dfrm.final$col.perc.b)*woe.dfrm.final$woe/100	
	iv.total.final <- sum(woe.dfrm.final$iv.bins, na.rm=TRUE)


	## Add variable with corresponding WOE values for final binning
	
	# Save row order of input data as ID variable
	df$initial.order.id  <- 1:nrow(df)
	
	# Add final binned (numerical) variable with WOE values (via left join with WOE table)
	df <- merge(df, woe.dfrm.final[,5,drop=FALSE], by.x=colnames(df)[ncol(df)-1], by.y=0, all.x=TRUE)
	colnames(df)[ncol(df)] <- paste(pred.var,".binned.woe",sep="")
	
	# Restore initial column and row order and get rid of initial.order.id and row names
	df <- df[order(df$initial.order.id), ]
	df <- subset(df, select=c(2:(ncol(df)-2),1,ncol(df)))
	row.names(df) <- NULL


	## Save final binning solution via look-up-table for deployment
	levels(df[,pred.var]) <- c(levels(df[,pred.var]), "Missing")   # add factor level 'Missing'
	df[,pred.var][is.na(df[,pred.var])] <- "Missing"   # replace NA with string 'Missing'
	look.up.table <- aggregate(df[,ncol(df)], list(df[,pred.var], df[,ncol(df)-1]), mean, na.rm=TRUE)
	look.up.table <- cbind.data.frame(look.up.table, iv.total.final)   # add column with final total Information Value
	colnames(look.up.table)[3] <- "woe"
	look.up.table <- merge(look.up.table, woe.dfrm.final[ , -which(names(woe.dfrm.final) %in% c("woe"))], by.x=2, by.y=0)
	look.up.table <- look.up.table[order(look.up.table$woe, look.up.table$Group.2),]   # sort by woe value and merged bin name

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



#' @title Binning via Fine and Coarse Classing
#'
#' @description
#' \code{woe.binning} generates a supervised fine and coarse classing of numeric
#' variables and factors with respect to a dichotomous target variable. Its parameters
#' provide flexibility in finding a binning that fits specific data characteristics
#' and practical needs.
#'
#' @section Binning of Numeric Variables:
#' Numeric variables (continuous and ordinal) are binned by merging initial classes with
#' similar frequencies. The number of initial bins results from the \emph{min.perc.total}
#' parameter: min.perc.total will result in trunc(1/min.perc.total) initial bins,
#' whereby \emph{trunc} is needed to guarantee bins with similar frequencies.
#' For example \emph{min.perc.total=0.07} will cause trunc(14.3)=14 initial classes.
#' Next, if \emph{min.perc.class}>0, bins with sparse target classes will be merged with
#' the next upper bin, and in case of the last bin with the next lower one. NAs have
#' their own bin and will not be merged with others. Finally nearby bins with most similar
#' weight of evidence (WOE) values are joined step by step until the information value
#' (IV) decreases more than specified by a percentage value (\emph{stop.limit} parameter)
#' or until two bins are reached.
#' @section Binning of Factors:
#' Factors (categorical variables) are binned by merging factor levels. As a start sparse
#' levels (defined via the \emph{min.perc.total} and \emph{min.perc.class} parameters)
#' are merged to a \sQuote{miscellaneous} level: if possible, respective levels (including
#' sparse NAs) are bundled as \sQuote{misc. level pos.} (associated with positive WOE
#' values), respectively as \sQuote{misc. level neg.} (associated with negative WOE
#' values). In case a misc. level contains only NAs it will be named \sQuote{Missing}.
#' Afterwards levels with similar WOE values are joined step by step until the information
#' value (IV) decreases more than specified by a percentage value (\emph{stop.limit} parameter)
#' or until two bins are reached.
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
#' woe.binning(df, target.var, pred.var, min.perc.total,
#'             min.perc.class, stop.limit, abbrev.fact.levels, event.class)
#'
#' @return
#' \code{woe.binning} generates an object containing the information necessary
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
#' classes before any merging is applied. For example \emph{min.perc.total=0.05}
#' (5\%) will result in 20 initial classes. For factors the original
#' levels with a percentage below this limit are collected in a \sQuote{miscellaneous}
#' level before the merging based on the \emph{min.perc.class} and on the
#' WOE starts. Increasing the \emph{min.perc.total} parameter will avoid
#' sparse bins. Accepted range: 0.01-0.2; default: 0.05.
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
#' Stops WOE based merging of the predictor's classes/levels in case the
#' resulting information value (IV) decreases more than \emph{x}\% (e.g. 0.05 = 5\%)
#' compared to the preceding binning step. \emph{stop.limit=0} will skip any
#' WOE based merging. Increasing the \emph{stop.limit} will simplify the binning
#' solution and may avoid overfitting. Accepted range: 0-0.5; default: 0.1.
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
#' binning <- woe.binning(df, 'creditability', 'duration.in.month',
#'                        min.perc.total=0.05, min.perc.class=0.01,
#'                        stop.limit=0.1, event.class='bad')
#'
#' # Bin a single factor
#' binning <- woe.binning(df, 'creditability', 'purpose',
#'                        min.perc.total=0.05, min.perc.class=0, stop.limit=0.1,
#'                        abbrev.fact.levels=50, event.class='bad')
#'
#' # Bin two variables (one numeric and one factor)
#' # with default parameter settings
#' binning <- woe.binning(df, 'creditability', c('credit.amount','purpose'))
#'
#' # Bin all variables of the data frame (apart from the target variable)
#' # with default parameter settings
#' binning <- woe.binning(df, 'creditability', df)
#'
#' @importFrom stats aggregate
#' @importFrom stats embed
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @export

##### This function calls the actual binning function above for every specified predictor variable that needs to be binned. #####

woe.binning <- function(df, target.var, pred.var, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, event.class) {


	#### Warning message and defaults in case parameters are not specified
	if ( missing(df)==TRUE || missing(target.var)==TRUE || missing(pred.var)==TRUE ) { warning("Incorrect specification of data frame and/or variables.") }	
	if ( missing(min.perc.total)==TRUE ) { min.perc.total=0.05 }
	if ( min.perc.total<0.01 || min.perc.total>0.2 || !is.numeric(min.perc.total) ) {
		warning("Incorrect parameter specification; accepted min.perc.total parameter range is 0.01-0.2. Parameter was set to default (0.05).")
		min.perc.total=0.05
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
	binning <- lapply(pred.var, function(x) woe.binning.2(df, target.var, x, min.perc.total, min.perc.class, stop.limit, abbrev.fact.levels, bad, good))

	#### Read names and IV total values in the list and put them together with the binning tables
	names.of.pred.var <- lapply(pred.var, function(x) x)
	iv.total.list <- lapply(binning, function(x) colMeans(x[4]))
	binning <- matrix(c(names.of.pred.var, binning, iv.total.list),ncol=3)

	#### Sort via IV total
	binning <- binning[rev(sort.list(as.numeric(binning[,3]))),]
		
	binning

		
}