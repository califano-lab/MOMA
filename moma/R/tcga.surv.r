
#' @title get.clin Parse the full clinical merged file from GDAC-Firehose and make a dlpyr tibble out of it
#' @import tibble
#' @importFrom dplyr inner_join select contains 
#' @importFrom readr read_tsv
#' @import survival
#' @import grDevices
#' @import graphics
#' @import RColorBrewer
#' @param clin.file Clinical Merge file from TCGA
#' @return Returns a tibble data.frame of the clinical data. Rows are patient.bcr_patient_barcode IDs
#' @export
get.clin <- function(clin.file=NULL) {

	if (is.null(clin.file)) {
		clin.dir <- list.files(pattern="gdac..*Merge.*")
		if (length(clin.dir)==0) { print ("Error: can't find clinical data directory!"); q(); }

		tissue <- strsplit(strsplit(clin.dir, '_')[[1]][2], '\\.')[[1]][1]
		clin.file <- paste0(clin.dir, '/', tissue, '.clin.merged.txt')
	}

	data <- t(readr::read_tsv(clin.file, col_names=FALSE))
	patient.barcodes <- as.character(data[,which(data[1,]=="patient.bcr_patient_barcode")])
	rownames(data) <- patient.barcodes
	colnames(data) <- data[1,]
	clinical <- tibble::as_tibble(data)

	clinical
}


#' @title Add clustering solutions to the tibble dataset 'survObj' and return
#' @param data Clinical survival data in tibble format, processed from TCGA GDAC / Firehose
#' @param clustering Clustering solution with numeric values as cluster ID, named by samples
#' @return The same tibble with an additional 'clusters' variable added
#' @export
tibble.add_clusters <- function(data, clustering) {

	# Make a joint vector
	data$samples <- toupper(data$patient.bcr_patient_barcode)

	clusters <- data.frame(samples=names(clustering), cluster=clustering)
	res <- dplyr::inner_join(data, clusters, by=c("samples"="samples"))
	if (length(res$cluster)==0) {
		# likely that clustering names are not patient ids but full sample ids: remove the last -01 and it should work
		names(clustering) <- unlist(lapply(strsplit(names(clustering), '-'), function(x) paste(x[1:3], collapse='-')))
		clusters <- data.frame(samples=names(clustering), cluster=clustering)
		res <- dplyr::inner_join(data, clusters, by=c("samples"="samples"))
		if (length(res$cluster)==0) {
			print ("Failed to add clustering solution in tibble.add_clusters, sample names didn't match!")
			q();
		}
	}
	res
}


#' @title Add a survival object to the tibble dataset 'survObj' and return
#' @param data Clinical survival data in tibble format, processed from TCGA GDAC / Firehose
#' @return The original tibble with a 'survObj' added to the data frame
#' @export
#' @return A copy of the input tibble with a '$survObj' attached
tibble.survfit <- function(data) {

	data$days.death <- as.numeric(apply(as.matrix(dplyr::select(data, dplyr::contains("days_to_death"))), 1, function(x) max(na.omit(x))))
	data$days.last_followup <- as.numeric(apply(as.matrix(dplyr::select(data, dplyr::contains("days_to_last_followup"))), 1, function(x) max(na.omit(x))))
	# vital status summary
	data$vital_stat_collapsed <- unlist(apply(as.matrix(dplyr::select(data, dplyr::contains("vital_status"))), 1, function(x) {

		x <- na.omit(x)
		if (length(x)==0) {
			status <- NA		
			return (status)
		}

		status <- NA
		if (any(x=="dead")) {
			status <- "dead"
		} else if (all(x=="alive")) {
			status <- "alive"
		}
		status

	}))

	# there may be multiple intersections in TCGA data: it's buggy!! Take the largest number here, let it overwrite the last 
	data$days.last_observed <- as.numeric(apply(cbind(data$days.death, data$days.last_followup), 1, function(x) ifelse(is.na(x[1]), x[2], x[1]) ))

	# Make a joint vector
	data$vital_status.binary <- rep(0, length(data$patient.vital_status))
	data$vital_status.binary[data$vital_stat_collapsed=="dead"] <- 1

	# built in check: do all those with a days to death date also have a vital status of 1??
	if (!all(data$vital_status.binary[which(!is.na(data$days.death))]==1)) {
		print ("Error: patient with death days column has vital status of 'alive'!")
		idx <- which(!is.na(data$days.death))
		print (cbind(data$samples[idx], data$vital_status.binary[idx], data$days.death[idx], data$days.last_followup[idx]))
		q();
	}
	
	data$survObj <- survival::Surv(data$days.last_observed, data$vital_status.binary)

	data
}

#' @title progression free survival
#' @param data Clinical survival data in tibble format, processed from TCGA GDAC / Firehose
#' @return A copy of the input tibble with a $survObj attached
#' @export
tibble.survfit.progression_free <- function(data) {

	data$days.death <- as.numeric(apply(as.matrix(dplyr::select(data, dplyr::contains("days_to_death"))), 1, function(x) max(na.omit(x))))
	data$days.last_followup <- as.numeric(apply(as.matrix(dplyr::select(data, dplyr::contains("days_to_last_followup"))), 1, function(x) max(na.omit(x))))
	data$days.new_tumor <- as.numeric(apply(as.matrix(dplyr::select(data, dplyr::contains("days_to_new_tumor_event"))), 1, function(x) max(na.omit(x))))

	# vital status summary
	data$vital_stat_collapsed <- unlist(apply(as.matrix(dplyr::select(data, dplyr::contains("vital_status"))), 1, function(x) {

		x <- na.omit(x)
		if (length(x)==0) {
			status <- NA		
			return (status)
		}

		status <- NA
		if (any(x=="dead")) {
			status <- "dead"
		} else if (all(x=="alive")) {
			status <- "alive"
		}
		status

	}))

	# there may be multiple intersections in TCGA data: it's buggy!! Take the days to death, or days to last followup if no death date
	data$days.last_remission_free <- as.numeric(apply(cbind(data$days.death, data$days.new_tumor, data$days.last_followup), 1, function(x) {
		if(!all(is.na(x[1:2]))) {
			# the smaller of days to death or remission
			return (min(x[1:2]))
		} else {
			return (x[3])
		}
	}))

	# Make a joint vector
	data$remission_status_binary <- rep(0, length(data$patient.vital_status))
	data$remission_status_binary[data$vital_stat_collapsed=="dead"] <- 1
	# set remission status to 1 if any new tumor found
	data$remission_status_binary[which(!is.na(data$days.new_tumor))] <- 1
	
	# built in check: do all those with a days to death date also have a vital status of 1??
	if (!all(data$vital_status.binary[which(!is.na(data$days.death))]==1)) {
		print ("Error: patient with death days column has vital status of 'alive'!")
		idx <- which(!is.na(data$days.death))
		print (cbind(data$samples[idx], data$vital_status.binary[idx], data$days.death[idx], data$days.last_followup[idx]))
		q();
	}
	
	data$survObj <- survival::Surv(data$days.last_remission_free, data$remission_status_binary)

	data
}


#' @title Get survival based on clustering
#' @param clustering Clustering solution with numeric values as cluster ID, named by samples
#' @param clinical.tibble Clinical data from TCGA GDAC / Firehose in tibble format
#' @param progression.free.surv Compute progression-free survival only? (default=FALSE)
#' @return A list containing the survObj, and associated p-values between the best-worst surviving clustrers, as well as overall stratification based on a COX PH model.
#' @export
tibble.survfit_select <- function(clustering, clinical.tibble, progression.free.surv=FALSE) {

	names(clustering) <- unlist(lapply(strsplit(names(clustering), '-'), function(x) paste(x[1:3], collapse='-')))
	clin.tibble <- tibble.add_clusters(clinical.tibble, clustering)

	if (!progression.free.surv) {
		clin.tibble <- tibble.survfit(clin.tibble)
	} else {
		clin.tibble <- tibble.survfit.progression_free(clin.tibble)
	}
	# Find the best and worst surviving clusters based on the statistics
	SurvDiff <- survival::survdiff(survObj ~ cluster, data=clin.tibble, rho=0)
	pval.overall <- 1 - pchisq(SurvDiff$chisq, length(SurvDiff$n) - 1)

	survival.scores <- SurvDiff$obs/SurvDiff$exp
	best.surv.cluster <- strsplit(names(SurvDiff$n)[survival.scores==min(survival.scores)], '=')[[1]][2]
	worst.surv.cluster <- strsplit(names(SurvDiff$n)[survival.scores==max(survival.scores)], '=')[[1]][2]

	# Redo the survival analysis difference between just these clusters
	clin.subset <- dplyr::filter(clin.tibble, cluster %in% as.numeric(c(best.surv.cluster, worst.surv.cluster)))
	if (!progression.free.surv) {
		clin.subset <- tibble.survfit(clin.subset)
	} else {
		clin.subset <- tibble.survfit.progression_free(clin.subset)
	}
	SurvDiff.best_worst <- survival::survdiff(survObj ~ cluster, data=clin.subset, rho=0)
	pval.best_worst <- 1 - pchisq(SurvDiff.best_worst$chisq, length(SurvDiff$n) - 1)
	list(survObj = clin.tibble$survObj, pval.best_v_worst = pval.best_worst, pval.overall = pval.overall, 
		best.clus = best.surv.cluster, worst.clus = worst.surv.cluster)
}


#' @title Get the set of optimal clusterings from the Silhouette results. Find the one that optimizes survival
#' @param search.results Multi-k clustering solution, including analytical cluster silhouette scores
#' @param clinical.tibble Clinical data from TCGA GDAC / Firehose in tibble format
#' @param tissue TCGA tissue type
#' @param progression.free.surv Compute progression-free survival only? (default=FALSE)
#' @return A list containing the best clustering solution and associated survival separation p-value. Also includes the best and worst surviving clusters
#' @export
get.best.clustering.supervised <- function(search.results, clinical.tibble, tissue, progression.free.surv=FALSE) {

	search.results

	if (progression.free.surv) {
		print ("Fitting as progression/disease - free survival...")
	}

	scores <- unlist(lapply(search.results, function(x) x$reliability))
	names(scores) <- unlist(lapply(search.results, function(x) x$k))
	ks <- unlist(lapply(search.results, function(x) x$k))
	clustering <- search.results[[which(scores==max(scores))]]
	print (paste("Best clustering k = ", clustering$k))

	# compare to get set of equivalent scores
	pvals <- unlist(lapply(2:15, function(k) {
		cs <- search.results[[(k-1)]]
		pval.1 <- ks.test(cs$element.reliability, clustering$element.reliability)$p.value
		pval.1
	}))
	names(pvals) <- 2:15
	fwers <- sort(p.adjust(pvals, method='bonferroni'))
	# those where we can't reject the null hypothesis are considered equivalent
	equiv.clusters <- as.numeric(names(fwers[which(fwers > 0.05)]))
	print (paste("Considering statistically equivalent clustering solutions (analytical score) ", equiv.clusters))

	# for each clustering, do survival analysis	
	survival.results <- lapply(equiv.clusters, function(k) {
		#print (paste("Calculating survival separation on clustering k = ", k))
		clustering <- search.results[[(k-1)]]$clustering
		res <- tibble.survfit_select(clustering, clinical.tibble, progression.free.surv)
		res
	})
	names(survival.results) <- equiv.clusters
	pvals.overall <- sort(sapply(survival.results, function(x) x$pval.overall), decreasing=FALSE)
	best.k <- as.numeric(names(pvals.overall[1]))

	pvals.best_v_worst <- sort(sapply(survival.results, function(x) x$pval.best_v_worst), decreasing=FALSE)

	#add clustering to 'alt.clusters' variable in tibble

	clustering <- search.results[[(best.k-1)]]$clustering
	surv.res <- survival.results[[as.character(best.k)]]
	# surv.res$best.clus
	# surv.res$worst.clus
	# redo surv objets
	
	survplot.best.v.worst(clinical.tibble, clustering, 'survival.supervised.png', paste(toupper(tissue), " supervised clustering, best vs worst survival clusters"))

	list(clustering=clustering, clustering.k=best.k, pval.overall=pvals.overall[1], pval.best_v_worst=pvals.best_v_worst[1], clus.choices=list('best'=surv.res$best.clus, 'worst'=surv.res$worst.clus))
}


#' @title Get the best and worst clustering solutions, plot the Kaplan-Meyer of those populations
#' @param clinical.tibble Clinical data from TCGA GDAC / Firehose in tibble format
#' @param clustering Clustering solution with numeric values as cluster ID, named by samples
#' @param output Output .png file for Kaplan-Meyer plot
#' @param title.print Plot title
#' @return Nothing, writes the plot to the supplied .png output file
#' @export
survplot.best.v.worst <- function(clinical.tibble, clustering, output, title.print) {
	
	names(clustering) <- unlist(lapply(strsplit(names(clustering), '-'), function(x) paste(x[1:3], collapse='-')))
	clin.tibble <- tibble.add_clusters(clinical.tibble, clustering)
	clin.tibble <- tibble.survfit(clin.tibble)
	# Find the best and worst surviving clusters based on the statistics
	SurvDiff <- survival::survdiff(survObj ~ cluster, data=clin.tibble, rho=0)
	#print (SurvDiff)
	survival.scores <- SurvDiff$obs/SurvDiff$exp
	best.surv.cluster <- strsplit(names(SurvDiff$n)[survival.scores==min(survival.scores)], '=')[[1]][2]
	worst.surv.cluster <- strsplit(names(SurvDiff$n)[survival.scores==max(survival.scores)], '=')[[1]][2]

	# Redo the survival analysis difference between just these clusters
	clin.subset <- dplyr::filter(clin.tibble, cluster %in% as.numeric(c(best.surv.cluster, worst.surv.cluster)))
	clin.subset <- tibble.survfit(clin.subset)
	SurvDiff <- survival::survdiff(survObj ~ cluster, data=clin.subset, rho=0)
	pval <- 1 - pchisq(SurvDiff$chisq, length(SurvDiff$n) - 1)

	#import::from(survminer, ggsurvplot)

	fit <- survival::survfit(survObj ~ cluster, data=clin.subset)
	colors <- colorRampPalette(brewer.pal(8,"Dark2"))(2)
	png(output, width=700, height=500)
	par(mar=c(c(10,10,10,10)))
	plot(fit,
		conf.type="log-log", 
		col=colors,
		lty=1:2, 
		lwd=4,
		mark.time=TRUE,
		cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2,
		xlab="Survival Time (Days)", ylab="Survival Probability")
#
	pval.print <- formatC(pval, format="e", digits=1)
	#title(main=title.print)
#
	mtext(pval.print, pos=2)
	legend.labels <- unique(clin.subset$cluster)
	legend("top", legend=legend.labels, col=colors, lty=1:2, lwd=3, horiz=FALSE)
	dev.off()
}


