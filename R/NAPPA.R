# The NAPPA package is copyright (C) 2013 AztraZeneca
# This program is free software; you can redistribute it and/pr modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License or any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details
# You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# Author(s): Mark Wappett, Mike Dymond, Tom Liptrot and Chris Harbron

NAPPA <- function(data=data, tissueType=tissueType, sampleNumber=sampleNumber)
{
	temprow <- matrix(c(rep.int(NA,length(data))),nrow=2,ncol=length(data))
	newrow <- data.frame(temprow)
	colnames(newrow) <- colnames(data)
	data <- rbind(newrow, data)
	lengthHK <- length(which(data[,1] == "Housekeeping"))
	if(lengthHK == 0)
	{
		err <- "No genes have been indicated as housekeepers"
	}

	# Report problems with positive controls and return error if below 5000 on average
	pcTEST <- as.numeric(data[22,4:dim(data)[2]])
	w1 <- which(pcTEST < 5000)
	if(length(w1) > 0)
	{
		lowINDICES <- paste(w1, collapse=", ")
		err <- paste(c("Lanes ", lowINDICES, " have counts of < 5000 - this may indicate a problem with the codeset for that sample.  You may want to remove these samples"), collapse=" ")
	}

	# Prepare matric for testing housekeeper performance
	w1 <- which(data[,1] == "Housekeeping")
	hkCOUNTtest <- data[w1,]
	hkCOUNTtestVEC <- numeric(0)
	iTH <- 0
	for(i in 4:dim(data)[2])
	{
		iTH <- iTH + 1
		hkCOUNTtest[,i] <- as.numeric(hkCOUNTtest[,i])
		m1 <- mean(hkCOUNTtest[,i])
		hkCOUNTtestVEC[iTH] <- m1
	}

	# Report problems with housekeepers and return error if below 50 on average
	w1 <- which(hkCOUNTtestVEC < 50)
	if(length(w1) > 0)
	{
		lowINDICES <- paste(w1, collapse=", ")
		err <- paste(c("Housekeeping genes in lane ", lowINDICES, " are less than 50 on average.  This may indicate a problem with the sample amount/ codeset. The housekeeper normalisation step is likely to fail.  You may want to remove these samples"), collapse=" ")
		print(err)
	}
	n2 <- as.character(t(data[4,])[4:dim(data)[2]])
	n1 <- c("Code.Class", "Name", "Accession")
	n <- c(n1, n2)
	all <- data[22:dim(data)[1],]
	colnames(all) <- n
	all2 <- all
	for(i in 4:dim(all)[2])
	{
		all2[,i] <- as.numeric(all[,i]) / as.numeric(data[14,i])
	}

	# Multiply back by 280
	for(i in 4:dim(all2)[2])
	{
		all2[,i] <- all2[,i]*280
	}
	whichNeg <- which(all2[,1] == "Negative")
	negative <- subset(all2, all2[,1] == "Negative")
	all2 <- all2[-whichNeg,]
	iTH <- 0
	negIND <- numeric(0)
	for(i in 4:dim(negative)[2])
	{
		iTH <- iTH + 1
		m1 <- mean(negative[,i])
		negIND[iTH] <- m1
	}

	# Seperate out genes of different classes
	pos1 <- subset(all2, all2[,1] == "Positive")
	hk1 <- subset(all2, all2[,1] == "Housekeeping")
	end1 <- subset(all2, all2[,1] == "Endogenous")
	columnNames <- colnames(end1)[4:dim(end1)[2]]

	# By Lane, characterise and remove 0 values, and keep expression data in seperate indexed lists of data frames
	zeroIndexList <- list()
	nnIndexList <- list()
	nnGeneList <- list()
	zeroGeneList <- list()
	iTH <- 0
	for(i in 4:dim(end1)[2])
	{
		iTH <- iTH + 1
		w1 <- which(end1[,i] <= 0)
		w2 <- which(end1[,i] > 0)
		zeroIndexList[[iTH]] <- w1
		nnIndexList[[iTH]] <- w2
		t1 <- as.data.frame(end1[,i], stringsAsFactors=FALSE)
		if(length(w1) > 0)
		{
			t2 <- as.data.frame(t1[w1,], stringsAsFactors=FALSE)
			t1 <- t1[-w1,]
			t1 <- as.data.frame(t1, stringsAsFactors=FALSE)
			metaData <- end1[-w1,1:3]
			metaData2 <- end1[w1,1:3]
			mT <- cbind(metaData, t1,w2)
			mT2 <- cbind(metaData2, t2,w1)
			nnGeneList[[iTH]] <- mT
			zeroGeneList[[iTH]] <- mT2
			names(nnGeneList[[iTH]]) <- c("Code.Class", "Name", "Accession", "Expression")
			nnGeneList[[iTH]][,4] <- as.numeric(nnGeneList[[iTH]][,4])
			rownames(nnGeneList[[iTH]]) <- nnGeneList[[iTH]][,2]
		}
		else
		{
			metaData <- end1[,1:3]
			mT <- cbind(metaData, t1,w2) 
			nnGeneList[[iTH]] <- mT
			names(nnGeneList[[iTH]]) <- c("Code.Class", "Name", "Accession", "Expression")
			nnGeneList[[iTH]][,4] <- as.numeric(nnGeneList[[iTH]][,4])
			rownames(nnGeneList[[iTH]]) <- nnGeneList[[iTH]][,2]
		}
	}
	
		subtract.bg <- function(count, b){
	if(count <=20 | count<=2*b) (count - b*(1-(b^count)/factorial(count)/sum(b^(0:count)/factorial(count))))
	else(count - b)
	}

	# Subtract background from endogenous genes
	iTH <- 0
	nnGeneList2 <- nnGeneList
	for(i in 1:length(nnGeneList))
	{
		iTH <- iTH + 1
		table1 <- nnGeneList[[i]]
		for(j in 1:dim(table1)[1])
		{
			table1[j,4] <- subtract.bg(table1[j,4], negIND[iTH])
		}
		nnGeneList2[[i]] <- table1
	}

	# Subtract background from housekeeping genes
	iTH <- 0
	for(i in 4:dim(hk1)[2])
	{
		iTH <- iTH + 1
		for(j in 1:dim(hk1)[1])
		{
			hk1[j,i] <- subtract.bg(hk1[j,i], negIND[iTH])
		}
	}
	# For each lane, calculate the E2 value, and divide all the Endogenous genes by this value (including housekeepers)
	E2vec <- numeric(0)
	iTH <- 3
	nnGeneList3 <- nnGeneList2
	for(i in 1:length(nnGeneList2))
	{
		iTH <- iTH + 1
		table1 <- nnGeneList2[[i]]
		E2 <- (sum(pos1[1:4,iTH]) / (128+32+8+2))
		E2vec[iTH] <- E2
		table1[,4] <- table1[,4] / E2
		hk1[,iTH] <- hk1[,iTH] / E2
		table1[,4] <- log2(table1[,4])
		hk1[,iTH] <- log2(hk1[,iTH])
		nnGeneList3[[i]] <- table1
	}

	# Retrieve Housekeeper genes
	hkVEC <- numeric(0)
	iTH <- 0
	for(i in 4:dim(hk1)[2])
	{
		iTH <- iTH + 1
		m1 <- mean(hk1[,i])
		hkVEC[iTH] <- m1
	}
	nano_tool_data <- list(nnGeneList3 = nnGeneList3, hk_post_e2 = hk1)
	genes <- unique(unlist(llply(nano_tool_data$ nnGeneList3, rownames)))

	geneMeansTab <- laply(nano_tool_data$ nnGeneList3, function(.x){
			i_all_genes <- match( genes, rownames(.x))
			.x$Expression[i_all_genes]
			})
	colnames(geneMeansTab) <- genes
	meanVec = colMeans(geneMeansTab[1:sampleNumber,], na.rm = TRUE)
	minVec <- apply(geneMeansTab, 2, function(x) min(x, na.rm=TRUE))
		problemGenes <- data.frame(gene=character(0), pvalue=numeric(0), stringsAsFactors=FALSE)
	ind3 <- 0
	for(i in 1:dim(geneMeansTab)[2])
	{
		gene_t_test <- t.test(geneMeansTab[1:sampleNumber,i], geneMeansTab[,i])
		if(gene_t_test$p.value < 0.05)
		{
			ind3 <- ind3 + 1
			problemGenes[ind3,1] <- colnames(geneMeansTab)[i]
			problemGenes[ind3,2] <- gene_t_test$p.value
		}
		else
		{
			# Do Nothing
		}	
	}

	# New function and housekeeper normalisation
	SigmoidHousekeeperNormalisation <- function(nano_tool_data, gene_means = NULL, tissue_type = c('tumour', 'cells'), xmid = -4.2, scal = 1.7){
		#get expresion by inferring slope using a calibration curve and independent gene means
		#gene means is a vector giving average expression levels for each gene after the E2 correction stage (either from a large set of samples or from one 'refernce' samnple.

		genes <- unique(unlist(llply(nano_tool_data$ nnGeneList3, rownames)))

		# this forms an expression matrix form the data held in nnGeneList3 afer the E2 step
		exp_post_e2 <- laply(nano_tool_data$ nnGeneList3, function(.x){
			i_all_genes <- match( genes, rownames(.x))
			.x$Expression[i_all_genes]
			})
		
		colnames(exp_post_e2) <- genes
	
		if(is.null(gene_means)) {
			gene_means = colMeans(exp_post_e2, na.rm = TRUE)
			warning('No gene_means supplied so using means of this data set')
		}
	
		#This forms a vector of avergae houskeeper values
		hk <- colMeans(nano_tool_data$hk_post_e2[,-(1:3)], na.rm = TRUE)
	
		#if(any(names(gene_slopes) != genes)) stop("genes != gene_slopes ") # check this is right!
	
		#This infers the correction slope using the sigmoid calibration curve
		#The default values of xmid and scal have been obtained by lookng at hundreds of tumour samples across hundreds of genes

		#sets sigmoid according to tissue type
		sigmoid_params <- switch(tissue_type, 
			tumour = c(-4.2, 1.7),
			cells = c(-0.79, 1.87))
		
		xmid = sigmoid_params [1]
		scal = sigmoid_params [2]
		
		gene_slopes <- 1 / (1 + exp( (xmid - gene_means) / scal))
	
		#forms a matrix for the normalised expression
		normalised_expression <- matrix(NA, nrow = nrow(exp_post_e2), ncol  = ncol(exp_post_e2))
	
		#Normalieses expression levels using the sigmoid infered slopes
		for(i in 1:ncol(exp_post_e2)){
			gene_i_count <- exp_post_e2[,i]
			slope = gene_slopes[i]
			normalised_expression[,i] <- gene_i_count - slope * hk
		}
	
		#Returns the normlised expression		
		colnames(normalised_expression) <- genes
		t(normalised_expression)
	
	}

	# Run the housekeeper normalisation step!!!
	outMat1 <- SigmoidHousekeeperNormalisation(nano_tool_data, gene_means= meanVec, tissueType)
	colnames(outMat1) <- columnNames
	# Set NA values to minimums and reconstruct the expression matrix
	for(i in 1:dim(outMat1)[1])
	{
		w1 <- which(is.na(outMat1[i,]))
		if(length(w1) > 0)
		{
			for(j in 1:length(w1))
			{
				outMat1[i,w1[j]] <- minVec[i]
			}
		}	
		else
		{
			# Do Nothing
		}
	}

	# Re-scale data (so most of it isn't negative)
	outMat2 <- outMat1
	for(i in 1:dim(outMat2)[2])
	{
		for(j in 1:dim(outMat2)[1])
		{
			n1 <- outMat2[j,i] + 10
			outMat2[j,i] <- n1
		}
	}
	rownames(hk1) <- hk1[,2]
	hk1 <- hk1[,-c(1:3)]
	for(i in 1:dim(hk1)[2])
	{
		for(j in 1:dim(hk1)[1])
		{
			n1 <- hk1[j,i] + 10
			hk1[j,i] <- n1
		}
	}
	outList <- list(Endogenous=outMat2, Housekeeping=hk1, WarningGenes=problemGenes)
	return(outList)
}
 
