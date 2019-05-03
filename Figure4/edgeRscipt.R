library(edgeR)
# define variables
cpmCutoff <- 1
featuresToRemove <- c("alignment_not_unique","ambiguous", "no_feature","not_aligned", "too_low_aQual")
donor <- "donor"
# import sample info see other files in folder
target <- read.table("sampleinfo.txt", header=TRUE, sep="\t", na.strings="")
target[,"group"] <- as.factor(target[,"group"])
target[,"group"] <- relevel(target[,"group"],ref="CCR8mICOSm")
target <- target[order(target[,"group"]),]
rownames(target) <- as.character(target[,1])
# set labels and files (count files)
labels <- as.character(target[,1])
files <- as.character(target[,2])
# creare rawCount df 
rawCounts <- read.table(file.path(files[1]), sep="\t", quote="\"", header=FALSE, skip=0,stringsAsFactors=FALSE)
rawCounts <- rawCounts[,c(1, 2)]
colnames(rawCounts) <- c("Id", labels[1])
# add others samples rawCounts 
for (i in 2:length(files)){
  tmp <- read.table(file.path("/Users/simonepuccio/Documents/HumanitasProjects/SP010/", files[i]), sep="\t", quote="\"", header=FALSE, skip=0, stringsAsFactors=FALSE)
  tmp <- tmp[,c(1, 2)]
  colnames(tmp) <- c("Id", labels[i])
  rawCounts <- merge(rawCounts, tmp, by="Id", all=TRUE)
}
# clean rawCounts
rawCounts[is.na(rawCounts)] <- 0
counts <- as.matrix(rawCounts[,-1])
rownames(counts) <- rawCounts[,1]
counts <- counts[order(rownames(counts)),]
# 
cat("\nFeatures removed:\n")
for (f in setdiff(featuresToRemove,"")){
  match <- grep(f, rownames(counts))
  if (length(match)>0){
    cat(rownames(counts)[match],sep="\n")
    counts <- counts[-match,]
  }
}
cat("\nTop of the counts matrix:\n")
print(head(counts))
cat("\nBottom of the counts matrix:\n")
print(tail(counts))
#### Filtering 
minReplicates <- min(table(target[,"group"]))
fcounts <- counts[rowSums(cpm(counts) >= cpmCutoff) >= minReplicates,]
cat("Number of features discarded by the filtering:\n")
cat(nrow(counts)-nrow(fcounts),"\n")
####### creare design
design <- formula(paste("~", ifelse(!is.null(donor), paste(donor,"+"), ""), "group"))
dge <- DGEList(counts=fcounts, remove.zeros=TRUE)
dge$design <- model.matrix(design, data=target)
cat("\nDesign of the statistical model:\n")
cat(paste(as.character(design),collapse=" "),"\n")	
# 
dge <- calcNormFactors(dge, method="TMM")
cat("\nNormalization factors:\n")
print(dge$samples$norm.factors)
# estimating dispersions
dge <- estimateGLMCommonDisp(dge, dge$design)
dge <- estimateGLMTrendedDisp(dge, dge$design)
dge <- estimateGLMTagwiseDisp(dge, dge$design)
# statistical testing: perform all the comparisons between the levels of varInt
fit <- glmFit(dge, dge$design)
cat(paste("Coefficients of the model:",paste(colnames(fit$design),collapse="  ")),"\n")
colsToTest <- grep("group",colnames(fit$design))
namesToTest <- paste0(gsub("group","",colnames(fit$design)[colsToTest]),"_vs_","CCR8mICOSm")
results <- list()
# testing coefficients individually (tests againts the reference level)
for (i in 1:length(colsToTest)){
  cat(paste0("Comparison ",gsub("_"," ",namesToTest[i]),": testing coefficient ",colnames(fit$design)[colsToTest[i]]),"\n")
  lrt <- glmLRT(fit, coef=colsToTest[i])
  results[[namesToTest[i]]] <- topTags(lrt,n=nrow(dge$counts),adjust.method="BH",sort.by="none")$table
}
# defining contrasts for the other comparisons (if applicable)
if (length(colsToTest)>=2){
  colnames <- gsub(varInt,"",colnames(fit$design))
  for (comp in combn(length(colsToTest),2,simplify=FALSE)){ 
    contrast <- numeric(ncol(dge$design))
    contrast[colsToTest[comp[1:2]]] <- c(-1,1)
    namecomp <- paste0(colnames[colsToTest[comp[2]]],"_vs_",colnames[colsToTest[comp[1]]])
    cat(paste0("Comparison ",gsub("_"," ",namecomp),": testing contrast (",paste(contrast,collapse=", "),")"),"\n")
    lrt <- glmLRT(fit, contrast=contrast)
    results[[namecomp]] <- topTags(lrt,n=nrow(dge$counts),adjust.method="BH",sort.by="none")$table
  }
}

# compute normalization and et baseMean
res <- results
tmm <- dge$samples$norm.factors
N <- colSums(dge$counts)
f <- tmm * N/mean(tmm * N)
normCounts <- round(scale(dge$counts, center=FALSE, scale=f))
base <- data.frame(Id=rownames(counts), counts)
names(base) <- c("Id", colnames(counts))
norm.bm <- data.frame(Id=rownames(normCounts),normCounts)
names(norm.bm) <- c("Id", paste0("norm.",colnames(normCounts)))
norm.bm$baseMean <- round(apply(scale(dge$counts, center=FALSE, scale=f),1,mean),2)
for (cond in levels(target[,"group"])){
  norm.bm[,cond] <- round(apply(as.data.frame(normCounts[,target[,"group"]==cond]),1,mean),0)
}
base <- merge(base,norm.bm,by="Id",all=TRUE)
# results = table with pvalue, FDR and logFC
results
# normCounts = table with normalize count
normCounts
# rawCounts
dge$counts

