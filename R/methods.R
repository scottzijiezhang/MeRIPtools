#################################################################################################################################################3
#' @title callPeakFisher
#' @param MeRIP The MeRIP object from countReads function.
#' @param min_counts The minimal number of reads present in a bin to be called a peak.
#' @param peak_cutoff_fdr The cutoff of fdr of fisher's exact test to call peak.
#' @param peak_cutoff_oddRatio The minimal oddRatio of fisher's exact test to call peak.
#' @param threads The number of threads to use. Default uses 1 threads.
#' @export
callPeakFisher <- function(MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05 , peak_cutoff_oddRatio = 1, threads = 1){
  if(! is(MeRIP, "MeRIP") ){
    stop("The input MeRIP must be a MeRIP dataset!")
  }else if( is(MeRIP, "MeRIP.Peak") ){
    cat(paste0("Input is an object of MeRIP.Peak, will override the peakCallResult by current call of fisher exact test!\n"))
  }else{
    cat("Performing fisher exact test on MeRIP dataset...\n")
  }
  input <- as.matrix(MeRIP@reads[,1:length(MeRIP@samplenames)])
  ip <- as.matrix(MeRIP@reads[,(1+length(MeRIP@samplenames)):(2*length(MeRIP@samplenames))])
  colnames(input) <- colnames(ip) <- MeRIP@samplenames
  ## check if geneBins already exist
  geneBins <- geneBins(MeRIP)

  ## number of genes to call peaks
  batch_id_list <- unique(geneBins$gene)
  num_batch_ids <- length(batch_id_list)
  cat("Calling peaks for ",num_batch_ids, "genes... \n")
  registerDoParallel( cores = threads)
  start_time <- Sys.time()
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to call peaks in continuous bins...\n"))
  peak_call_batches <- foreach( i = 1:num_batch_ids, .combine = rbind)%dopar%{

    idx_batch <- which(geneBins$gene == batch_id_list[i])
    batch_input <- input[idx_batch,]
    batch_ip <- ip[idx_batch,]
    overall_input <- round( apply(batch_input, 2, median, na.rm = TRUE)  )
    overall_ip <- round( apply(batch_ip, 2, median, na.rm = TRUE)  )

    ## loop through all sample
    fisher_exact_test_p <- NULL
    fisher_exact_test_oddRatio <- NULL
    for(j in 1:length(overall_input) ){
      fisher_result <- t( mapply(.fisher_exact_test, batch_ip[,j], batch_input[,j], overall_ip[j], overall_input[j]) )
      fisher_exact_test_p <- cbind(fisher_exact_test_p,fisher_result[,1])
      fisher_exact_test_oddRatio <- cbind(fisher_exact_test_oddRatio,fisher_result[,2] )
    }

    above_thresh_counts <- ( (batch_input + batch_ip) >= min_counts )

    fisher_exact_test_fdr <- matrix(1,nrow = nrow(fisher_exact_test_p),ncol = ncol(fisher_exact_test_p))
    if(sum(rowSums(above_thresh_counts)> (length(overall_input)/2))>1){
      fisher_exact_test_fdr[rowSums(above_thresh_counts)> (length(overall_input)/2) ,] <- apply(fisher_exact_test_p[which(rowSums(above_thresh_counts)> (length(overall_input)/2)) ,] , 2, p.adjust, method = 'fdr')
    }

    fisher_exact_test_peak <- (fisher_exact_test_fdr < peak_cutoff_fdr &
                                 fisher_exact_test_oddRatio > peak_cutoff_oddRatio &
                                 above_thresh_counts)

    fisher_exact_test_peak
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to call peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  colnames(peak_call_batches) <- MeRIP@samplenames

  ## output data
  data.out <- as(MeRIP,"MeRIP.Peak")
  data.out@geneBins <- geneBins
  data.out@peakCallResult <- peak_call_batches
  data.out@peakCalling <- "fisher's exact test"
  return(data.out)

}


# helper function 
.fisher_exact_test <- function(IP, input, IP_overall, input_overall, pseudo_count=1){

  test.m <- matrix(c(IP, input, IP_overall, input_overall), nrow = 2, byrow = FALSE,
                   dimnames = list(c("IP", "input"), c("bin", "overall")))

  # add a pseudo_count (1) to avoid zeros in denominators when calculating the odds ratio
  test.m <- test.m + pseudo_count

  fisher_test <- fisher.test(test.m, alternative="greater")

  fisher_result <- data.frame(pvalue = NA, odds_ratio = NA)

  fisher_result$pvalue <- fisher_test$p.value
  # fisher_result$odds_ratio <- fisher_test$estimate

  # use sample odds ratio (unconditional MLE) rather than the conditional Maximum Likelihood Estimate from Fisher's exact test
  a <- test.m[1,1] # IP in bin
  b <- test.m[1,2] # IP overall
  c <- test.m[2,1] # input in bin
  d <- test.m[2,2] # input overall

  # fisher_result$odds_bin <- a/c

  fisher_result$odds_ratio <- (a*d)/(b*c)

  return(fisher_result)
}

################################################################################################################################################################
#' @export
#' @rdname IP.files
setMethod("IP.files", signature("MeRIP"), function(object){
  object@bamPath.ip
})

#' @export
#' @rdname Input.files
setMethod("Input.files", signature("MeRIP"), function(object){
  object@bamPath.input
})

#' @export
#' @rdname counts
setMethod("counts", signature("MeRIP.Peak"), function(object){
  object@reads
})


#######################################################################################################################################################################################################################

#' @title reportJointPeak
#' @param MeRIPdata The MeRIP.Peak object containing peak calling result
#' @param joint_threshold Define the number of sample required to have consistent peak in a locus to call joint peak.
#' @return MeRIP.Peak object with jointPeaks data
#' @export
reportJointPeak <- function(MeRIPdata, joint_threshold = 2,threads = 1){
  ## check parameter
  if(joint_threshold > length(MeRIPdata@samplenames)){
    stop("joint_threshold cannot be larger than the total number of samples!")
  }
  ## check data format
  if(! is( MeRIPdata,"MeRIP.Peak") ){
    if( is(MeRIPdata, "MeRIP") ){
      cat("No peak calling has performed on input MeRIP object, perform default fisher's exact test to call peak with default parameters...\n")
      MeRIPdata <- callPeakFisher(MeRIPdata, threads = threads )
    }else{
          stop("The input needs to be an MeRIP class object!")
        }
  }else if(  MeRIPdata@jointPeak_threshold > 0  &  MeRIPdata@jointPeak_threshold == joint_threshold ){
    cat(paste0("Reporting joint peak at joint threshold "), MeRIPdata@jointPeak_threshold,"\n")
    geneBins <- geneBins(MeRIPdata)
    peak_id_pairs <- MeRIPdata@jointPeak_id_pairs

    num_peaks <- nrow(peak_id_pairs)
    geneGRList <- MeRIPdata@geneModel
    peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])

    if (num_peaks == 0){return(data.frame())
    }else {
      start_time <- Sys.time()
      registerDoParallel(cores = threads)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- reduce ( geneGRList[peakGenes[p]][[1]] )

        peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),MeRIPdata@binSize )

        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = 0,
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=",")
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }
    
    data.out <- MeRIPdata
    data.out@jointPeaks <- merged.report
    data.out@jointPeak_id_pairs <- peak_id_pairs
    ##Filter out duplicated peaks due to duplicated gene names
    data.out <- filter(data.out, !duplicated( paste(data.out@jointPeaks$chr, data.out@jointPeaks$start, data.out@jointPeaks$end,data.out@jointPeaks$strand, sep = ":") ) )
    return( data.out )
    
  }else if( MeRIPdata@jointPeak_threshold == 0  |  MeRIPdata@jointPeak_threshold != joint_threshold ){
    cat(paste0("Reporting joint peak at joint threshold "), joint_threshold ,"\n")
    
    ## Get joint peak
    geneBins <- geneBins(MeRIPdata)
    ## set logic vector for joint peak
    ID <- (rowSums(MeRIPdata@peakCallResult) >= joint_threshold)
    num_lines <- length(ID)

    # start ids of checkpoints
    ## find peak-starting checkpoints (either from nonpeak to peak, or peak in a different batch)
    start_id <- which((ID[2:num_lines]-ID[1:num_lines-1]==1) |
                        ((geneBins$gene[2:num_lines]!=geneBins$gene[1:num_lines-1]) & (ID[2:num_lines] == TRUE)) )
    start_id <- start_id + 1 # add 1 since ID was counted from 2 to num_lines
    if ( ID[1]==TRUE ) { start_id <- c(1,start_id) } # if the first checkpoint bin is peak

    # end ids of checkpoints
    ## find peak-ending checkpoints (either from peak to nonpeak, or peak in a different batch)
    end_id <- which((ID[1:num_lines-1]-ID[2:num_lines]==1) |
                      ((geneBins$gene[1:num_lines-1]!=geneBins$gene[2:num_lines]) & (ID[1:num_lines-1] == TRUE)) )
    if (ID[num_lines]==TRUE) {end_id <- c(end_id,num_lines)} # if the last checkpoint bin is peak

    peak_id_pairs <- cbind(start_id, end_id)
    num_peaks <- nrow(peak_id_pairs)
    geneGRList <- MeRIPdata@geneModel
    peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])
    #geneBins <- .getGeneBins(geneGRList,peakGenes,MeRIPdata@binSize )

    if (num_peaks == 0){return(data.frame())
    }else{
      start_time <- Sys.time()
      registerDoParallel(cores = threads)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- reduce ( geneGRList[peakGenes[p]][[1]] )

        peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),MeRIPdata@binSize )

        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = 0,
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=",")
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }

    data.out <- MeRIPdata
    data.out@jointPeaks <- merged.report
    data.out@jointPeak_id_pairs <- peak_id_pairs
    data.out@jointPeak_threshold <- joint_threshold
    ## reset the data associated with jointPeaks because new joint peaks has been changed.
    if(MeRIPdata@jointPeak_threshold != joint_threshold & MeRIPdata@jointPeak_threshold >0){
      cat(paste0("joint threshold was previous set at ",MeRIPdata@jointPeak_threshold,".\nWill remove joint-peak read count,test statistics etc. to avoid inconsistency with new joint peaks. Please re-fetch joint peak counts!\n"))
      data.out@jointPeak_ip <-  new("matrix")
      data.out@jointPeak_input <-  new("matrix")
      data.out@norm.jointPeak_ip <-  new("matrix")
      data.out@jointPeak_adjExpr <-  new("matrix")
      data.out@test.est <-  new("matrix")
      data.out@test.method <-  "none"
    }
    ##Filter out duplicated peaks due to duplicated gene names
    data.out <- filter(data.out, !duplicated( paste(data.out@jointPeaks$chr, data.out@jointPeaks$start, data.out@jointPeaks$end,data.out@jointPeaks$strand, sep = ":") ) )
    return( data.out )
    
  }


}


.peakExons <- function(peak,y){
  exonID <- peak$start <= y$end & peak$end >= y$start
  if(sum(exonID) == 1){
    return(data.frame(start = peak$start, end = peak$end, width = peak$end - peak$start))
  }else if(sum(exonID) > 1){
    peakexon <- y[exonID,]
    peakexon[1,"start"] <- peak$start
    peakexon[sum(exonID),"end"] <- peak$end
    return(data.frame(start = peakexon$start, end = peakexon$end, width = peakexon$end - peakexon$start + 1))
  }
}


.getGeneBins <- function(geneGRList,geneNames,binSize,threads = 1){
  
  registerDoParallel(cores = threads)
  geneBins <-foreach(i = 1:length(geneNames), .combine = rbind)%dopar%{
    geneName = geneNames[i]
    geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons
    
    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel) )
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    #creat a corresponding map from RNA to DNA
    RNA2DNA = 1:exon.length
    pointer = 1
    for (j in 1:no.exon){
      RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
      pointer = pointer + df.geneModel$width[j]
    }
    RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates
    #creat center points of continuous window
    if(exon.length <= binSize){
      slidingStart= exon.length/2
      mapping = data.frame(start = RNA2DNA[slidingStart-exon.length/2+1], end = RNA2DNA[slidingStart + exon.length/2]  )
    }else{
      slidingStart= round(seq(from = binSize/2, to = (exon.length - binSize/2), length.out = ceiling(exon.length/binSize) ) )
      mapping = data.frame(start = RNA2DNA[slidingStart - binSize/2 +1], end = RNA2DNA[slidingStart + binSize/2 ]  )
    }
    mapping$chr = as.character(dna.range$seqnames)
    mapping$strand = as.character(dna.range$strand)
    cbind(data.frame(geneName,slidingStart),mapping)
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  return(geneBins)
}

.getPeakBins <- function(geneGRList,geneName,slidingStarts,binSize){
  
  geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons
  
  # DNA location to gene location conversion
  df.geneModel= as.data.frame(geneModel) ##data frame of gene model
  dna.range = as.data.frame(range(geneModel) )
  df.geneModel$end = df.geneModel$end - dna.range$start + 1
  df.geneModel$start = df.geneModel$start - dna.range$start + 1
  DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
  no.exon = dim(df.geneModel)[1]
  for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
  exon.length = sum(DNA2RNA)
  #creat a corresponding map from RNA to DNA
  RNA2DNA = 1:exon.length
  pointer = 1
  for (j in 1:no.exon){
    RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
    pointer = pointer + df.geneModel$width[j]
  }
  RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates
  #creat center points of continuous window
  if(exon.length <= binSize){
    slidingStart= exon.length/2
    mapping = data.frame(start = RNA2DNA[slidingStarts[1]-exon.length/2+1], end = RNA2DNA[slidingStarts[2] + exon.length/2]  )
  }else{
    mapping = data.frame(start = RNA2DNA[slidingStarts[1] - binSize/2 +1], end = RNA2DNA[slidingStarts[2] + binSize/2 ]  )
  }
  
  mapping$chr = as.character(dna.range$seqnames)
  return(mapping[,c("chr","start","end")])
  
}

########################################################################################################################################################


#' @title jointPeakCount
#' @description Extract read count of joint peaks for each sample and store the read counts in the MeRIP.Peak object. 
#' @param MeRIPdata The data list as output of callPeakFisher()
#' @return Returns the MeRIPdata with joint peak read count for input and IP stored.
#' @export
jointPeakCount <- function(MeRIPdata ){
  if( is( MeRIPdata,"MeRIP.Peak") ){
    if( nrow(MeRIPdata@jointPeak_id_pairs) >0 & MeRIPdata@jointPeak_threshold> 0){}else{stop("Please report JointPeak first!")}
  }else{
    stop("The MeRIPdata needs to be a MeRIP.Peak object!")
  }
  
  geneBins <- geneBins(MeRIPdata)

  peak_id_pairs <- MeRIPdata@jointPeak_id_pairs
  

  ip <- MeRIPdata@reads[,(1+length(MeRIPdata@samplenames)):(2*length(MeRIPdata@samplenames))]
  input <- MeRIPdata@reads[,1:length(MeRIPdata@samplenames)]
  
  joint_peak_ip <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = ip) )
  
  joint_peak_input <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = input) )
  
  data.out <- MeRIPdata
  data.out@jointPeak_ip <- joint_peak_ip
  data.out@jointPeak_input <- joint_peak_input
  
  return(data.out)
  
}

#' @export
#' @title extract geneBins
#' @description geneBins extractor
setMethod("geneBins", signature("MeRIP"), function(object){
  if( nrow(object@geneBins) > 0 ){
    return(object@geneBins)
  }else{
    ## split gene and bin names
    aa <- strsplit(rownames(object@reads), ",")
    gene.name <- unlist(lapply(aa, function(x){
      return(x[1])
    }))
    bin.name <- unlist(lapply(aa, function(x){
      return(x[2])
    }))
    geneBins <- data.frame(gene=gene.name,bin=as.integer(bin.name))
    rownames(geneBins) <- rownames(object@reads)
    return(geneBins)
  }
})
#' @export
setMethod("geneBins", signature("MeRIP.Peak"), function(object){
  callNextMethod()
})

##########################################################################################################################################################

#' @title normalizeLibrary
#' @description Normalized the input as RNA-seq data and normalize IP by enrichment. Specifically, we normalize ip libraries sizes so that the geometry mean of enrichment are the same.
#' @param object MeRIP.Peak object.
#' @import DESeq2
#' @export
setMethod("normalizeLibrary", signature("MeRIP.Peak"), function(object){
  
  ## load data from input
  jointPeak_ip <- object@jointPeak_ip
  jointPeak_input <- object@jointPeak_input
  input <- object@reads[,1:length(object@samplenames)]
  colnames(input) <- object@samplenames
  geneBins <- geneBins(object)
  
  ## Get input geneSum (gene level quantification)
  geneSum <- NULL
  for(i in 1:ncol(input) ){
    y <- input[,i]
    gene.sum <- tapply(y,geneBins$gene,sum)
    geneSum <- cbind(geneSum,gene.sum)
  }
  colnames(geneSum) <- object@samplenames
  
  size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
  
  norm.jointPeak_input <-t( t(jointPeak_input) / size.input )
  geneSum.norm <- t ( t(geneSum)/size.input)
  
  ## Get the gene level input count for corresponding peaks
  geneCounts.peak <- geneSum.norm[geneBins[rownames(norm.jointPeak_input),"gene"],]
  enrich <- as.data.frame(jointPeak_ip/geneCounts.peak)
  enrich <- enrich[!apply(enrich,1, function(x){any(is.na(x)) | any(is.infinite(x))}),]
  
  size.enrich.deseq2 <- DESeq2::estimateSizeFactorsForMatrix(enrich[,1:length(object@samplenames)])
  
  norm.jointPeak_ip <-t( t(jointPeak_ip)/size.enrich.deseq2 )
  sizeFactor <- data.frame(input=size.input,ip=size.enrich.deseq2)
  
  object@geneSum <- geneSum.norm
  object@norm.jointPeak_ip <- norm.jointPeak_ip
  object@sizeFactor <- sizeFactor
  return(object)
  
})


#######################################################################################################################################################

## a helper function to remove 0 from vector
.noZero <- function(x){sapply(x,max,1)}

#' @title adjustExprLevel
#' @param object The MeRIP.Peak object that has been normalized for library size.
#' @param adjustBy By default, adjust post-IP count by INPUT geneSum. Can also choose "pos" to use current position count to adjust for expression level.
#' @return  object The MeRIP.Peak object now with IP-count adjusted for expression level.
#' @export
setMethod("adjustExprLevel", signature("MeRIP.Peak"), function(object, adjustBy = "geneSum" ){
  if( nrow(object@norm.jointPeak_ip)<0 ){
    stop("Please normalize library size before running expression level adjustment!")
  }else if(adjustBy == "geneSum"){
    geneSum <- geneExpression(object)
    geneSum <- t(apply(geneSum,1,.noZero))
    gene.size <- t( apply(geneSum,1,function(x){x/mean(x)}) )
    gene.size.factor <- gene.size[jointPeak(object)$name,]
    ip_norm_geneSum <- object@norm.jointPeak_ip/gene.size.factor
    ip_norm_geneSum <- round(ip_norm_geneSum)
    object@jointPeak_adjExpr <- ip_norm_geneSum
    return(object)
  }else if(adjustBy == "pos"){
    norm.jointPeak_input <-t( t(object@jointPeak_input) / object@sizeFactor$input )
    norm.jointPeak_input <- t(apply(norm.jointPeak_input,1,.noZero))
    pos.size <-  t( apply(norm.jointPeak_input,1,function(x){x/mean(x)}) )
    ip_norm_pos <- object@norm.jointPeak_ip/pos.size
    ip_norm_pos <- round(ip_norm_pos)
    object@jointPeak_adjExpr <- ip_norm_pos
    return(object)
  }else{
    stop("Must specify adjustBy = \"geneSum\" or by \"pos\"...")
  }
})

############################################################################################################################################
#' @title consistentPeak
#' @param object The MeRIP.Peak object contain peak calling result.
#' @param samplenames The samplenames to be reported for consistent peaks.
#' @param joint_threshold Define the number of sample required to have consistent peak in a locus to call consistent peak in a group.
#' @return Peak consistent across specified samples at specified joint_threshod. If joint_threshold not specified, report consistent peaks across all samples specified. If no sample specified, report consistent peak across all samples.
#' @export
setMethod("consistentPeak", signature("MeRIP.Peak"), function(object, samplenames = NULL, joint_threshold = NA, threads = 1){
  if(!is.null(samplenames) & ! all(samplenames %in% object@samplenames ) ){
    stop("Not all samples specified are in the input dataset!")
  }else if( !is.null(samplenames) ){
    ## define samples to be reported
    sample_report_id <- match(samplenames, object@samplenames )
    
    if(is.na(joint_threshold)){
      joint_threshold <- length(samplenames)
      cat("Reporting peak concsistent in all samples for\n",paste(samplenames,collapse = " "),"\n")
    }else if(joint_threshold > length(samplenames) ){
      joint_threshold <- length(samplenames)
      cat("Reporting peak concsistent in all samples for\n",paste(samplenames,collapse = " "),"\n")
    }else{
      cat("Reporting joint peak consistant in at least ",joint_threshold," samples among \n",paste(samplenames,collapse = " "),"\n")
    }
  }else{ 
    sample_report_id <- 1:length( object@samplenames )
    joint_threshold <- ifelse(is.na(joint_threshold),length( object@samplenames ), joint_threshold )
    }
  
  if( nrow(object@peakCallResult)>0 & object@peakCalling != "none" ){
    
    ## Get joint peak
    geneBins <- geneBins(object)
    
    ## set logic vector for joint peak
    if(length(sample_report_id) >1 ){
      ID <- (rowSums(object@peakCallResult[,sample_report_id]) >= joint_threshold)
    }else if(length(sample_report_id) == 1){
      ID <- object@peakCallResult[,sample_report_id]
    }
    
    num_lines <- length(ID)
    
    # start ids of checkpoints
    ## find peak-starting checkpoints (either from nonpeak to peak, or peak in a different batch)
    start_id <- which((ID[2:num_lines]-ID[1:num_lines-1]==1) |
                        ((geneBins$gene[2:num_lines]!=geneBins$gene[1:num_lines-1]) & (ID[2:num_lines] == TRUE)) )
    start_id <- start_id + 1 # add 1 since ID was counted from 2 to num_lines
    if ( ID[1]==TRUE ) { start_id <- c(1,start_id) } # if the first checkpoint bin is peak
    
    # end ids of checkpoints
    ## find peak-ending checkpoints (either from peak to nonpeak, or peak in a different batch)
    end_id <- which((ID[1:num_lines-1]-ID[2:num_lines]==1) |
                      ((geneBins$gene[1:num_lines-1]!=geneBins$gene[2:num_lines]) & (ID[1:num_lines-1] == TRUE)) )
    if (ID[num_lines]==TRUE) {end_id <- c(end_id,num_lines)} # if the last checkpoint bin is peak
    
    ## get peak id-pairs that defines the bins to be merged as one peak
    peak_id_pairs <- cbind(start_id, end_id)
    num_peaks <- nrow(peak_id_pairs)
    geneGRList <- object@geneModel
    peakGenes <- as.character(geneBins[peak_id_pairs[,1],"gene"])
    
    ## Get raw readcount for peaks.
    ip <- object@reads[,(length(object@samplenames)+sample_report_id)]
    input <- object@reads[,sample_report_id]
    ## standardize library size
    all.size <- colSums( cbind(input,ip) )/mean( colSums( cbind(input,ip) ) )
    
    if(length(sample_report_id) >1 ){
      ## get Summed read count
      ip_sum <- round( rowSums( t( t(ip)/all.size[ length(sample_report_id)+1:length(sample_report_id) ] ) ) )
      input_sum <- round(rowSums( t( t(input)/all.size[ 1:length(sample_report_id) ] ) ) )
    }else if(length(sample_report_id) == 1){
      ## No need to sum for single sample
      ip_sum <- round( ip/all.size[2] )
      input_sum <- round( input/all.size[1] )
    }
    
    
    ## Get peak ip count
    joint_peak_ip <- apply(peak_id_pairs,1,function(x,y){
      if(x[1]==x[2]){
        return(y[x[1]:x[2]])
      }else{
        sum(y[x[1]:x[2]])
      }
    },y = ip_sum)
    
    ## Get peak input count
    joint_peak_input <- apply(peak_id_pairs,1,function(x,y){
      if(x[1]==x[2]){
        return(y[x[1]:x[2]])
      }else{
        sum(y[x[1]:x[2]])
      }
    },y = input_sum)
    
    ## Get peak gene ip median
    joint_peak_ip_median <-  round( tapply(ip_sum,geneBins$gene,median)[peakGenes]  )
    ## Get peak gene input median
    joint_peak_input_median <- round( tapply(input_sum,geneBins$gene,median)[peakGenes]  )
    
    ## Compute fisher's test for each peak
    peak_test <- do.call( rbind.data.frame, apply(cbind(joint_peak_ip,joint_peak_input,joint_peak_ip_median,joint_peak_input_median), 1, function(x){.fisher_exact_test(x[1],x[2],x[3],x[4])}) )
    
    if (num_peaks == 0){return(data.frame())
    }else {
      start_time <- Sys.time()
      registerDoParallel(cores = threads)
      cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
      cat(paste("Using",getDoParWorkers(),"thread(s) to report merged report...\n"))
      merged.report<- foreach( p = 1:num_peaks, .combine = rbind)%dopar%{
        peak_row_id <- peak_id_pairs[p,]
        geneExons <- reduce ( geneGRList[peakGenes[p]][[1]] )
        
        peak <- .getPeakBins(geneGRList,peakGenes[p],c(geneBins$bin[peak_row_id[1]],geneBins$bin[peak_row_id[2]]),object@binSize )
        
        peakE <- .peakExons(peak,as.data.frame(geneExons))
        data.frame(chr=peak$chr,
                   start = peak$start,
                   end = peak$end,
                   name = peakGenes[p],
                   score = peak_test$pvalue[p],
                   strand = as.character(strand(geneExons))[1],
                   thickStart = peak$start,
                   thickEnd = peak$end,
                   itemRgb=0,
                   blockCount = nrow(peakE),
                   blockSizes = paste(peakE$width,collapse=","),
                   blockStarts = paste(peakE$start - replicate(nrow(peakE),peakE$start[1]),collapse=","),
                   enrichmentScore = peak_test$odds_ratio[p]
        )
      }
      rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
      end_time <- Sys.time()
      cat(paste("Time used to report peaks:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
    }
    
    return( merged.report )
    
  } else{
    stop("Please run callPeakFisher() before reporting peaks...")
  }
  
})


########################################################################################################################################################
#' @title Wraper function to use QNB package to test for differential peaks.
#' @description 
#' @name QNBtest
#' @param object The MeRIP.Peak object
#' @import QNB
#' @export
setMethod("QNBtest", signature("MeRIP.Peak"), function(object){
  if( nrow(object@variate) != length(object@samplenames)  ){
    stop(" Predictor variable lengthen needs to match the sample size! If you haven't set the predictor variable, please set it by variable(object) <- data.frmae(group = c(...)) ")
  }else if(length(unique(variable(object)[,1])) > 2 ){
    stop("The levels of predictor variable needs to be two!")
  }
  
  Ctl_id <- which( object@variate[,1] %in%  unique(object@variate[,1])[1] )
  Treat_id <-which( object@variate[,1] %in%  unique(object@variate[,1])[2] )
  QNB.res <- QNB::qnbtest(control_ip = object@jointPeak_ip[,Ctl_id],
                          treated_ip = object@jointPeak_ip[,Treat_id],
                          control_input = object@jointPeak_input[,Ctl_id],
                          treated_input = object@jointPeak_input[,Treat_id],plot.dispersion = FALSE)
  colnames(QNB.res) <- c("p.treated","p.control", "log2.RR" ,"beta",  "p_value", "q","padj" )
  
  object@test.est <- QNB.res
  object@test.method <- "QNB test"
  cat("The test was performed by R package QNB, please cite the QNB paper:\nLian Liu (2017). QNB: Differential RNA Methylation Analysis for Count-Based Small-Sample Sequencing Data with a Quad-Negative Binomial
  Model\nif you used this result in you publication!")
  return(object)
})


######################################################################################################################################

#' @title plot distribution of peaks on gene annotation
#' @name peakDistribution
#' @description plot distribution of peaks on gene annotation
#' @param object The MeRIP.Peak object
#' @param saveName the file name to save ditribution plot
#' @export
setMethod("peakDistribution", signature("MeRIP.Peak"), function(object){
  ## collapes the peak to the center
  x <- jointPeak(object)
  for(i in 1:nrow(x)){
    start = x[i,2]
    end = x[i,3]
    blocks = x[i,10]
    blockStart = as.numeric( unlist(strsplit(as.character(x[i,12]),",")) )
    blockSize = as.numeric( unlist(strsplit(as.character(x[i,11]),",")) )
    if(blocks <2){
      x[i,2] =  x[i,3] = round(start + sum(blockSize)/2 )
    } else{
      pointer = 2
      while(sum(blockSize[1:pointer]) < sum(blockSize)/2 ){pointer = pointer+1}
      x[i,2] =  x[i,3] = round(start + blockStart[pointer] + sum(blockSize)/2 - sum(blockSize[1:(pointer-1)]) )
    }
    x[i,10] = 1
  }
  
  ################
  gr.peak <- reduce( makeGRangesFromDataFrame(x) )
  txdb <- makeTxDbFromGFF(file = object@gtf,format = "gtf")
  
  cds <-  cdsBy(txdb,by = "tx")
  n.cds <- length( which(countOverlaps(gr.peak,cds)>0) )
  
  fiveUTR <- fiveUTRsByTranscript(txdb)
  n.fiveUTR <- length( which(countOverlaps(gr.peak,fiveUTR)>0) )
  
  threeUTR <- threeUTRsByTranscript(txdb)
  n.threeUTR <-  length( which(countOverlaps(gr.peak,threeUTR)>0) )
  
  exon <- exonsBy(txdb,by = "tx")
  all_mRNA <- unique(c(names(fiveUTR),names(threeUTR),names(cds)))
  name_ncRNA <- setdiff(names(exon),all_mRNA)
  ncRNA <- exon[name_ncRNA]
  n.ncRNA <-  length( which(countOverlaps(gr.peak,ncRNA)>0) )
  
  slices <- c(n.fiveUTR,n.cds,n.threeUTR,n.ncRNA)
  pct <- round(slices/sum(slices)*100)
  lbls <- c("5'UTR","CDS","3'UTR", "ncRNA")
  lbls <- paste(lbls, pct,sep = "\n") # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  distr <- data.frame(Annotation = factor(c("5'UTR","CDS","3'UTR", "ncRNA"),levels = c("5'UTR","CDS","3'UTR", "ncRNA")),
             fraction = pct)
  ggplot(distr,aes( x = factor(1), y = fraction, fill = Annotation)) + geom_bar(width = 1,stat = "identity") +coord_polar( theta = "y")+theme_minimal()+
    theme(axis.title = element_blank(), axis.text = element_blank(),axis.title.y = element_blank(),panel.grid=element_blank(), legend.title = element_text(face = "bold"),legend.text = element_text(face = "bold")) + 
    geom_text(aes(y = rev( fraction/2 + c(0, cumsum(fraction)[-length(fraction)])),label = lbls), size=4)
})

######################################################################################################################################
#' @export
#' @title Prepare coverage plot
#' @description import GTF into the MeRIP object for plot
#' @param object The MeRIP object
#' @param gtf optional gtf file if the stored path to gtf file has changed.
setMethod("PrepCoveragePlot", signature("MeRIP.Peak"), function(object , gtf = NULL){
  ## check gtf slot
  if( file.exists(object@gtf) ){
    object@GTF <- rtracklayer::import(object@gtf, format = "gtf")
    return(object)
  }else if(! is.null(gtf) ){
    object@GTF <- rtracklayer::import( gtf, format = "gtf")
    object@gtf <- gtf
    cat(paste0("assigning new path to gtf file: ",gtf," \n"))
    return(object)
  }else{
    stop("The gtf file doesn't exist! Please suply path to gtf file in by PrepCoveragePlot(MeRIP, gtf = \"path/to/gtf\" )")
  }
})

###################################################################################################################################
#' @export
setMethod("extractIP", signature("MeRIP.Peak"), function(object, normalized = FALSE, adjusted = FALSE ){
  if(nrow(object@jointPeak_ip)>0){
    if( adjusted ){
      if(nrow(object@jointPeak_adjExpr)>0){
        return( object@jointPeak_adjExpr )
      }else{stop("IP read count has not yet been adjusted for expression variation. Call adjustExprLevel() first!")}
    }else if(normalized ){
      if(nrow(object@norm.jointPeak_ip)>0){
        return( object@norm.jointPeak_ip )
      }else{stop("IP read count has not yet been normalzied. Call normalizedPeak() first!")}
    }else{
      return(object@jointPeak_ip)
    }
  }else{
    stop("JointPeaks counts not reported, please call jointPeakCount() to extract JointPeak count first.\n")
  }
})

#' @export
setMethod("extractInput", signature("MeRIP.Peak"), function(object){
  if(nrow(object@jointPeak_input)>0){
    return(object@jointPeak_input)
  }else{
    stop("JointPeaks counts not reported, please call jointPeakCount() to extract JointPeak count first.\n")
  }
})

#######################################################################################################################################
#' @title plotGeneCov
#' @param object The data list from countReads and other analysis.
#' @param geneName The gene symbol to be ploted.
#' @param GTF The GRanges object containing gtf annotation. Can obtain by rtracklayer::import("file.gtf", format= "gtf").
#' @param libraryType Specify whether the library is the same or opposite strand of the original RNA molecule. Default is "opposite".
#' @param center Specify the method to calculate average coverage of each group. Could be mean or median.
#' @param ZoomIn c(start,end) The coordinate to zoom in at the gene to be ploted.
#' @param adjustExprLevel logical parameter. Specify whether normalize the two group so that they have similar expression level.
#' @export
setMethod("plotGeneCov", signature("MeRIP.Peak"), function(object, geneName, libraryType = "opposite", center = mean,ZoomIn = NULL, adjustExprLevel = F ){
  if( nrow(variable(object))>0 ){
  X <- factor(variable(object)[,1])
  plotGeneCoverage(IP_BAMs = IP.files(object),
                   INPUT_BAMs =Input.files(object),
                   size.IP = object@sizeFactor$ip,
                   size.INPUT = object@sizeFactor$input,
                   X, geneName,
                   geneModel = object@geneModel,
                   libraryType, center  ,object@GTF, ZoomIn, adjustExprLevel, plotSNP = NULL  )+
    theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),legend.title =  element_text(hjust = 0.5,size = 13,face = "bold"),legend.text =  element_text(size = 12,face = "bold"))
}else{
  plotGeneCoverage(IP_BAMs = IP.files(object),
                   INPUT_BAMs = Input.files(object),
                   size.IP = object@sizeFactor$ip,
                   size.INPUT = object@sizeFactor$input,
                   rep("All samples",length(object@samplenames)), geneName,
                   geneModel = object$geneModel,
                   libraryType, center, object@GTF ,ZoomIn, adjustExprLevel, plotSNP = NULL  )+
    theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),legend.position="none" )
}
})


#####################################################################################################################################3
.genotypeDosage <- function(x){
  split <- sapply(x,strsplit,"[|]")
  return( suppressWarnings( unlist(lapply(split,function(x){sum(as.integer(x),na.rm = T)})) ) )
}

.getGenotype <- function(vcf.gz, SNPID){
  tmp <- system2(command = "zcat", args = paste0(vcf.gz," |grep '",SNPID,"' |cat"),stdout = T)
  geno <- strsplit(tmp,split = "\t")[[1]]
  return(geno)
}

#' @name plotSNPpeakPairs
#' @param genotypeFile The path to gzipped vcf file where genotype is available.
#' @param SNPID rsID of the SNP to be ploted
#' @param geneName The name (as defined in gtf file) of the gene you want to plot
#' @param libraryType "opposite" for mRNA stranded library, "same" for samll RNA library
#' @param adjustExprLevel Logic parameter determining whether adjust coverage so that input are at "same" expression level.
#' @export
setMethod("plotSNPpeakPairs", signature("MeRIP.Peak"), function(object, genotypeFile,SNPID, geneName ,libraryType = "opposite", center = mean ,ZoomIn=NULL, adjustExprLevel = TRUE){
  
  if( length(object@GTF) == 0 ){ stop("Please run PrepCoveragePlot(object) first.") }
  
  geno <- .getGenotype(genotypeFile,SNPID)
  
  dosage_geno <- .genotypeDosage(geno[-c(1:9)])
  
  snp_chr <- geno[1]
  snp_loc <- geno[2]
  snp_ref <- geno[4]
  snp_alt <- geno[5]
  
  X <- gsub("0",paste0(snp_ref,snp_ref),as.character(dosage_geno))
  X <- gsub("1",paste0(snp_ref,snp_alt),as.character(X))
  X <- gsub("2",paste0(snp_alt,snp_alt),as.character(X))
  
  X <- factor(X,levels=c(paste0(snp_ref,snp_ref),paste0(snp_ref,snp_alt),paste0(snp_alt,snp_alt))[c(paste0(snp_ref,snp_ref),paste0(snp_ref,snp_alt),paste0(snp_alt,snp_alt))%in%X])
  
  suppressMessages(plotGeneCoverage(IP_BAMs = IP.files(object), 
                                    INPUT_BAMs = Input.files(object), 
                                    size.IP =  object@sizeFactor$ip ,
                                    size.INPUT = object@sizeFactor$input ,
                                    X = X, geneName = geneName, 
                                    geneModel = object@geneModel, 
                                    libraryType = libraryType, 
                                    center = center ,
                                    GTF = object@GTF,
                                    ZoomIn = ZoomIn, adjustExprLevel = adjustExprLevel, plotSNP = data.frame(loc = as.numeric(snp_loc), anno = paste0(snp_ref,"/",snp_alt) )  )+
                     ggtitle(paste0("SNPID: ",SNPID,"    m6A peak on: ",geneName))+theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),legend.title =  element_text(hjust = 0.5,size = 13,face = "bold"),legend.text =  element_text(size = 12,face = "bold"))+
                     scale_y_continuous(expand = c(0.02,0,-0.01,0) )
                   )
})


#######################################################################################################################################################

#' @title getTPM from geneSum
#' @name geneExpressionTMP
#' @param object The MeRIP.Peak object
#' @param meanFragmentLength The mean length of RNA fragment (insert of RNA library). Default is 150bp.
#' @param normalize Logical indicating whether normalized TPM or raw TPM should be returned.
#' @export
setMethod("geneExpressionTMP", signature("MeRIP.Peak") , function(object, meanFragmentLength = 150, normalize = T){
  input <- object@reads[,1:length(object@samplenames)]
  gene.name <- geneBins(object)$gene
  geneSum <- geneExpression(object)
  colnames(geneSum) <- object@samplenames
  
  genes <- rownames(geneSum)
  
  cat("calculating gene length...\n")
  geneLength <- sapply(genes,function(yy){
    sum( as.data.frame( object@geneModel[[yy]] )$width )
  })
  
  
  effLength <- geneLength - meanFragmentLength
  effLength <- sapply(effLength,max,1) ## remove effective length smaller than 0.
  
  cat(paste0("computing TPM from read counts using mean fragment length = ",meanFragmentLength,".\n"))
  
  rate <- apply(geneSum,2,function(aa){aa/effLength})
  totalCounts <-colSums(rate)
  
  tpm <- t( t(rate)/totalCounts ) *1e6
  
  size.tpm <- DESeq2::estimateSizeFactorsForMatrix(tpm)
  tpm_norm <- t(t(tpm)/ size.tpm)
  
  if(normalize){
    return(tpm_norm)
  }else{
    return(tpm)
  }
})

###########################################################################################################################################################

#' @title plotTPM
#' @param TPM Dataframe of gene TPM
#' @param geneName The name of genes to be ploted.
#' @param group Categorical info for each sample.
#' @param logCount where to plot count at log scale
#' @export
plotTPM <- function(TPM,geneName,group,logCount = FALSE, facet_grid = FALSE){
  if(length(geneName) ==1){
    temp <- as.data.frame(t(TPM[geneName,] ) )
  }else{
    temp <- as.data.frame(TPM[geneName,] )
  }
  
  if(logCount){
    temp <- log(temp)
    colnames(temp) <- paste0(group,1:length(group))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(group)[1]
    for(i in 2:length(group)){
      temp_melt$Group[grep(unique(group)[i],temp_melt$variable)] <- unique(group)[i]
    }
    
    axis.font <- element_text(face = "bold", color = "black")
    if(facet_grid){
      ggplot(temp_melt, aes(x= Group,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log TPM")+facet_grid(.~ name)+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_blank(),
                           axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                           legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                           axis.text.x =element_blank() ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                           plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.5,family = "arial"),
                           axis.ticks.x = element_blank(),
                           strip.text.x = element_text(size = 15,face = "bold") )+
        ggtitle("Gene expression level")
    }else{
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log TPM")+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_text(size=20, face="bold", hjust=0.5,family = "arial"),
                           axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                           legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                           axis.text.x = element_text(size = 15,face = "bold",family = "arial",colour = "black") ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                           plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.5,family = "arial"))+
        ggtitle("Gene expression level")
    }
    
  }else{
    colnames(temp) <- paste0(group,1:length(group))
    temp$name <- factor(geneName,levels = geneName)
    temp_melt <- reshape2::melt(temp,id.vars = "name")
    temp_melt$Group <- unique(group)[1]
    for(i in 2:length(group)){
      temp_melt$Group[grep(unique(group)[i],temp_melt$variable)] <- unique(group)[i]
    }
    axis.font <- element_text(face = "bold", color = "black")
    if(facet_grid){
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="TPM")+facet_grid(.~ name)+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_blank(),
                           axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                           legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                           axis.text.x = element_blank() ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                           plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.4,family = "arial"),
                           axis.ticks.x = element_blank(),
                           strip.text.x = element_text(size = 15,face = "bold") )+
        ggtitle("Gene expression level")
    }else{
      ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="TPM")+
        theme(axis.title =axis.font, axis.text = axis.font)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           axis.line = element_line(colour = "black",size = 1),
                           axis.title.x=element_text(size=20, face="bold", hjust=0.5,family = "arial"),
                           axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                           legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                           axis.text.x = element_text(size = 15,face = "bold",family = "arial",colour = "black") ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                           plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.4,family = "arial"))+
        ggtitle("Gene expression level")
    }
    
  }
}

##########################################################################################################################################

#' @title Perform inferential test using Poisson Random effect model in RADAR 
#' @description 
#' @name RADARtest
#' @param object The MeRIP.Peak object
#' @export
setMethod("RADARtest", signature("MeRIP.Peak"), function(object, exclude = NULL, maxPsi = 100){
  
  if( nrow(object@variate) != length(object@samplenames)  ){
    stop(" Predictor variable lengthen needs to match the sample size! If you haven't set the predictor variable, please set it by variable(object) <- data.frmae(group = c(...)) ")
  }else if(length(unique(variable(object)[,1])) > 2 & !is.numeric(variable(object)[,1]) ){
    stop("The levels of predictor variable needs to be two!")
  }
  
  if(!is.null(exclude) & all( exclude %in% samplenames(object)) ){
    object <- select(object, setdiff(samplenames(object), exclude ))
  }
  
  allY <- object@jointPeak_adjExpr
  psi <- 10 # start point
  
  ## convert predictor variable if it is not numeric
  if( is.numeric(variable(object)[,1]) ){
    X <- variable(object)[,1]
  }else{
    tmp <- rafalib::as.fumeric(as.character(variable(object)[,1])) - 1
    names(tmp) <- as.character(variable(object)[,1])
    cat("The predictor variable has been converted:\n")
    print(tmp)
    X <- rafalib::as.fumeric(as.character(variable(object)[,1])) - 1 # convert categorical variable into numerical variable.
  }
  
  if( ncol(variable(object)) == 1 ){
    
    cat("running PoissonGamma test at single beta mode\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ])
      model1 <- glm(Y ~ X, family = poisson(link = 'log'))
      coef <- model1$coefficients
      mu2 <- coef[1]
      beta <- coef[2]
      est <- try(unlist(PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    rownames(all.est) <- rownames(allY)[all.id]
    
    
  }else if(ncol(variable(object)) > 1){
    if(! any(sapply(2:ncol(variable(object)),function(x) is.numeric(variable(object)[,x])) ) ){stop("Please convert all covariates into numerical variables. Discrete variable should be converted to binary variable.")}
    X.all <- as.matrix( cbind(X,variable(object)[,2:ncol(variable(object))])  )
    colnames(X.all) <- colnames(variable(object))
    design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste(colnames(X.all), sep = "",collapse = " + ")) )
    cat("running PoissonGamma test at multi-beta mode...\n")
    pb <- txtProgressBar(min = 1, max = nrow(allY), style = 3) ##creat a progress bar to track loop progress
    all.est <- NULL
    all.id <- NULL
    for(kk in 1:nrow(allY)){
      Y <- unlist(allY[kk, ] )
      aa <- unlist(summary( lm( design.multiBeta ) )$coefficients[, 1])
      mu2 <- aa[1]
      beta <- aa[2:(ncol(X.all)+1 )]
      est <- try(unlist(PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 10, down = 0.1,psi_cutoff = maxPsi)))
      if(class(est) != "try-error"){
        all.est <- rbind(all.est, est)
        all.id <- c(all.id, kk)
      }
      setTxtProgressBar(pb, kk) # update progress bar
    }
    
    rownames(all.est) <- rownames(allY)[all.id] ## assign window names to test statistics
    colnames(all.est) <- gsub("3","",colnames(all.est))
  }
  
  fdr <- p.adjust(all.est[,"p_value"],method = "fdr" )
  
  object@test.est <- cbind(all.est,fdr)
  object@test.method <- "PoissonGamma test"
  cat("\n")
  return(object)
})


##########################################################################################################################################

#' @export
setMethod("samplenames", signature("MeRIP.Peak"), function(object ){
  object@samplenames
})
#' @export
setMethod("samplenames", signature("MeRIP"), function(object ){
  object@samplenames
}) 
#' @export
setMethod("samplenames<-", signature("MeRIP.Peak"), function(object ,value){
  object@samplenames <- value
  object
})
#' @export
setMethod("samplenames<-", signature("MeRIP"), function(object ,value){
  object@samplenames <- value
  object
})

###########################################################################################################################
#' @export
setMethod("variable", signature("MeRIP.Peak"), function(object ){
  object@variate
})
#' @export
setMethod("variable<-", signature("MeRIP.Peak"), function(object ,value){
  object@variate <- value
  object
})

#' @title geneExpression
#' @description Extract gene level expression (RNAseq) data in normalized read counts
#' @param object The MeRIP object
#' @import DESeq2
#' @export
setMethod("geneExpression", signature("MeRIP.Peak"), function(object ){
  if(nrow(object@geneSum)> 0 ){
    return(object@geneSum)
  }else{
    input <- object@reads[,1:length(object@samplenames)]
    colnames(input) <- object@samplenames
    geneBins <- geneBins(object)
    ## Get input geneSum (gene level quantification)
    geneSum <- NULL
    for(i in 1:ncol(input) ){
      y <- input[,i]
      gene.sum <- tapply(y,geneBins$gene,sum)
      geneSum <- cbind(geneSum,gene.sum)
    }
    colnames(geneSum) <- object@samplenames
    
    size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
    
    norm.jointPeak_input <-t( t(jointPeak_input) / size.input )
    geneSum.norm <- t ( t(geneSum)/size.input)
    return(geneSum.norm)
  }
})

##################################################################################################################################
#' @title extractor for RNAseq data
#' @name geneExpression
#' @description Extract gene level expression (RNAseq) data in normalized read counts
#' @param object The MeRIP object
#' @import DESeq2
#' @export
setMethod("geneExpression", signature("MeRIP.Peak"), function(object ){
  if(nrow(object@geneSum)> 0 ){
    return(object@geneSum)
  }else{
    input <- object@reads[,1:length(object@samplenames)]
    colnames(input) <- object@samplenames
    geneBins <- geneBins(object)
    ## Get input geneSum (gene level quantification)
    geneSum <- NULL
    for(i in 1:ncol(input) ){
      y <- input[,i]
      gene.sum <- tapply(y,geneBins$gene,sum)
      geneSum <- cbind(geneSum,gene.sum)
    }
    colnames(geneSum) <- object@samplenames
    
    size.input <- DESeq2::estimateSizeFactorsForMatrix(geneSum)
    
    norm.jointPeak_input <-t( t(jointPeak_input) / size.input )
    geneSum.norm <- t ( t(geneSum)/size.input)
    return(geneSum.norm)
  }
})

#########################################################################################################################
#' @export
setMethod("jointPeak", signature("MeRIP.Peak"), function(object){
  if(nrow(object@jointPeaks)>0){
    return(object@jointPeaks)
  }else{
    stop("JointPeaks not reported, please call reportJointPeak() to reportJointPeak first.\n")
  }
})

#########################################################################################################################################
#' @export
setMethod("filter", signature("MeRIP.Peak"), function(object, i){
  
  if( object@jointPeak_threshold == 0 ){
    stop("joint peak not reported yet, cannot filter peaks!")
  }else if( nrow(object@jointPeaks) != nrow(object@jointPeak_id_pairs) ){
    stop("jointPeaks number doesn't match jointPeak_id_pairs!")
  }
  
  if( is.logical(i) & length(i) == nrow(object@jointPeaks) ){
    id <- which(i)
    
    if(nrow(object@jointPeak_ip)>0 ){
      if(nrow(object@jointPeak_ip) == nrow(object@jointPeak_input) & nrow(object@jointPeak_ip) ==  nrow(object@jointPeaks) ){
        object@jointPeak_ip <- object@jointPeak_ip[id,]
        object@jointPeak_input <- object@jointPeak_input[id,]
      }else{ stop("JointPeaks number doesn't match the joint peak input/ip read count dimension!")}
    }
    if(nrow(object@norm.jointPeak_ip)>0 ){
      if( nrow(object@norm.jointPeak_ip) == nrow(object@jointPeaks) ){
        object@norm.jointPeak_ip <- object@norm.jointPeak_ip[id,]
      }else{stop("Normalized jointPeak ip counts doesn't match the dimension of jointPeaks!")}
    }
    if(nrow(object@jointPeak_adjExpr)>0 ){
      if( nrow(object@jointPeak_adjExpr) == nrow(object@jointPeaks) ){
        object@jointPeak_adjExpr <- object@jointPeak_adjExpr[id,]
      }else{stop("Normalized jointPeak ip counts doesn't match the dimension of jointPeaks!")}
    }
    object@jointPeaks <-object@jointPeaks[id,]
    object@jointPeak_id_pairs <- object@jointPeak_id_pairs[id,]
  }else{
    stop("The filtering criteria needs to be the same length as number of peaks")
  }
  return(object)
})

############################################################################################################################

#' @export
#' @title subset MeRIPdata
#' @name select
#' @description subset dataset by samples.
#' @param object The MeRIP object
#' @param samples The samplenames to be subset or the index number of samples to be subset.
#' @return an MeRIP object of selected samples.
setMethod("select", signature("MeRIP"),function(object , samples ){
  
  id <- if( is.integer(samples) & all(samples > 0) & all(samples < length(samplenames(object))) ){
    samples
  }else if( all(samples %in% samplenames(object)) ){
    match(samples , samplenames(object))
  }
  
  newOb <- new(Class = "MeRIP",
               reads = counts(object)[,c(id,id + length(samplenames(object)))],
               binSize = object@binSize,
               geneModel = object@geneModel,
               gtf = object@gtf,
               bamPath.input = Input.files(object)[id],
               bamPath.ip = IP.files(object)[id],
               samplenames = samplenames(object)[id],
               geneBins = geneBins(object),
               GTF = object@GTF)
  
  newOb@geneSum <- if( nrow(object@geneSum) > 1 ){object@geneSum[,id]}else{object@geneSum}
})


#' @export
#' @name select
#' @description subset dataset by samples.
#' @param object The MeRIP.Peak object
#' @param samples The samplenames to be subset or the index number of samples to be subset.
#' @return an MeRIP.Peak object of selected samples.
setMethod("select", signature("MeRIP.Peak"),function(object , samples ){
  
  id <- if( is.integer(samples) & all(samples > 0) & all(samples < length(samplenames(object))) ){
    samples
  }else if( all(samples %in% samplenames(object)) ){
    match(samples , samplenames(object))
  }else{
    stop("Please specify valide samplenames or index to subset the dataset.")
  }
  
  newOb <- new(Class = "MeRIP.Peak",
               reads = counts(object)[,c(id,id + length(samplenames(object)))],
               binSize = object@binSize,
               geneModel = object@geneModel,
               gtf = object@gtf,
               bamPath.input = Input.files(object)[id],
               bamPath.ip = IP.files(object)[id],
               samplenames = samplenames(object)[id],
               geneBins = geneBins(object),
               GTF = object@GTF,
               peakCalling = object@peakCalling,
               jointPeak_threshold = object@jointPeak_threshold,
               test.method = object@test.method ,
               jointPeak_id_pairs = object@jointPeak_id_pairs,
               jointPeaks = object@jointPeaks
  )
  
  newOb@geneSum <- if( nrow(object@geneSum) > 1 ){object@geneSum[,id]}else{object@geneSum}
  newOb@peakCallResult <- if(nrow(object@peakCallResult) > 1){object@peakCallResult[,id]}else{object@peakCallResult}
  newOb@jointPeak_ip <- if(nrow(object@jointPeak_ip) > 1 ){object@jointPeak_ip[,id]}else{object@jointPeak_ip}
  newOb@jointPeak_input <- if(nrow(object@jointPeak_input) > 1 ){object@jointPeak_ip[,id]}else{object@jointPeak_input}
  newOb@norm.jointPeak_ip <- if(nrow(object@norm.jointPeak_ip) > 1 ){object@norm.jointPeak_ip[,id]}else{object@norm.jointPeak_ip}
  newOb@sizeFactor <- if(nrow(object@sizeFactor) > 1 ){object@sizeFactor[id,]}else{object@sizeFactor}
  newOb@variate <- if(nrow(object@variate) > 1 ){ if(ncol(object@variate)>1){
    object@variate[id,]
  }else{
    newVar <- data.frame(object@variate[id,])
    colnames(newVar) <- colnames(object@variate)
    newVar
  } }else{
    object@variate
    }
  newOb@jointPeak_adjExpr <- if(nrow(object@jointPeak_adjExpr) > 1 ){object@jointPeak_adjExpr[,id]}else{object@jointPeak_adjExpr}
  if(nrow(object@test.est) > 1 ){cat("Inferential test is not inherited because test result changes when samples are altered!\nPlease re-do test.\n")}
  return(newOb)
})


#############################################################################################################################

#' @export
#' @name results
#' @title export results
#' @description The extractor for final test result.
#' @param object The MeRIP.Peak object.
#' @return joint peaks (with test result) in a data.frame.
setMethod("results", signature("MeRIP.Peak"), function(object){
  if(object@test.method != "none"){
    return(cbind(jointPeak(object),
                 as.data.frame(object@test.est)
    )
    )
  }else if(nrow(object@jointPeaks)> 1){
    cat("No inferential test has been performed. Returning jointPeaks!\n")
    return( jointPeak(object) )
  }else{
    stop("There is no result to be returned.")
  }
})

