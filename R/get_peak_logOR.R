#' @title get_peak_logOR
#' @description Compute peak logOR, adjust GC bias and IP efficiency for m6A QTL analysis
#' @param MeRIPdata The MeRIP.Peak object containing peak calling result.
#' @param vcf_file The vcf file for genotype. The chromosome position must be sorted!!
#' @param BSgenome The BSgenome object. This needs to match the genome version of the gtf files. 
#' @param AdjustGC Logic. Choose whether explicitly adjust GC bias (default: TRUE). 
#' @param AdjIPeffi Logic. Choose whether explicitly adjust overall IP efficiency (default: TRUE). 
#' @import vcfR
#' @import BSgenome
#' @export
get_peak_logOR <- function( MeRIPdata, vcf_file = NULL, BSgenome = BSgenome.Hsapiens.UCSC.hg19, AdjustGC = TRUE, AdjIPeffi = TRUE, thread = 1 ){
  
  ## check input
  if( !is(MeRIPdata, "MeRIP.Peak") ){
    stop("The input MeRIPdata needs to be an MeRIP.Peak object!")
  }else if( !nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_ip) & nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_input) ){
    stop("The peak counts matrix dimension must match the dimension of jointPeaks!")
  }
  
  ## select samples matching the vcf file
  if(!is.null(vcf_file)){
    ## check samples in genotype files and in MeRIP.Peak object
    tmpVcf <- tempfile(fileext = ".vcf.gz")
    system(paste0("zcat ",vcf_file," | awk 'NR==1 {print $0} (/^#CHROM/){print $0}'| gzip > ",tmpVcf))
    tmp.vcf <- try(read.vcfR( file =tmpVcf, verbose = F ) , silent = T)
    ##########################################################################################################
    ### This is to handle a wired error in the read.vcfR function. ###########################################
    if(class(tmp.vcf) == "try-error"){                                                                      ##
      system(paste0("zcat ",vcf_file," | awk '/^#/ {print $0}'| gzip >",tmpVcf))                            ##
      tmp.vcf <-read.vcfR( file =tmpVcf, verbose = F )                                                      ##
    }                                                                                                       ##
    ##########################################################################################################
    unlink(tmpVcf) # remove the temp file to free space.
    genotypeSamples <- colnames(tmp.vcf@gt)[-c(1)]
    
    if(length(intersect(genotypeSamples,samplenames(MeRIPdata))) == 0 ){
      stop("The samplenames must match in VCF file and in MeRIP.Peak object! We found no overlap between sample names in these two files!")
    }else if(length(intersect(genotypeSamples,samplenames(MeRIPdata))) != length(samplenames(MeRIPdata) )){
      cat("The samples in the VCF don't totally match samples in the MeRIP.Peak object; ")
      cat("Only samples in MeRIP.Peak object overlapping samples in VCF file will be analyzed in QTL mapping!\nSubsetting samples...\n")
      MeRIPdata <- MeRIPtools::select(MeRIPdata, intersect(genotypeSamples,samplenames(MeRIPdata))  )
      cat(paste0(paste(intersect(genotypeSamples,samplenames(MeRIPdata)),collapse = " "), "\n(",length(intersect(genotypeSamples,samplenames(MeRIPdata))),") samples will be analyzed!"))
    }else{
      ## make sure the order of samples aligned between phenotype and genotype
      MeRIPdata <- select(MeRIPdata, intersect(genotypeSamples,samplenames(MeRIPdata))  )
    }
    
  }else{
    cat("No VCF file input. Use all samples in the MeRIPdata data. ")
  }
  
  
  ### Preprocess
  ## Filter out peaks with zero count
  MeRIPdata <- filter(MeRIPdata, !apply(extractInput(MeRIPdata), 1, function(x) any(x == 0 )) )
  cat("Peaks with zero read count in input data have been removed.\n")
  
  T0 <- colSums(counts(MeRIPdata)[,1:length(MeRIPdata@samplenames)] )
  T1 <- colSums(counts(MeRIPdata)[,(length(MeRIPdata@samplenames)+1) : (2*length(MeRIPdata@samplenames)) ] )
  
  ## Filter out peaks with OR < 1
  enrichFlag <- apply( t( t(extractIP(MeRIPdata))/T1 )/ t( t( extractInput(MeRIPdata) )/T0 ),1,function(x){sum(x>1)> MeRIPdata@jointPeak_threshold})
  MeRIPdata <-  filter(MeRIPdata, enrichFlag )
  cat(paste0("Peaks with odd ratio > 1 in more than ",MeRIPdata@jointPeak_threshold," samples will be retained.\n",nrow(jointPeak(MeRIPdata))," peaks remaining for QTL mapping.\n"))
  
  ## Estimate IP efficiency
  OR <- t( apply(extractIP(MeRIPdata),1,.noZero)/T1 )/ t( t( extractInput(MeRIPdata) )/T0 )
  colnames(OR) <- MeRIPdata@samplenames
  logOR <- log(OR)
  
  logOR.id <- which( rowMeans(logOR) < quantile( rowMeans(logOR), 0.95 ) & rowMeans(logOR) > quantile( rowMeans(logOR), 0.05) )# remove two tails
  K_IPe_ij <- apply(logOR[logOR.id,], 2, function(x){
    
    fit <- lm(y~m, data = data.frame(y = x, m=rowMeans(logOR)[logOR.id] ))
    y.est <- predict(fit, newdata =  data.frame(m = rowMeans(logOR)))
    return( y.est - rowMeans(logOR) )
  })
  
  ## Estimate GC bias offset
  if(AdjustGC){
    cat("Computing GC content for peaks\n")
    ## GC content correction
    peak.gr <- .peakToGRangesList( jointPeak(MeRIPdata))
    cat("...")
    
    registerDoParallel( thread )
    peakSeq <- foreach( i = 1:length(peak.gr), .combine = c )%dopar%{
      paste( getSeq( BSgenome , peak.gr[[i]] ,as.character =T ) , collapse = "")
    }
    peakGC <- sapply( peakSeq, function(x){ sum( str_count(x, c("G","g","C","c"))  )/nchar(x) } )
    
    cat("...")
    
    peakGC_l <- round(peakGC,digits = 2)
    peakGC_l[which(peakGC_l<0.2)] <-median(peakGC_l[which(peakGC_l<0.2)] ) # combine some bins at low GC due to low number of peaks
    peakGC_l[which(peakGC_l>0.84)] <-median(peakGC_l[which(peakGC_l>0.84)] ) # combine some bins at high GC due to low number of peaks
    l <- sort(unique(peakGC_l))
    
    if(AdjIPeffi){ 
      y <- log(OR) - K_IPe_ij
    }else{
      y <- log(OR)
    }
    colnames(y) <- MeRIPdata@samplenames
    
    b.l <- tapply(rowMeans( y ), peakGC_l ,  median) 
    bil <- apply( y, 2, tapply, peakGC_l, median )
    bi. <- apply(  y , 2, median )
    b.. <- median(  y )
    Fil <- as.data.frame( as.matrix(bil) - as.vector(b.l) ) - ( bi. - b.. )
    Fij <- foreach( ii = 1:length(MeRIPdata@samplenames), .combine = cbind)%dopar%{
      GC_fit <- lm(Fil[,ii] ~ poly(l,4) )
      predict(GC_fit, newdata = data.frame(l = peakGC) )
    }
    colnames(Fij) <- MeRIPdata@samplenames
    cat("...\n")
  }
  
  ## Adjust logOR with GC and IP efficiency
  if(AdjustGC & AdjIPeffi){
    cat("Adjust logOR with GC and IP efficiency ... \n")
    logOR_adjusted <- log(OR) -  K_IPe_ij - Fij
  }else if( AdjustGC & !AdjIPeffi){
    cat("Adjust logOR with GC ... \n")
    logOR_adjusted <- log(OR) - Fij
  }else if(AdjIPeffi){
    cat("Adjust logOR with IP efficiency ... \n")
    logOR_adjusted <- log(OR) - K_IPe_ij
  }else{
    logOR_adjusted <- log(OR)
  }
  
  ## parse bed12 file to get peak coordinates
  peak_bed <- jointPeak(MeRIPdata)
  peak_bed$PEAK <- paste0(peak_bed$chr,":", peak_bed$start,"-",peak_bed$end,"_",peak_bed$name,"_",peak_bed$strand )
  
  peak_logOR <- data.frame(peak_bed[, c("chr", "start", "end", "PEAK")], logOR_adjusted)
  
  return(list(logOR = peak_logOR, MeRIPdata = MeRIPdata))
  
}


