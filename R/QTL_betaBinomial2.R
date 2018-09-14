#' @title QTL_BetaBin2
#' @param MeRIPdata The MeRIP.Peak object
#' @param vcf_file The vcf file for genotype. The chromosome position must be sorted!!
#' @param BSgenome The BSgenome object. This needs to match the genome version of the gtf files. 
#' @param testWindow Integer. Test SNPs in <testWindow> bp window flanking the peak.
#' @param Chromosome The chromsome to run QTL test.
#' @param Range The position range on a chromosome to test.
#' @param Covariates The matrix for covariates to be included in the test.
#' @import stringr
#' @import vcfR
#' @import BSgenome
#' @import gamlss
#' @import gamlss.dist
#' @import broom
#' @export
QTL_BetaBin2 <- function( MeRIPdata , vcf_file, BSgenome = BSgenome.Hsapiens.UCSC.hg19,testWindow = 100000, Chromosome, Range = NULL, Covariates = NULL, AdjustGC = TRUE, PCsToInclude = 0 , thread = 1 ){
   
  ##check input
  if( !is(MeRIPdata, "MeRIP.Peak") ){
    stop("The input MeRIPdata needs to be an MeRIP.Peak object!")
  }else if( !nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_ip) & nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_input) ){
    stop("The peak counts matrix dimension must match the dimension of jointPeaks!")
  }else if( nrow(MeRIPdata@sizeFactor) == 0){
    cat("The input data has not performed library normalization. Normalizing library...\n")
    MeRIPdata <- normalizeLibrary(MeRIPdata)
  }
  ### Preprocess
  ##fitler out peaks with zero count
  MeRIPdata <- filter(MeRIPdata, !apply(extractInput(MeRIPdata), 1, function(x) any(x == 0 )) )
  
  
   ##filter out peaks with OR  < 1
  enrichFlag <- apply(  t( t(extractIP(MeRIPdata))/MeRIPdata@sizeFactor$ip ) /  t( t( extractInput(MeRIPdata) )/MeRIPdata@sizeFactor$input ) ,1,function(x){sum(x>1)> MeRIPdata@jointPeak_threshold})
  MeRIPdata <-  filter(MeRIPdata, enrichFlag )
  
  norm.IP <- t( apply(extractIP(MeRIPdata),1,.noZero)/MeRIPdata@sizeFactor$ip )
  norm.Input <- t( t( extractInput(MeRIPdata) )/MeRIPdata@sizeFactor$input )
  
  ## estimate IP efficiency
  OR <- norm.IP/norm.Input
  colnames(OR) <- MeRIPdata@samplenames
  
  if(AdjustGC){
    cat("Computing GC content for peaks...\n")
    ## GC content correction
    peak.gr <- .peakToGRangesList( jointPeak(MeRIPdata))
    registerDoParallel( thread )
    peakGC <- foreach( i = 1:length(peak.gr), .combine = c)%dopar%{
      peakSeq <- paste( getSeq( BSgenome , peak.gr[[i]] ,as.character =T ) , collapse = "")
      sum( str_count(peakSeq, c("G","g","C","c"))  )/nchar(peakSeq)
    }
    peakGC_l <- round(peakGC,digits = 2)
    peakGC_l[which(peakGC_l<0.2)] <-median(peakGC_l[which(peakGC_l<0.2)] ) # combine some bins at low GC due to low number of peaks
    peakGC_l[which(peakGC_l>0.84)] <-median(peakGC_l[which(peakGC_l>0.84)] ) # combine some bins at high GC due to low number of peaks
    l <- sort(unique(peakGC_l))
    if(AdjIPeffi){y <- (log( OR ) - K_IPe)}else{y <-  log(OR) } 
    colnames(y) <- MeRIPdata@samplenames
    b.l <- tapply(rowMeans( y ), peakGC_l ,  median) 
    bil <- apply( y, 2, tapply, peakGC_l, median )
    bi. <- apply(  y , 2, median )
    b.. <- median(  y )
    Fil <- as.data.frame( as.matrix(bil) - as.vector(b.l) ) - ( bi. - b.. )
    Fij <- foreach( ii = 1:length(MeRIPdata@samplenames), .combine = cbind)%do%{
      GC_fit <- lm(Fil[,ii] ~ poly(l,4) )
      predict(GC_fit, newdata = data.frame(l = peakGC) )
    }
    colnames(Fij) <- MeRIPdata@samplenames
  }
  ##Principal components
  if(PCsToInclude > 0 & PCsToInclude <= length(MeRIPdata@samplenames) ){
    cat("Computing Principal components...\n")
    if( AdjustGC ){
      PCs <- prcomp(t( (log(OR) - Fij)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),] ))$x
    }else{
      PCs <- prcomp(t( log(OR)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),]  ))$x
    }
  }else if(PCsToInclude > length(MeRIPdata@samplenames) ){
    stop("The number of PCs needs to be no larger than the sample size!")
  }
 
  
  ## set ranges on the chromosome that can be tested
  con <- pipe(paste0("zcat ",vcf_file," | awk '!/^#/ {print $2}' | tail -n1"))
  vcfRange <- c( read.table(gzfile(vcf_file), nrows = 1)[,2] ,
                scan( con , quiet = T ))
  close(con)
  ## update test Range if necessary
  if(!is.null(Range)){
    vcfRange <- intersect(IRanges(vcfRange[1],vcfRange[2]),IRanges(Range[1],Range[2]) )
  }else{
    vcfRange <- IRanges(vcfRange[1],vcfRange[2])
  }

  ## parse bed12 file
  peak_bed <- jointPeak(MeRIPdata)
  peak_bed.gr <- makeGRangesFromDataFrame(peak_bed, keep.extra.columns = T)
  test.id <- which(peak_bed$chr == Chromosome & (( peak_bed$end+testWindow )  > start(vcfRange) ) & (( peak_bed$start-testWindow )  < end(vcfRange) ) )
  
  peak_bed.gr <- peak_bed.gr[test.id ]
  Y1 <- norm.IP[test.id,]
  Y0 <- norm.Input[test.id,]
  if(AdjustGC){FIj <- Fij[test.id,] }
  
  ## ditermine study design according to parameters
  variables <-  ""
  if( AdjustGC ){ variables <- paste(variables, "offset(Fj)", sep = " +") } 
  if(! is.null(Covariates)  ){ 
    colnames(Covariates)
    variables <-  paste(variables, paste(colnames(Covariates),collapse = " + "), sep = "+")
  }
  if( PCsToInclude > 0 ){
    variables <-  paste(variables, paste("PC", 1:PCsToInclude, sep = "",collapse = " + "), sep = "+")
  }
  design <- formula( paste0("cbind(Y1i , Y0i) ~" , variables," + G ") )
  
  cat(paste0("Start beta-binomial test for ",length(peak_bed.gr)," peaks and SNPs in ",round(testWindow/1000,digits = 1),"kb flanking each peaks.\n"))
  ## test each peak
  startTime <- Sys.time()
  registerDoParallel(thread)
  testResult <- foreach( i = 1:length(peak_bed.gr), .combine = rbind  )%dopar% {

    ## get the range where SNPs are available
    testRange <- GenomicRanges::intersect(IRanges(start(peak_bed.gr[i])-testWindow, end(peak_bed.gr[i])+testWindow ),vcfRange)

    ## Test association if there is SNP available for this peak
    if(length(testRange)==1){
      ## read genotype
      system(paste0("zcat ",vcf_file," | awk 'NR==1 {print $0} (/^#CHROM/){print $0}(!/^#/ && $2 > ",start(testRange)," && $2 < ",end(testRange)," ) {print $0}'| gzip > ~/tmp",i,".vcf.gz"))
      geno.vcf <-try(  read.vcfR( file =paste0("~/tmp",i,".vcf.gz"), verbose = F ) , silent = T)
      ####################################################################
      ### This is to handle a wired error in the read.vcfR function.
      if(class(geno.vcf) == "try-error"){
        system(paste0("zcat ",vcf_file," | awk '/^#/ {print $0} (!/^#/ && $2 > ",start(testRange)," && $2 < ",end(testRange)," ) {print $0}'| gzip > ~/tmp",i,".vcf.gz"))
        geno.vcf <-read.vcfR( file =paste0("~/tmp",i,".vcf.gz"), verbose = F )
      }
      ####################################################################
      file.remove(paste0("~/tmp",i,".vcf.gz"))
      ## filter biallelic snps
      geno.vcf <- geno.vcf[is.biallelic(geno.vcf),]

      ## get genotype as Dosage
      tmp_geno <- extract.gt(geno.vcf, element = 'GT' )
      geno <- t( apply( tmp_geno ,1, .genoDosage ) )
      colnames(geno) <- colnames(tmp_geno)
      ## filter out any genotype that has MAF<0.05
      MAF <- apply(geno,1,function(x) !any(table(x)>0.95*ncol(geno)) )
      geno <- geno[MAF,] 
      geno.vcf <- geno.vcf[MAF,]
      
      if(AdjustGC){Fj <- FIj[i,]}
      
      tmp_est <- as.data.frame(matrix(nrow = nrow(geno),ncol = 4),row.names = rownames(geno) )
      for( ii in 1:nrow(geno) ){
        if(AdjustGC){fit_data <- data.frame(Y0i = Y0[i,], Y1i = Y1[i,], T1, T0, K_IPe, Fj , G = geno[ii,])}else{fit_data <- data.frame(Y0i = Y0[i,], Y1i = Y1[i,], T1, T0, K_IPe , G = geno[ii,])}
        ## add PCs to data
        if(PCsToInclude > 1 ){ fit_data <- cbind(fit_data,data.frame(PCs[,1:PCsToInclude]) )}else if(PCsToInclude == 1){fit_data <- cbind(fit_data,data.frame(PC1 = PCs[,1]) ) }
        ## add covariates
        if(! is.null(Covariates) ){fit_data <- cbind(fit_data,as.data.frame(Covariates) )  }
        
        fit <- try( gamlss( design ,data = fit_data , family = BB(mu.link = "logit")) )
        if(class(fit)[1]!= "try-error"){
          est <- tidy(fit)
          tmp_est[ii,] <- data.frame(beta =  est[est$term == "G","estimate"], 
                                     theta = 1/exp(est[est$parameter == "sigma","estimate"]), 
                                     pvalue = est[est$term == "G","p.value"], 
                                     p.theta = est[est$parameter == "sigma","p.value"] )  
        }else{
            tmp_est[ii,] <- data.frame(beta = NA, theta = NA, pvalue =NA, p.theta = NA ) 
          }
        
      }
      colnames(tmp_est) <- c("beta","theta","pvalue","p.theta")
        
      ## calculate distance with respect to transcript(gene) strand
      distance <- if(as.character(strand(peak_bed.gr[i])) == "+"){
        as.integer(geno.vcf@fix[,"POS"])- round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2)
      }else{
        round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2) - as.integer(geno.vcf@fix[,"POS"])
      }
      report <- cbind(data.frame(
        SNP = paste(geno.vcf@fix[,"CHROM"],geno.vcf@fix[,"POS"],sep = ":"),
        SNPID = rownames(tmp_est),
        REF = geno.vcf@fix[,"REF"],
        ALT = geno.vcf@fix[,"ALT"],
        PEAK = paste0(Chromosome,":",peak_bed.gr[i]$thickStart,"-",peak_bed.gr[i]$thickEnd,"_",peak_bed.gr[i]$name,"_",strand( peak_bed.gr[i]) ),
        DISTANCE = distance
      ),tmp_est)
      report[!is.na(report$beta),]
    }

  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  endTime <- Sys.time()
  cat(paste("Time used to test association: ",difftime(endTime, startTime, units = "mins")," mins... \n"))

  return(testResult)
}


## Helper function to convert genotype into dosage. 
.genoDosage <- function(x){
  return( stringr::str_count(x,"1") )
}

