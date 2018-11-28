#' @title QTL_BetaBin
#' @param MeRIPdata The MeRIP.Peak object
#' @param vcf_file The vcf file for genotype. The chromosome position must be sorted!!
#' @param BSgenome The BSgenome object. This needs to match the genome version of the gtf files. 
#' @param testWindow Integer. Test SNPs in <testWindow> bp window flanking the peak.
#' @param Chromosome The chromsome to run QTL test.
#' @param Range The position range on a chromosome to test.
#' @param Covariates The matrix for covariates to be included in the test.
#' @param AdjustGC Logic. Choose whether explicitly adjust GC bias.
#' @param AdjIPeffi Logic. Choose whether explicitly adjust overall IP efficiency
#' @param normalizeGenotype Logic. Choose whether genotype is normalized to mean = 0, var = 1 before regression.
#' @import stringr
#' @import vcfR
#' @import BSgenome
#' @import gamlss
#' @import gamlss.dist
#' @import broom
#' @export
QTL_BetaBin <- function( MeRIPdata , vcf_file, BSgenome = BSgenome.Hsapiens.UCSC.hg19,testWindow = 100000, Chromosome, Range = NULL, Covariates = NULL, AdjustGC = TRUE, AdjIPeffi = TRUE , PCsToInclude = 0 , normalizeGenotype = FALSE, thread = 1 ){
  
  ##check input
  if( !is(MeRIPdata, "MeRIP.Peak") ){
    stop("The input MeRIPdata needs to be an MeRIP.Peak object!")
  }else if( !nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_ip) & nrow(MeRIPdata@jointPeaks) == nrow(MeRIPdata@jointPeak_input) ){
    stop("The peak counts matrix dimension must match the dimension of jointPeaks!")
  }
  
  ## check samples in genotype files and in MeRIP.Peak object
  tmpVcf <- tempfile(fileext = ".vcf.gz")
  system(paste0("zcat ",vcf_file," | awk 'NR==1 {print $0} (/^#CHROM/){print $0}'| gzip > ",tmpVcf))
  tmp.vcf <-try(  read.vcfR( file =tmpVcf, verbose = F ) , silent = T)
  ############################################################################################################################################################
  ### This is to handle a wired error in the read.vcfR function. #############################################################################################
  if(class(tmp.vcf) == "try-error"){                                                                                                                       ##
    system(paste0("zcat ",vcf_file," | awk '/^#/ {print $0}'| gzip >",tmpVcf))  ##
    tmp.vcf <-read.vcfR( file =tmpVcf, verbose = F )                                                                                                       ##
  }                                                                                                                                                         ##
  ############################################################################################################################################################
  unlink(tmpVcf) # remove the temp file to free space.
  genotypeSamples <- colnames(tmp.vcf@gt)[-c(1)]
  if(length(intersect(genotypeSamples,samplenames(MeRIPdata))) == 0 ){
    stop("The samplenames must match in VCF file and in MeRIP.Peak object! We found no overlap between sample names in these two files!")
  }else if(length(intersect(genotypeSamples,samplenames(MeRIPdata))) != length(samplenames(MeRIPdata) )){
    cat("The samples in the VCF don't totally match samples in the MeRIP.Peak object; ")
    cat("Only samples in MeRIP.Peak object overlapping samples in VCF file will be analyzed in QTL mapping!\nSubsetting samples...\n")
    MeRIPdata <- select(MeRIPdata, intersect(genotypeSamples,samplenames(MeRIPdata))  )
    cat(paste0(paste(intersect(genotypeSamples,samplenames(MeRIPdata)),collapse = " "), "\n(",length(intersect(genotypeSamples,samplenames(MeRIPdata))),") samples will be analyzed!"))
  }else{
    ## make sure the order of samples aligned between phenotype and genotype
    MeRIPdata <- select(MeRIPdata, intersect(genotypeSamples,samplenames(MeRIPdata))  )
  }
  
  ### Preprocess
  ##fitler out peaks with zero count
  MeRIPdata <- filter(MeRIPdata, !apply(extractInput(MeRIPdata), 1, function(x) any(x == 0 )) )
  cat("Peaks with zero read count in input data have been removed.\n")
  
  T0 <- colSums(counts(MeRIPdata)[,1:length(MeRIPdata@samplenames)] )
  T1 <- colSums(counts(MeRIPdata)[,(length(MeRIPdata@samplenames)+1) : (2*length(MeRIPdata@samplenames)) ] )
  ##filter out peaks with OR < 1
  enrichFlag <- apply( t( t(extractIP(MeRIPdata))/T1 )/ t( t( extractInput(MeRIPdata) )/T0 ),1,function(x){sum(x>1)> MeRIPdata@jointPeak_threshold})
  MeRIPdata <-  filter(MeRIPdata, enrichFlag )
  cat(paste0("Peaks with odd ratio > 1 in more than ",MeRIPdata@jointPeak_threshold," samples will be retained.\n",nrow(jointPeak(MeRIPdata))," peaks remaining for QTL mapping.\n"))
  ## estimate IP efficiency
  OR <- t( apply(extractIP(MeRIPdata),1,.noZero)/T1 )/ t( t( extractInput(MeRIPdata) )/T0 )
  colnames(OR) <- MeRIPdata@samplenames
  OR.id <- which( rowMeans(OR) < quantile( rowMeans(OR), 0.95 ) & rowMeans(OR) > quantile( rowMeans(OR), 0.1) )# remove two tails
  K_IPe <- log( apply(OR[OR.id,], 2, function(x){coef(lm(x~rowMeans(OR)[OR.id]) )[2]}) ) # fit linear moedel to obtain the IP efficiency offset
  
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
    if(AdjIPeffi){y <- (log( OR ) - K_IPe)}else{y <-  log(OR) }
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
  ##Principal components
  if(PCsToInclude > 0 & PCsToInclude <= length(MeRIPdata@samplenames) ){
    cat("Computing Principal components.\n")
    if(AdjustGC & AdjIPeffi){
      PCs <- prcomp(t( (log(OR ) -  K_IPe - Fij)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),] ) )$x
    }else if( AdjustGC & !AdjIPeffi){
      PCs <- prcomp(t( (log(OR ) - Fij)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),] ))$x
    }else if(AdjIPeffi){
      PCs <- prcomp(t( (log(OR ) - K_IPe)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),] ))$x
    }else{
      PCs <- prcomp( t( log(OR)[apply(extractIP(MeRIPdata),1,function(x) all(x!=0)),] ))$x
    }
  }else if(PCsToInclude > length(MeRIPdata@samplenames) ){
    stop("The number of PCs needs to be no larger than the sample size!")
  }
 
  cat("The following message can be ignored...\n----------------\n")
  ## set ranges on the chromosome that can be tested
  con1 <- pipe(paste0("zcat ",vcf_file," | grep -w '",Chromosome,"' |awk '!/^#/ {print $2}' | head -n1"))
  con2 <- pipe(paste0("zcat ",vcf_file," | grep -w '",Chromosome,"' |awk '!/^#/ {print $2}' | tail -n1"))
  vcfRange <- c( scan( con1 , quiet = T ) ,
                scan( con2 , quiet = T ))
  close(con1)
  close(con2)
  cat("----------------\n")
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
  Y1 <- extractIP(MeRIPdata)[test.id,]
  Y0 <- extractInput(MeRIPdata)[test.id,]
  if(AdjustGC){FIj <- Fij[test.id,] }
  
  ## ditermine study design according to parameters
  variables <-  "offset(log(T1/T0))"
  if( AdjIPeffi ){ variables <- paste(variables, "offset(K_IPe)", sep  = " +") }
  if( AdjustGC ){ variables <- paste(variables, "offset(Fj)", sep = " +") } 
  if(! is.null(Covariates)  ){ 
    colnames(Covariates)
    variables <-  paste(variables, paste(colnames(Covariates),collapse = " + "), sep = "+")
  }
  if( PCsToInclude > 0 ){
    variables <-  paste(variables, paste("PC", 1:PCsToInclude, sep = "",collapse = " + "), sep = "+")
  }
  design <- formula( paste0("cbind(Y1i , Y0i) ~" ,"  G + ", variables) )
  
  cat(paste0("Start beta-binomial regression for ",length(peak_bed.gr)," peaks and SNPs in ",round(testWindow/1000,digits = 1),"kb flanking each peaks on chromosome ",Chromosome,".\n"))
  if(AdjustGC){cat("Will correct sample specific GC bias\n")}
  ## test each peak
  startTime <- Sys.time()
  registerDoParallel(thread)
  testResult <- foreach( i = 1:length(peak_bed.gr), .combine = rbind  )%dopar% {

    ## get the range where SNPs are available
    testRange <- GenomicRanges::intersect(IRanges(start(peak_bed.gr[i])-testWindow, end(peak_bed.gr[i])+testWindow ),vcfRange)

    ## Test association if there is SNP available for this peak
    if(length(testRange)==1){
      
      ## Use unix command line to accesss the genotype from vcf file. This is for fast data access.
      tmpVcf <- tempfile(fileext = ".vcf.gz")
      system(paste0("zcat ",vcf_file," | awk 'NR==1 {print $0} (/^#CHROM/){print $0}(!/^#/ && $2 > ",start(testRange)," && $2 < ",end(testRange)," ) {print $0}'| gzip > ",tmpVcf))
      geno.vcf <-try(  read.vcfR( file =tmpVcf, verbose = F ) , silent = T)
      ############################################################################################################################################################
      ### This is to handle a wired error in the read.vcfR function. #############################################################################################
      if(class(geno.vcf) == "try-error"){                                                                                                                       ##
        system(paste0("zcat ",vcf_file," | awk '/^#/ {print $0} (!/^#/ && $2 > ",start(testRange)," && $2 < ",end(testRange)," ) {print $0}'| gzip >",tmpVcf))  ##
        geno.vcf <-read.vcfR( file =tmpVcf, verbose = F )                                                                                                       ##
      }                                                                                                                                                         ##
      ############################################################################################################################################################
      unlink(tmpVcf) # remove the temp file to free space.
      
      ## check if genotype available for this peak
      if( nrow(geno.vcf@fix) == 0 ){ return(NULL) }
      
      ## Get genotype as dosage format and filter for MAF
      if( unique(geno.vcf@gt[,"FORMAT"]) == "DS" ){
        ## Directly extract dosage
        geno <- if(nrow(geno.vcf@fix) == 1){ 
          t(apply(extract.gt(geno.vcf, element = 'DS' ),2,as.numeric ) )
        }else{
          apply(extract.gt(geno.vcf, element = 'DS' ),2,as.numeric )
          }
        rownames(geno) <- geno.vcf@fix[,"ID"]
        ## filter out any genotype that has MAF<0.05
        MAF <- apply(geno,1,function(x) !any(table(round(x) )>0.95*ncol(geno)) )
        geno <- geno[MAF,] 
        geno.vcf <- geno.vcf[MAF,]
      }else if(unique(geno.vcf@gt[,"FORMAT"]) == "GT"){
        geno.vcf <- geno.vcf[is.biallelic(geno.vcf),]
        ## get genotype as Dosage
        tmp_geno <- extract.gt(geno.vcf, element = 'GT' )
        geno <- t( apply( tmp_geno ,1, .genoDosage ) )
        colnames(geno) <- colnames(tmp_geno)
        ## filter out any genotype that has MAF<0.05
        MAF <- apply(geno,1,function(x) !any(table(x)>0.95*ncol(geno)) )
        geno <- geno[MAF,] 
        geno.vcf <- geno.vcf[MAF,]
      }
      if( nrow(geno) == 0 ){ return(NULL) } ## skip this iteration if no genotype left after filtering
      
      ## normalized genotype is necessary
      if(normalizeGenotype){
        geno <- t( apply(geno,1,function(x){ (x - mean(x) )/sd(x) }) )
      }
      
      if(AdjustGC){Fj <- FIj[i,]}
      
      
      ## Test QTLs for peak.j
      tmp_est <- as.data.frame(matrix(nrow = nrow(geno),ncol = 5),row.names = rownames(geno) )
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
                                     std.err = est[est$term == "G","std.error"],
                                     pvalue = est[est$term == "G","p.value"], 
                                     theta = 1/exp(est[est$parameter == "sigma","estimate"]),
                                     p.theta = est[est$parameter == "sigma","p.value"] ) 
        }else{
            tmp_est[ii,] <- data.frame(beta = NA, std.err = NA,  pvalue =NA,theta = NA, p.theta = NA ) 
          }
        
      }
      colnames(tmp_est) <- c("beta","std.err","pvalue","theta","p.theta")
        
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
  cat(paste("Time used to test association: ",round(difftime(endTime, startTime, units = "mins"),digits = 1)," mins. \n"))
  cat(paste0(nrow(testResult)," SNP-peak pair tested.\n"))
  return(testResult)
}


## Helper function to convert genotype into dosage. 
.genoDosage <- function(x){
  return( stringr::str_count(x,"1") )
}

