#' @title QTL_PoissonGamma
#' @param pheno The phenotype data matrix. Needs to be IP read count that has been normalized for expression level.
#' @param vcf_file The vcf file for genotype. The chromosome position must be sorted!!
#' @param peak_bed The peak file in BED12 format that needs to correspond to phenotype data matrix.
#' @param testWindow Integer. Test SNPs in <testWindow> bp window flanking the peak.
#' @param Chromosome The chromsome to run QTL test.
#' @param Range The position range on a chromosome to test.
#' @param Covariates The matrix for covariates to be included in the test.
#' @param maxPsi The max estimation for the random effect parameter Psi.
#' @import stringr
#' @import vcfR
#' @export
QTL_PoissonGamma <- function( pheno, vcf_file, peak_bed, testWindow = 100000, Chromosome, Range = NULL, Covariates = NULL, maxPsi = 100, thread = 1 ){

  ##check input
  if(nrow(pheno) != nrow(peak_bed) ){
    stop("The number of row of phenotype is not equal to the number of peaks!")
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
  colnames(peak_bed) <- c("chr","start","end","name","score","strand","thickStart","thickEnd","RGB","numBlock","blockSize","blockStart")
  peak_bed.gr <- makeGRangesFromDataFrame(peak_bed, keep.extra.columns = T)
  test.id <- which(peak_bed$chr == Chromosome & (( peak_bed$end+testWindow )  > start(vcfRange) ) & (( peak_bed$start-testWindow )  < end(vcfRange) ) )
  peak_bed.gr <- peak_bed.gr[test.id ]
  phenoY <- pheno[test.id,]

  ## test each peak
  startTime <- Sys.time()
  registerDoParallel(thread)
  testResult <- foreach( i = 1:length(peak_bed.gr), .combine = rbind  )%dopar% {

    ## get the range where SNPs are available
    testRange <- intersect(IRanges(start(peak_bed.gr[i])-testWindow, end(peak_bed.gr[i])+testWindow ),vcfRange)

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

      ## Determine whether to include covariates
      if( is.null(Covariates) ){
        Y <- round( phenoY[i,] )
        psi <- 10
        tmp_est <- t( apply(geno,1,function(X){
          model1 <- glm(Y ~ X, family = poisson(link = 'log'))
          coef <- model1$coefficients
          mu2 <- coef[1]
          beta <- coef[2]
          est <- try(unlist(PoissionGamma(Y, X, beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
          if(class(est) != "try-error"){ return(est) }
        })
        )

        ## Post process estimations
        tested.id <- match(rownames(tmp_est),geno.vcf@fix[,"ID"])
        ## calculate distance with respect to transcript(gene) strand
        distance <- if(as.character(strand(peak_bed.gr[i])) == "+"){
          as.integer(geno.vcf@fix[tested.id,"POS"])- round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2)
        }else{
          round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2) - as.integer(geno.vcf@fix[tested.id,"POS"])
          }
        report <- data.frame(
          SNP = paste(geno.vcf@fix[tested.id,"CHROM"],geno.vcf@fix[tested.id,"POS"],sep = ":"),
          SNPID = rownames(tmp_est),
          REF = geno.vcf@fix[tested.id,"REF"],
          ALT = geno.vcf@fix[tested.id,"ALT"],
          PEAK = paste0(Chromosome,":",peak_bed.gr[i]$thickStart,"-",peak_bed.gr[i]$thickEnd,"_",peak_bed.gr[i]$name,"_",strand( peak_bed.gr[i]) ),
          DISTANCE = distance,
          beta = tmp_est[,"beta"],
          pvalue = tmp_est[,"p_value"],
          psi = tmp_est[,"psi"]
        )

      }else{
        ## Check covariates
        if(!is.numeric(Covariates)){stop("Please convert covariates into numerical variables...")}
        Y <- unlist( round( phenoY[i,] ) )
        psi <- 10
        ## Test against each genotype
        tmp_est <- t( apply(geno,1,function(X1){
          X.all <- cbind(X1,Covariates) # new design matrix
          colnames(X.all) <- paste("X",1:ncol(X.all),sep = "")
          design.multiBeta <- formula( paste( "log(Y+1) ~ ",paste("X.all[,", 1:ncol(X.all),"]", sep = "",collapse = " + ")) )
          ## Run multi-beta PoissonGamma
          aa <- unlist(summary( lm( design.multiBeta ) )$coefficients[, 1])
          mu2 <- aa[1]
          beta <- aa[2:(ncol(X.all)+1 )]
          est <- try(unlist(PoissonGamma::PoissionGamma_multiple_beta(Y, X.all, beta, psi, mu2, gamma = 0.25, steps = 10, down = 0.1,psi_cutoff = maxPsi)))
          if(class(est) != "try-error"){ return(est) }
          })
        )

        ## Post process estimations
        tested.id <- match(rownames(tmp_est),geno.vcf@fix[,"ID"])
        ## calculate distance with respect to transcript(gene) strand
        distance <- if(as.character(strand(peak_bed.gr[i])) == "+"){
          as.integer(geno.vcf@fix[tested.id,"POS"])- round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2)
        }else{
          round((peak_bed.gr[i]$thickStart+peak_bed.gr[i]$thickEnd)/2) - as.integer(geno.vcf@fix[tested.id,"POS"])
        }
        report <- data.frame(
          SNP = paste(geno.vcf@fix[tested.id,"CHROM"],geno.vcf@fix[tested.id,"POS"],sep = ":"),
          SNPID = rownames(tmp_est),
          REF = geno.vcf@fix[tested.id,"REF"],
          ALT = geno.vcf@fix[tested.id,"ALT"],
          PEAK = paste0(Chromosome,":",peak_bed.gr[i]$thickStart,"-",peak_bed.gr[i]$thickEnd,"_",peak_bed.gr[i]$name,"_",strand( peak_bed.gr[i]) ),
          DISTANCE = distance,
          beta = tmp_est[,"beta1"],
          pvalue = tmp_est[,"p_value3"],
          psi = tmp_est[,"psi"]
        )

      }
      report
    }

  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  endTime <- Sys.time()
  cat(paste("Time used to test association: ",difftime(endTime, startTime, units = "mins")," mins... \n"))

  return(testResult)
}



.genoDosage <- function(x){
  return( stringr::str_count(x,"1") )
}

