
#' @export
MeRIP <- setClass("MeRIP",
         representation( reads = "matrix",
                         binSize = "numeric",
                         geneModel = "GRangesList",
                         gtf = "character",
                         bamPath.input = "character",
                         bamPath.ip = "character",
                         samplenames = "character",
                         geneBins = "data.frame",
                         geneSum = "matrix",
                         GTF = "GRanges"
                         ),
         validity = function(object){
           errors <- character()
           ## Check basic element
           if( all(dim(object@reads) <= 1) ){
             errors <- c(errors, paste0("read count table is a required data for MeRIP instance!"))
           }
           if( length(object@binSize)==0){
             errors <- c(errors, paste0("bin size to slice the transcript is a required parameter for MeRIP instance!"))
           }
           if( length(object@geneModel)<=0 ){
             errors <- c(errors, paste0("GeneModel is a required data for MeRIP instance!"))
           }
           if( length(object@bamPath.input)<=0){
             errors <- c(errors, paste0("path to the input BAM files are the required data for MeRIP instance!"))
           }
           if( length(object@bamPath.ip)<=0){
             errors <- c(errors, paste0("path to the IP BAM files are the required data for MeRIP instance!"))
           }
           if( length(object@samplenames)<=0){
             errors <- c(errors, paste0("Sample names are the required data for MeRIP instance!"))
           }
           ## match data dimension with num of samples
           if(ncol(object@reads) !=  2 * length(object@samplenames) ){
             errors <- c(errors, paste0("The number of colnumns of read count is ",ncol(object@reads),". The number of samples is ",length(object@samplenames),". The number of colnumns of read count should be 2x the number of samples!"))
           }
           if(length(object@bamPath.input) != length(object@samplenames) ){
             errors <- c(errors, paste0("The number of input bam file path(es) should equal to the number of samplenames!"))
           }
           if(length(object@bamPath.ip) != length(object@samplenames) ){
             errors <- c(errors, paste0("The number of IP bam file path(es) should equal to the number of samplenames!"))
           }
           
           if (length(errors) == 0) TRUE else errors
         },
         prototype(geneBins = data.frame(gene = character(), bin = character() ), geneSum = matrix(),GTF = GRanges())
         )

#' @export
MeRIP.Peak <- setClass("MeRIP.Peak",representation( peakCallResult = "matrix",
                                                    jointPeak_id_pairs = "matrix",
                                                    jointPeaks = "data.frame",
                                                    jointPeak_ip = "matrix",
                                                    jointPeak_input = "matrix",
                                                    norm.jointPeak_ip = "matrix",
                                                    sizeFactor = "data.frame",
                                                    variate = "data.frame",
                                                    jointPeak_adjExpr = "matrix",
                                                    test.est = "matrix",
                                                    peakCalling = "character",
                                                    jointPeak_threshold = "numeric",
                                                    test.method = "character"),
                       contains = "MeRIP",
                       prototype(peakCalling = "none", jointPeak_threshold = 0, test.method = "none")
                       )



#' @export
setMethod("show",signature("MeRIP"), function(object){
 summary(object)
})
#' @export
setMethod("show",signature("MeRIP.Peak"), function(object){
  summary(object)
})


#' @export
setMethod("summary", signature("MeRIP"), function(object){
  cat(paste0("MeRIP dataset of ",length(object@samplenames)," samples.\n"))
  cat("The total read count for Input and IP samples are (Million reads):\n")
  totoReads <- rbind("Input" = round(colSums(object@reads[,1:length(object@samplenames)])/1e6, digits = 2),
                     "IP" = round(colSums(object@reads[,-c(1:length(object@samplenames))])/1e6, digits = 2))
  colnames(totoReads) <- object@samplenames
  print(totoReads)
})
#' @export
setMethod("summary", signature("MeRIP.Peak"), function(object){
  cat(paste0("MeRIP dataset of ",length(object@samplenames)," samples.\n"))
  cat("The total read count for Input and IP samples are (Million reads):\n")
  totoReads <- rbind("Input" = round(colSums(object@reads[,1:length(object@samplenames)])/1e6, digits = 2),
                     "IP" = round(colSums(object@reads[,-c(1:length(object@samplenames))])/1e6, digits = 2))
  colnames(totoReads) <- object@samplenames
  
  print(totoReads)
  if(object@peakCalling != "none"){cat(paste0("\nPeak calling done by ",object@peakCalling,".\n"))}
  if(nrow(object@jointPeaks)>0 & object@jointPeak_threshold> 0 ){
    cat(paste0(nrow(object@jointPeaks)," joint peak reported at threshold ",object@jointPeak_threshold," (requiring peaks to be called in at least ",object@jointPeak_threshold," samples).\n") )
  }
  if(nrow(object@geneSum)>0){ cat("Input gene level count available.\n") }
  if(nrow(object@variate)>0){ cat(paste0("There are ",ncol(object@variate)," predictor variables/covariates. Can access by function variable(MeRIPdata). \n"))}
  if(object@test.method != "none" & nrow(object@test.est)>0){ print(paste0("Differential peaks tested by ",object@test.method,".\n"))}
})


## A function to convert old version list data into S4 class MeRIP dataset
#' @export
makeMeRIPfromList <- function(x, gtf){
  if( all(c("reads","binSize","geneModel","bamPath.input","bamPath.ip","samplenames") %in% names(x) ) ){
    return(MeRIP(reads = x$reads, binSize = x$binSize, geneModel = x$geneModel, gtf = gtf ,bamPath.input = x$bamPath.input, bamPath.ip = x$bamPath.ip, samplenames = x$samplenames) )
  }else{
    stop("The list must have the most basic elements for MeRIP dataset! ")
      }

}



