#' @import DESeq2
#' @export
setGeneric("counts", getGeneric("counts", package = "DESeq2"))

#' @export
#' @rdname counts
setMethod("counts", signature("MeRIP"), function(object){
  object@reads
})

#' @export
#' @rdname Input.files
setGeneric("Input.files", function(object) {
  standardGeneric("Input.files")
})

#' @export
#' @rdname IP.files
setGeneric("IP.files", function(object) {
  standardGeneric("IP.files")
})

#' @export
setGeneric("geneBins", function(object) {
  standardGeneric("geneBins")
})

#' @export
setGeneric("jointPeak", function(object) {
  standardGeneric("jointPeak")
})

#' @export
setGeneric("filter", function(object, ... ) {
  standardGeneric("filter")
})

#' @export
setGeneric("extractInput", function(object) {
  standardGeneric("extractInput")
})

#' @export
setGeneric("extractIP", function(object, ...) {
  standardGeneric("extractIP")
})

#' @export
setGeneric("PrepCoveragePlot", function(object, ...) {
  standardGeneric("PrepCoveragePlot")
})

#' @export
setGeneric("normalizeLibrary",function(object, ...){
  standardGeneric("normalizeLibrary")
})

#' @export
setGeneric("adjustExprLevel",function(object, ...){
  standardGeneric("adjustExprLevel")
})

#' @export
setGeneric("geneExpression",function(object, ...){
  standardGeneric("geneExpression")
})

#' @export
setGeneric("consistentPeak",function(object, ...){
  standardGeneric("consistentPeak")
})

#' @export
setGeneric("variable",function(object){
  standardGeneric("variable")
})

#' @export
setGeneric("variable<-",function(object, value){
  standardGeneric("variable<-")
})

#' @export
setGeneric("samplenames",function(object){
  standardGeneric("samplenames")
})

#' @export
setGeneric("samplenames<-",function(object, value){
  standardGeneric("samplenames<-")
})

#' @export
setGeneric("QNBtest",function(object){
  standardGeneric("QNBtest")
})

#' @export
setGeneric("peakDistribution",function(object){
  standardGeneric("peakDistribution")
})

#' @export
setGeneric("plotGeneCov", function(object, geneName, libraryType, center,ZoomIn, adjustExprLevel ){ standardGeneric("plotGeneCov")})

#' @export
setGeneric("plotSNPpeakPairs",function(object, genotypeFile, SNPID, geneName, libraryType , center ,ZoomIn, adjustExprLevel ){standardGeneric("plotSNPpeakPairs")})

#' @export
setGeneric("geneExpressionTMP",function(object, meanFragmentLength = 150, normalize = T){
  standardGeneric("geneExpressionTMP")
})

#' @export
setGeneric("RADARtest",function(object,  exclude ,maxPsi){
  standardGeneric("RADARtest")
})

#' @export
setGeneric("select" ,function(object, samples){
  standardGeneric("select")
})

#' @export
setGeneric("results", function(object){standardGeneric("results")})

#' @export
setGeneric("betaBin",function(object){standardGeneric("betaBin")})





 