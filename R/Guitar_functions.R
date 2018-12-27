.peakToGRangesList <- function(a, header=FALSE) {

  # message
  print("Converting BED12 to GRangesList")
  print("It may take a few minutes")

  # read bed file
  #a = read.table(filepath,sep="\t",header=header,stringsAsFactors =FALSE)
  # mcols_info =a[,13:length(a[1,])]
  a = a[,1:12]

  # get transcripts
  no_tx = length(a[,1])
  tx_id = 1:no_tx;
  tx_name = paste("line_",1:no_tx,sep="")
  tx_chrom = a[,1]
  tx_strand = a[,6]
  tx_start = a[,2]+1
  tx_end = a[,3]
  transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
  head(transcripts)



  # get genes
  tx_name = tx_name
  gene_id = as.character(a[,4])
  gene_id[is.na(gene_id)]="NA"
  gene=data.frame(tx_name,gene_id)

  #
  splicing <- lapply(1:no_tx, .spliceSingleTrans, a=a, tx_start=tx_start)
  splicing <- .combineListOfSplicing(splicing)

  # make txdb
  peaks = suppressWarnings(
    makeTxDb(transcripts=transcripts,
             splicings=splicing,
             genes=gene))

  # generate GRangesList
  tx <- exonsBy(peaks, "tx",use.names=TRUE)
  #mcols(tx) <- a

  return(tx)
}

.combineListOfSplicing <- function(t){

  a <- paste("t[[",1:length(t),"]]", sep="")
  a <- paste(a,collapse =",")
  a <- paste("rbind(",a,")",sep="")
  c <- parse(text=a)
  b <- suppressWarnings(eval(c))

  return(b)
}

.spliceSingleTrans <- function(i,a,tx_start) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[[10]])

  # get start
  temp = as.integer(strsplit(as.character(tx[[12]]), ",")[[1]]) + tx_start[i]
  exon_start=temp

  # get end
  temp = as.integer(strsplit(as.character(tx[[11]]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2

  # get CDS
  cds_start = exon_start
  cds_end = exon_end

  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  return(splicing_tx)
}

# make Guitar Coordinates from TranscriptDb object
#' @import rtracklayer
.makeGuitarCoordsFromTxDb <- function(txdb,
                                     maximalAmbiguity = 3,
                                     minimalComponentLength = 100,
                                     minimalNcRNALength = 300,
                                     noBins=100){

  parameter = list()
  parameter$txdb <- txdb
  parameter$maximalAmbiguity <- maximalAmbiguity # whether overlap with another transcript
  parameter$minimalComponentLength <- minimalComponentLength # minimal length required for each component
  parameter$minimalNcRNALength <- minimalNcRNALength
  parameter$noBins <- noBins

  # prepare the bins
  component <- .extractComponent(parameter)
  side <- .get2sides(component)

  # print
  print("Building Guitar Coordinates. It may take a few minutes ...")
  utr3_bin <- .makeGuitarCoordsFromGRangesList(component[["utr3"]],parameter$noBins)
  utr5_bin <- .makeGuitarCoordsFromGRangesList(component[["utr5"]],parameter$noBins)
  cds_bin <- .makeGuitarCoordsFromGRangesList(component[["cds"]],parameter$noBins)
  ncRNA_bin <- .makeGuitarCoordsFromGRangesList(component[["ncRNA"]],parameter$noBins)
  mRNA_front_bin <- .makeGuitarCoordsFromGRangesList(side[["mrna_front"]],parameter$noBins)
  mRNA_back_bin <- .makeGuitarCoordsFromGRangesList(side[["mrna_back"]],parameter$noBins)
  ncRNA_front_bin <- .makeGuitarCoordsFromGRangesList(side[["ncrna_front"]],parameter$noBins)
  ncRNA_back_bin <- .makeGuitarCoordsFromGRangesList(side[["ncrna_back"]],parameter$noBins)
  print("Guitar Coordinates Built ...")

  # group together
  mcols(utr3_bin) <- data.frame(mcols(utr3_bin),comp="UTR3",category="mRNA")
  mcols(utr5_bin) <- data.frame(mcols(utr5_bin),comp="UTR5",category="mRNA")
  mcols(cds_bin) <- data.frame(mcols(cds_bin),comp="CDS",category="mRNA")
  mcols(mRNA_front_bin) <- data.frame(mcols(mRNA_front_bin),comp="Front",category="mRNA")
  mcols(mRNA_back_bin) <- data.frame(mcols(mRNA_back_bin),comp="Back",category="mRNA")

  mcols(ncRNA_bin) <- data.frame(mcols(ncRNA_bin),comp="lncRNA",category="lncRNA")
  mcols(ncRNA_front_bin) <- data.frame(mcols(ncRNA_front_bin),comp="Front",category="lncRNA")
  mcols(ncRNA_back_bin) <- data.frame(mcols(ncRNA_back_bin),comp="Back",category="lncRNA")


  GuitarCoords <- suppressWarnings(c(mRNA_front_bin,utr5_bin, cds_bin, utr3_bin, mRNA_back_bin,
                                     ncRNA_front_bin,ncRNA_bin,ncRNA_back_bin))

  return(GuitarCoords)}

.extractComponent <- function(parameter){

  txdb <- parameter$txdb
  # ambiguity filter
  exons <- exonsBy(txdb, by = "tx",use.names=TRUE)
  noTx <- length(exons)
  print(paste("total",noTx,"transcripts extracted ..."));

  temp <- countOverlaps(exons, exons)
  ambiguityFilteredTx <- names(exons[temp < (parameter$maximalAmbiguity+2)])
  noTxLeft <- length(ambiguityFilteredTx)
  print(paste("total",noTxLeft,"transcripts left after ambiguity filter ..."))
  exons <- exons[ambiguityFilteredTx]


  # extract important components
  cds <- cdsBy(txdb, by = "tx",use.names=TRUE)
  utr5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  utr3 <- threeUTRsByTranscript(txdb, use.names=TRUE)

  # extract mRNAs
  flag_utr5 <- (sum(width(utr5)) > parameter$minimalComponentLength)
  name_utr5 <- names(utr5)[flag_utr5]
  flag_utr3 <- (sum(width(utr3)) > parameter$minimalComponentLength)
  name_utr3 <- names(utr3)[flag_utr3]
  flag_cds <- (sum(width(cds)) > parameter$minimalComponentLength)
  name_cds <- names(cds)[flag_cds]
  name_mRNA <- intersect(intersect(name_utr5,name_utr3),name_cds)
  name_filtered_mRNA <- intersect(name_mRNA,names(exons))
  cds_filtered <- cds[name_filtered_mRNA]
  utr5_filtered <- utr5[name_filtered_mRNA]
  utr3_filtered <- utr3[name_filtered_mRNA]
  print(paste("total",length(cds_filtered),"mRNAs left after component length filter ..."))

  # extract mRNAs
  all_mRNA <- unique(c(names(utr5),names(utr3),names(cds)))
  name_ncRNA <- setdiff(names(exons),all_mRNA)
  ncRNA <- exons[name_ncRNA]
  flag_ncRNA <-
    (sum(width(ncRNA)) > parameter$minimalComponentLength) &
    (sum(width(ncRNA)) > parameter$minimalNcRNALength)
  name_ncRNA <- names(ncRNA)[flag_ncRNA]
  ncRNA_filtered <- ncRNA[name_ncRNA]
  print(paste("total",length(ncRNA_filtered),"ncRNAs left after ncRNA length filter ..."))

  # return the result
  comp <- list(cds=cds_filtered,utr3=utr3_filtered,utr5=utr5_filtered,ncRNA=ncRNA_filtered)
  return(comp)}

.get2sides <- function(component) {
  utr3 <- component[["utr3"]]
  utr5 <- component[["utr5"]]
  ncrna <- component[["ncRNA"]]

  mrna_front <- .getNeighborhood(comp=utr5, side=5)
  mrna_back <- .getNeighborhood(comp=utr3, side=3)
  ncrna_front <- .getNeighborhood(comp=ncrna, side=5)
  ncrna_back <- .getNeighborhood(comp=ncrna, side=3)

  result <- list(mrna_front=mrna_front,
                 mrna_back=mrna_back,
                 ncrna_front=ncrna_front,
                 ncrna_back=ncrna_back)
  return(result)
}

.makeGuitarCoordsFromGRangesList <- function(comp,
                                            noBins=100,
                                            collapseGene=FALSE,
                                            width = 51) {

  # get all the check points
  tx_length <- as.numeric(sum(width(comp)))
  checkpoints_interval <- tx_length/noBins

  # get transcript name and strand
  tx_name <- names(comp)
  granges <- unlist(comp)
  tx <- granges[tx_name]
  strand <- as.character(as.data.frame(strand(tx))[[1]])
  chr <- as.character(as.data.frame(seqnames(tx))[[1]])

  # get coordinates
  t <- lapply(X=1:noBins,
              FUN=.makeGuitarCoordsForSingleIndex,
              comp=comp,
              noBins=noBins,
              checkpoints_interval=checkpoints_interval,
              strand=strand,
              chr=chr,
              tx_name=tx_name)

  GuitarCoords <- .combineListOfGRanges(t)

  if (collapseGene) {
    temp <- split(GuitarCoords, mcols(GuitarCoords)$pos, drop=FALSE)
    mcols(temp) <- data.frame(pos=unique(mcols(GuitarCoords)$pos))
    GuitarCoords <- temp
  }

  # return the result
  return(GuitarCoords)
}

.combineListOfGRanges <- function(t){
  txt <- "c(t[[1]]"
  for (i in 2:length(t)) {
    txt <- paste(txt,",t[[",i,"]]",sep="")
  }
  txt <- paste(txt,")",sep="")
  c <- parse(text=txt)

  # suppressWarnings
  # GuitarCoords <- eval(c)
  GuitarCoords <- suppressWarnings(eval(c))

  return(GuitarCoords)
}
.makeGuitarCoordsForSingleIndex <- function(
  index, comp, noBins, checkpoints_interval,
  strand, chr, tx_name,
  binWidth=51) {

  # get checkpoints
  checkpoints <- ceiling(checkpoints_interval*(index-0.5))
  checkpoints_transcript <- GRanges(seqnames=tx_name,
                                    IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                    strand=strand)

  # convert to genomic coordinates
  checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)

  # resize
  # binWidth <- 4*ceiling(checkpoints_interval)+1
  checkRegion_genomic <- resize(x=checkPoints_genomic,
                                width=binWidth,
                                fix="center")


  # add annotation information
  mcols(checkRegion_genomic) <- data.frame(txid=tx_name, pos=(index-0.5)/noBins, interval=checkpoints_interval)

  return(checkRegion_genomic)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
.multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

.getNeighborhood <- function(comp, side=5, Width=1000) {
  # transcript info
  tx_name <- names(comp)
  granges <- unlist(comp)
  tx <- granges[tx_name]
  strand <- as.character(as.data.frame(strand(tx))[[1]])
  chr <- as.character(as.data.frame(seqnames(tx))[[1]])

  if (side == 5) {
    checkpoints <- rep(1,length(tx))
    checkpoints_transcript <- GRanges(seqnames=tx_name,
                                      IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                      strand=strand)
    # convert to genomic coordinates
    checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)
    # resize
    checkRegion_genomic <- resize(x=checkPoints_genomic, width=Width, fix="end")
  } else if (side == 3) {
    checkpoints <- sum(width(comp))
    checkpoints_transcript <- GRanges(seqnames=tx_name,
                                      IRanges(start=checkpoints, end=checkpoints, names=tx_name),
                                      strand=strand)
    # convert to genomic coordinates
    checkPoints_genomic <- mapFromTranscripts(checkpoints_transcript, comp)
    # resize
    checkRegion_genomic <- resize(x=checkPoints_genomic, width=Width, fix="start")
  }

  # convert to list
  names(checkRegion_genomic) <- rep("",length(tx))
  sidelist <- split(checkRegion_genomic, tx_name, drop=TRUE)
  sidelist <- sidelist[tx_name]
  mapped_chr <- as.character(as.data.frame(seqnames(checkRegion_genomic))[[1]])
  mcols(sidelist) <- data.frame(mapped_chr)

  return(sidelist)
}

combinedGuitarPlot <- function(ct, comLength = c(0.136,0.459, 0.405)){
  # stack
  adjust=1

  # plot the figure
  # extract information of mRNA and lncRNA
  ct$weight <- ct$count # as numeric
  ct1 <- ct[ct$category=="mRNA",] # mRNA
  ct2 <- ct[ct$category=="lncRNA",] # lncRNA

  # disable notes
  pos=Feature=weight=NULL

  # remove DNA
  id1 <- which(match(ct1$comp,c("Front","Back")) >0 )
  ct1 <- ct1[-id1,]
  id2 <- which(match(ct2$comp,c("Front","Back")) >0 )
  ct2 <- ct2[-id2,]

  # normalize feature
  featureSet <- as.character(unique(ct$Feature))
  for (i in 1:length(featureSet)) {
    id <- (ct1$Feature==featureSet[i])
    ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])

    id <- (ct2$Feature==featureSet[i])
    ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
  }

  p2 <-
    ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) +
    ggtitle("Distribution on lncRNA")  +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    xlab("") +
    ylab("Frequency") +
    geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
    annotate("text", x = 0.5, y = -0.2, label = "lncRNA")+
    annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    theme(legend.position="bottom")

  weight <- comLength/sum(comLength)
  names(weight) <- c("5'UTR","CDS","3'UTR")

  # density
  cds_id <- which(ct1$comp=="CDS")
  utr3_id <- which(ct1$comp=="UTR3")
  utr5_id <- which(ct1$comp=="UTR5")
  ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
  ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
  ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]

  # re-normalization
  featureSet <- as.character(unique(ct$Feature))
  for (i in 1:length(featureSet)) {
    id <- (ct1$Feature==featureSet[i])
    ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
  }

  # stratch
  x <- cumsum(weight)
  ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
  ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
  ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]

  p1 <-
    ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
    ggtitle("Distribution on mRNA") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
    xlab("") +
    ylab("Frequency") +
    geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
    annotate("text", x = x[1]/2, y = -0.2, label = "5'UTR") +
    annotate("text", x = x[1] + weight[2]/2, y = -0.2, label = "CDS") +
    annotate("text", x = x[2] + weight[3]/2, y = -0.2, label = "3'UTR") +
    theme(legend.position="bottom") +
    geom_vline(xintercept= x[1:2], linetype="dotted") +
    annotate("rect", xmin = 0, xmax = x[1], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    annotate("rect", xmin = x[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
    annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")

  .multiplot(p1, p2, cols=2)

}


#' @title GuitarPlotNew
#' @export
GuitarPlotNew <- function(gfeatures,
                          GuitarCoordsFromTxDb=NA,
                          txdb=NA,
                          genome=NA,
                          noBins=10,
                          saveToPDFprefix=NA,
                          returnCount=FALSE,
                          includeNeighborDNA=FALSE,
                          maximalFeatureAmbiguity=5,
                          rescaleComponent=TRUE,
                          fill=FALSE,
                          adjust=1){
  
  # make sure the Guitar coordinates are available
  suppressWarnings(
    if (is.na(GuitarCoordsFromTxDb)&is.na(txdb)&is.na(genome)) {
      stop("Must provide at least one of the three: GuitarCoords, txdb or genome")
    }
  )
  
  if ( suppressWarnings(is.na(GuitarCoordsFromTxDb)) ) {
    if (suppressWarnings(is.na(txdb))) {
      print("Downloading Transcriptome Information from UCSC ...")
      txdb <- suppressMessages(makeTxDbFromUCSC(genome=genome))
      print("Making Guitar Coordinates ...")
      GuitarCoordsFromTxDb <- suppressMessages(makeGuitarCoordsFromTxDb(txdb))
      GuitarCoords <- GuitarCoordsFromTxDb
    } else {
      print("Making Guitar Coordinates from provided TranscriptDb Object ...")
      GuitarCoordsFromTxDb <- makeGuitarCoordsFromTxDb(txdb, noBins=noBins)
      GuitarCoords <- GuitarCoordsFromTxDb
    }
  } else {
    print("Using provided Guitar Coordinates")
    GuitarCoords <- GuitarCoordsFromTxDb
  }
  
  # Generate named List
  noGroup <- length(gfeatures)
  group_names <- names(gfeatures)
  m <- gfeatures
  if (is.null(group_names)) {
    group_names <- paste("item",1:noGroup)
  }
  
  print("resolving ambiguious features ...")
  for (i in 1:noGroup) {
    temp = .countGuitarDensity(
      gfeatures[[i]],
      GuitarCoords,
      maximalFeatureAmbiguity)
    temp = cbind(temp,Feature=group_names[i])
    m[[i]] =temp
  }
  ct=.combineListOfDataFrame(m)
  ct[[4]] <- as.character(ct[[4]])
  ct <- ct[ct$count>0,]  # remove redundant info
  
  # organize output
  if (fill==FALSE) {
    .makeFigure_nofill(ct,
                       GuitarCoordsFromTxDb,
                       includeNeighborDNA,
                       rescaleComponent,
                       saveToPDFprefix,
                       adjust)
  } else {
    .makeFigure_fill(ct,
                     GuitarCoordsFromTxDb,
                     includeNeighborDNA,
                     rescaleComponent,
                     saveToPDFprefix,
                     adjust)
  }
  
  
  # return the result
  if (returnCount) {return(ct)}
}

#' @import GenomicRanges
.countGuitarDensity <- function(peak,GuitarCoords,maximalFeatureAmbiguity) {
  
  # count overlaps
  n <- countOverlaps(GuitarCoords,peak)
  
  # normalize by overlaps
  tx <- as.character(mcols(GuitarCoords)[[1]])
  txlist <- split(x=GuitarCoords,f=tx)
  peak_overlap <- countOverlaps(peak,txlist)
  
  # split features based on ambiguity
  z <- matrix(rep(n,maximalFeatureAmbiguity),ncol=maximalFeatureAmbiguity)
  for (i in 1:maximalFeatureAmbiguity) {
    temp <- countOverlaps(GuitarCoords,peak[peak_overlap==i])
    z[,i] = temp/i
  }
  n <- rowSums(z)
  q <- data.frame(mcols(GuitarCoords),count=as.numeric(n))
  return(q)
  
}
.combineListOfDataFrame <- function(t){
  if (length(t)==1) {
    return(t[[1]])
  } else {
    txt <- "rbind(t[[1]]"
    for (i in 2:length(t)) {
      txt <- paste(txt,",t[[",i,"]]",sep="")
    }
    txt <- paste(txt,")",sep="")
    c <- parse(text=txt)
    newframe <- eval(c)
    return(newframe)
  }
}





.makeFigure_nofill <- function(ct,
                               GuitarCoordsFromTxDb,
                               includeNeighborDNA,
                               rescaleComponent,
                               saveToPDFprefix,adjust=adjust) {
  
  # extract information of mRNA and lncRNA
  ct$weight <- ct$count # as numeric
  ct1 <- ct[ct$category=="mRNA",] # mRNA
  ct2 <- ct[ct$category=="lncRNA",] # lncRNA
  
  # save(ct1,ct2, includeNeighborDNA,rescaleComponent,saveToPDFprefix,GuitarCoordsFromTxDb, file = "data.Rdata")
  d <- mcols(GuitarCoordsFromTxDb)
  
  # disable notes
  pos=Feature=weight=NULL
  
  # plot
  if (includeNeighborDNA) {
    if (rescaleComponent==FALSE) {
      
      # normalize feature
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
        
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
      }
      
      # adjust position mRNA
      pos_adjust <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))-1
      ct1$pos <- ct1$pos + pos_adjust
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=5*weight))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature), colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = 1.5, y = -0.2, label = "5'UTR") +
        annotate("text", x = 2.5, y = -0.2, label = "CDS") +
        annotate("text", x = 0.5, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = 4.5, y = -0.2, label = "Tail (1kb)") +
        annotate("text", x = 3.5, y = -0.2, label = "3'UTR")  +
        geom_vline(xintercept=1:4, linetype="dotted") +
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 3, xmax = 4, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black") +
        xlim(0,5) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
      
      # adjust position lncRNA
      pos_adjust <- match(ct2$comp,c("Front","lncRNA","Back"))-1
      ct2$pos <- ct2$pos + pos_adjust
      p2 <-
        ggplot(ct2, aes(x=pos, group=Feature, weight=3*weight)) +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature), colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = 1.5, y = -0.2, label = "lncRNA") +
        annotate("text", x = 0.5, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = 2.5, y = -0.2, label = "Tail (1kb)")  +
        geom_vline(xintercept=1:2, linetype="dotted") +
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        xlim(0,3) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    } else {
      
      if (rescaleComponent) {
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"mRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=5)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","5'UTR","CDS","3'UTR","Tail")
        w1 <- weight
        
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"lncRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=3)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","lncRNA","Tail")
        w2 <- weight
        
        x1 <- cumsum(w1)
        x2 <- cumsum(w2)
      }
      
      
      # adjust position
      id <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))
      ct1$count <- ct1$count*w1[id]
      ct1$pos <- ct1$pos*w1[id] + c(0,x1)[id]
      
      id <- match(ct2$comp,c("Front","lncRNA","Back"))
      ct2$count <- ct2$count*w2[id]
      ct2$pos <- ct2$pos*w2[id] + c(0,x2)[id]
      
      
      # normalize
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$count[id]/sum(ct2$count[id])
      }
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = sum(x1[1:2])/2, y = -0.2, label = "5'UTR") +
        annotate("text", x = sum(x1[2:3])/2, y = -0.2, label = "CDS") +
        annotate("text", x = x1[1]/2, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = sum(x1[4:5])/2, y = -0.2, label = "Tail (1kb)") +
        annotate("text", x = sum(x1[3:4])/2, y = -0.2, label = "3'UTR")  +
        geom_vline(xintercept=x1[1:4], linetype="dotted") +
        annotate("rect", xmin = x1[1], xmax = x1[2], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[3], xmax = x1[4], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[2], xmax = x1[3], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black") +
        xlim(0,1) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
      
      
      p2 <-
        ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = sum(x2[1:2])/2, y = -0.2, label = "lncRNA") +
        annotate("text", x = x2[1]/2, y = -0.2, label = 'Promoter (1kb)') +
        annotate("text", x = sum(x2[2:3])/2, y = -0.2, label = "Tail (1kb)")  +
        geom_vline(xintercept=x2[1:2], linetype="dotted") +
        annotate("rect", xmin = x2[1], xmax = x2[2], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        xlim(0,1) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    }
    
  }
  
  if (includeNeighborDNA==FALSE) {
    
    # remove DNA
    id1 <- which(match(ct1$comp,c("Front","Back")) >0 )
    ct1 <- ct1[-id1,]
    id2 <- which(match(ct2$comp,c("Front","Back")) >0 )
    ct2 <- ct2[-id2,]
    
    # normalize feature
    featureSet <- as.character(unique(ct$Feature))
    for (i in 1:length(featureSet)) {
      id <- (ct1$Feature==featureSet[i])
      ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
      
      id <- (ct2$Feature==featureSet[i])
      ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
    }
    
    p2 <-
      ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) +
      ggtitle("Distribution on lncRNA")  +
      xlab("") +
      ylab("Frequency") +
      geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
      annotate("text", x = 0.5, y = -0.2, label = "lncRNA")+
      annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
      theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                         legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    
    if (rescaleComponent) {
      
      # normalization by length of components in mRNA
      # calculate relative length of each components
      temp <- unique(d[,c(1,3,4,5)])
      id1 <- which(match(temp$comp,c("Front","Back")) >0 )
      temp <- temp[-id1,] # remove DNA
      id1 <- which(match(temp$category,"mRNA") >0 )
      temp <- temp[id1,]
      temp <- matrix(temp$interval,ncol=3)
      temp <- temp/rowSums(temp)
      temp <- colSums(temp)
      temp <-temp/sum(temp)
      weight <- temp
      names(weight) <- c("5'UTR","CDS","3'UTR")
      
      
      # density
      cds_id <- which(ct1$comp=="CDS")
      utr3_id <- which(ct1$comp=="UTR3")
      utr5_id <- which(ct1$comp=="UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]
      
      # re-normalization
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      
      # stratch
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
      ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
        ggtitle("Distribution on mRNA") +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2) +
        annotate("text", x = x[1]/2, y = -0.2, label = "5'UTR") +
        annotate("text", x = x[1] + weight[2]/2, y = -0.2, label = "CDS") +
        annotate("text", x = x[2] + weight[3]/2, y = -0.2, label = "3'UTR") +
        geom_vline(xintercept= x[1:2], linetype="dotted") +
        annotate("rect", xmin = 0, xmax = x[1], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")+
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
      
      
      
    } else {
      
      pos_adjust <- match(ct1$comp,c("UTR5","CDS","UTR3"))-1
      ct1$pos <- ct1$pos + pos_adjust
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature,colour=factor(Feature), weight=3*weight))  +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature)),alpha=0.2) +
        annotate("text", x = 0.5, y = -0.2, label = "5'UTR") +
        annotate("text", x = 1.5, y = -0.2, label = "CDS") +
        annotate("text", x = 2.5, y = -0.2, label = "3'UTR") +
        geom_vline(xintercept=1:2, linetype="dotted") +
        annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")+
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
    }
  }
  
  suppressWarnings(
    if (is.na(saveToPDFprefix)) {
      # return the result
      print("no figure saved ...")
      .multiplot(p1, p2, cols=2)
    }  else {
      f1 <- paste(saveToPDFprefix,"_Guitar.pdf",sep="")
      
      pdf(file=f1,width=8, height=4)
      .multiplot(p1, p2, cols=2)
      dev.off()
      print(paste("Figures saved into",f1,"...", sep=" "))
    }
  )
  
  #suppressWarnings( .multiplot(p1, p2, cols=2))
  
}


.makeFigure_fill <- function(ct,
                             GuitarCoordsFromTxDb,
                             includeNeighborDNA,
                             rescaleComponent,
                             saveToPDFprefix,adjust=adjust) {
  
  # extract information of mRNA and lncRNA
  ct$weight <- ct$count # as numeric
  ct1 <- ct[ct$category=="mRNA",] # mRNA
  ct2 <- ct[ct$category=="lncRNA",] # lncRNA
  
  # save(ct1,ct2, includeNeighborDNA,rescaleComponent,saveToPDFprefix,GuitarCoordsFromTxDb, file = "data.Rdata")
  
  d <- mcols(GuitarCoordsFromTxDb)
  
  # disable notes
  pos=Feature=weight=NULL
  
  # plot
  if (includeNeighborDNA) {
    if (rescaleComponent==FALSE) {
      
      # normalize feature
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
        
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
      }
      
      # adjust position mRNA
      pos_adjust <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))-1
      ct1$pos <- ct1$pos + pos_adjust
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=5*weight))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature), colour=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = 1.5, y = -0.2, label = "5'UTR") +
        annotate("text", x = 2.5, y = -0.2, label = "CDS") +
        annotate("text", x = 0.5, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = 4.5, y = -0.2, label = "Tail (1kb)") +
        annotate("text", x = 3.5, y = -0.2, label = "3'UTR")  +
        geom_vline(xintercept=1:4, linetype="dotted") +
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 3, xmax = 4, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black") +
        xlim(0,5) +theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                      legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
      
      
      # adjust position lncRNA
      pos_adjust <- match(ct2$comp,c("Front","lncRNA","Back"))-1
      ct2$pos <- ct2$pos + pos_adjust
      p2 <-
        ggplot(ct2, aes(x=pos, group=Feature, weight=3*weight)) +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature), colour=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = 1.5, y = -0.2, label = "lncRNA") +
        annotate("text", x = 0.5, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = 2.5, y = -0.2, label = "Tail (1kb)")  +
        geom_vline(xintercept=1:2, linetype="dotted") +
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        xlim(0,3) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    } else {
      
      if (rescaleComponent) {
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"mRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=5)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","5'UTR","CDS","3'UTR","Tail")
        w1 <- weight
        
        temp <- unique(d[,c(1,3,4,5)])
        id1 <- which(match(temp$category,"lncRNA") >0 )
        temp <- temp[id1,]
        temp <- matrix(temp$interval,ncol=3)
        temp <- temp/rowSums(temp)
        temp <- colSums(temp)
        temp <- temp/sum(temp)
        weight <- temp
        names(weight) <- c("Promoter","lncRNA","Tail")
        w2 <- weight
        
        x1 <- cumsum(w1)
        x2 <- cumsum(w2)
      }
      
      
      # adjust position
      id <- match(ct1$comp,c("Front","UTR5","CDS","UTR3","Back"))
      ct1$count <- ct1$count*w1[id]
      ct1$pos <- ct1$pos*w1[id] + c(0,x1)[id]
      
      id <- match(ct2$comp,c("Front","lncRNA","Back"))
      ct2$count <- ct2$count*w2[id]
      ct2$pos <- ct2$pos*w2[id] + c(0,x2)[id]
      
      
      # normalize
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
        id <- (ct2$Feature==featureSet[i])
        ct2$weight[id] <- ct2$count[id]/sum(ct2$count[id])
      }
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
        scale_x_continuous(minor_breaks = seq(0, 5, 1)) +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = sum(x1[1:2])/2, y = -0.2, label = "5'UTR") +
        annotate("text", x = sum(x1[2:3])/2, y = -0.2, label = "CDS") +
        annotate("text", x = x1[1]/2, y = -0.2, label = "Promoter (1kb)") +
        annotate("text", x = sum(x1[4:5])/2, y = -0.2, label = "Tail (1kb)") +
        annotate("text", x = sum(x1[3:4])/2, y = -0.2, label = "3'UTR")  +
        geom_vline(xintercept=x1[1:4], linetype="dotted") +
        annotate("rect", xmin = x1[1], xmax = x1[2], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[3], xmax = x1[4], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x1[2], xmax = x1[3], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black") +
        xlim(0,1) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5))
      
      
      p2 <-
        ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) +
        ggtitle("Distribution on lncRNA")  +
        scale_x_continuous(minor_breaks = seq(6, 9, 1)) +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = sum(x2[1:2])/2, y = -0.2, label = "lncRNA") +
        annotate("text", x = x2[1]/2, y = -0.2, label = 'Promoter (1kb)') +
        annotate("text", x = sum(x2[2:3])/2, y = -0.2, label = "Tail (1kb)")  +
        geom_vline(xintercept=x2[1:2], linetype="dotted") +
        annotate("rect", xmin = x2[1], xmax = x2[2], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        xlim(0,1) +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    }
    
  }
  
  if (includeNeighborDNA==FALSE) {
    
    # remove DNA
    id1 <- which(match(ct1$comp,c("Front","Back")) >0 )
    ct1 <- ct1[-id1,]
    id2 <- which(match(ct2$comp,c("Front","Back")) >0 )
    ct2 <- ct2[-id2,]
    
    # normalize feature
    featureSet <- as.character(unique(ct$Feature))
    for (i in 1:length(featureSet)) {
      id <- (ct1$Feature==featureSet[i])
      ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
      
      id <- (ct2$Feature==featureSet[i])
      ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
    }
    
    p2 <-
      ggplot(ct2, aes(x=pos, group=Feature, weight=weight)) +
      ggtitle("Distribution on lncRNA")  +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
      xlab("") +
      ylab("Frequency") +
      geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2,position="fill") +
      annotate("text", x = 0.5, y = -0.2, label = "lncRNA")+
      annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
      theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                         legend.text = element_text(face = "bold"), legend.title = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5))
    
    if (rescaleComponent) {
      
      # normalization by length of components in mRNA
      # calculate relative length of each components
      temp <- unique(d[,c(1,3,4,5)])
      id1 <- which(match(temp$comp,c("Front","Back")) >0 )
      temp <- temp[-id1,] # remove DNA
      id1 <- which(match(temp$category,"mRNA") >0 )
      temp <- temp[id1,]
      temp <- matrix(temp$interval,ncol=3)
      temp <- temp/rowSums(temp)
      temp <- colSums(temp)
      temp <-temp/sum(temp)
      weight <- temp
      names(weight) <- c("5'UTR","CDS","3'UTR")
      
      
      # density
      cds_id <- which(ct1$comp=="CDS")
      utr3_id <- which(ct1$comp=="UTR3")
      utr5_id <- which(ct1$comp=="UTR5")
      ct1$count[utr5_id] <- ct1$count[utr5_id]*weight["5'UTR"]
      ct1$count[cds_id] <- ct1$count[cds_id]*weight["CDS"]
      ct1$count[utr3_id] <- ct1$count[utr3_id]*weight["3'UTR"]
      
      # re-normalization
      featureSet <- as.character(unique(ct$Feature))
      for (i in 1:length(featureSet)) {
        id <- (ct1$Feature==featureSet[i])
        ct1$weight[id] <- ct1$count[id]/sum(ct1$count[id])
      }
      
      # stratch
      x <- cumsum(weight)
      ct1$pos[utr5_id] <- ct1$pos[utr5_id]*weight["5'UTR"] + 0
      ct1$pos[cds_id] <- ct1$pos[cds_id]*weight["CDS"] + x[1]
      ct1$pos[utr3_id] <- ct1$pos[utr3_id]*weight["3'UTR"] + x[2]
      
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature, weight=weight))  +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature),colour=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = x[1]/2, y = -0.2, label = "5'UTR") +
        annotate("text", x = x[1] + weight[2]/2, y = -0.2, label = "CDS") +
        annotate("text", x = x[2] + weight[3]/2, y = -0.2, label = "3'UTR") +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5)) +
        geom_vline(xintercept= x[1:2], linetype="dotted") +
        annotate("rect", xmin = 0, xmax = x[1], ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[2], xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = x[1], xmax = x[2], ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")
      
      
    } else {
      
      pos_adjust <- match(ct1$comp,c("UTR5","CDS","UTR3"))-1
      ct1$pos <- ct1$pos + pos_adjust
      p1 <-
        ggplot(ct1, aes(x=pos, group=Feature,colour=factor(Feature), weight=3*weight))  +
        ggtitle("Distribution on mRNA") +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        xlab("") +
        ylab("Frequency") +
        geom_density(adjust=adjust,aes(fill=factor(Feature)),alpha=0.2,position="fill") +
        annotate("text", x = 0.5, y = -0.2, label = "5'UTR") +
        annotate("text", x = 1.5, y = -0.2, label = "CDS") +
        annotate("text", x = 2.5, y = -0.2, label = "3'UTR") +
        theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           legend.position = "none", plot.title = element_text(face = "bold",hjust = 0.5)) +
        annotate("rect", xmin = 0, xmax = 1, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 2, xmax = 3, ymin = -0.12, ymax = -0.08, alpha = .99, colour = "black")+
        annotate("rect", xmin = 1, xmax = 2, ymin = -0.16, ymax = -0.04, alpha = .2, colour = "black")
    }
  }
  
  suppressWarnings(
    if (is.na(saveToPDFprefix)) {
      # return the result
      print("no figure saved ...")
      .multiplot(p1, p2, cols=2)
    }  else {
      f1 <- paste(saveToPDFprefix,"_Guitar.pdf",sep="")
      
      pdf(file=f1,width=8, height=4)
      .multiplot(p1, p2, cols=2)
      dev.off()
      print(paste("Figures saved into",f1,"...", sep=" "))
    }
  )
  
  #suppressWarnings( .multiplot(p1, p2, cols=2))
  
}
