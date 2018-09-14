## the main function to generate fixed size random sampled windows on transcriptome given annotation GTF file
randomPeaks = function(size,gtf){
  library(GenomicFeatures)

  geneModel = gtfToGeneModel(gtf)
  no.genes = length(geneModel)
  bed12=data.frame() # initiate the bed12 data frame to store the peaks
  pb <-  txtProgressBar(min = 1, max = no.genes, style = 3)
  for(i in 1:no.genes){
    # DNA location to gene location conversion
    df.geneModel= as.data.frame(reduce(geneModel[i][[1]]) )##data frame of gene model
    dna.range = as.data.frame(range(geneModel[i]))
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    DNA2RNA=cumsum(DNA2RNA)*DNA2RNA

    #creat a corresponding map from RNA to DNA
    RNA2DNA = 1:exon.length
    pointer = 1
    for (j in 1:no.exon){
      RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
      pointer = pointer + df.geneModel$width[j]
    }
    RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates

    no.peak.to.sample = round(exon.length/1000)
    peaks.rna = ceiling(runif(no.peak.to.sample,size/2+1,exon.length-size/2-1))

    if(no.peak.to.sample > 0 ){
      peaks.dna = as.data.frame ( t(sapply(peaks.rna,function(x,RNA2DNA,size,strand,chrom){
        left = RNA2DNA[x-size/2]
        right = RNA2DNA[x + size/2]
        return(c(chrom,left,right,strand) )
      },RNA2DNA = RNA2DNA,size = size,strand = as.character(dna.range$strand), chrom =as.character(dna.range$seqnames))  ) )
      colnames(peaks.dna)=c("chr","start","end","strand")
      peak.gr = makeGRangesFromDataFrame( peaks.dna )
      for(j in 1:no.peak.to.sample){
        tmp = GenomicRanges::intersect(peak.gr[j],reduce(geneModel[i][[1]]) )
        bed.tmp = data.frame(matrix(nrow=1,ncol=12))
        colnames(bed.tmp)=c("chr","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSize","blockStart")
        bed.tmp["chr"] = unique(as.character(seqnames(tmp)))
        bed.tmp[c(2,3)] = as.data.frame(range(tmp))[1,c(2,3)]
        bed.tmp["strand"]= unique(as.character(strand(tmp)))
        bed.tmp["name"] = names(geneModel)[i]
        bed.tmp["score"] = 1
        bed.tmp[c("thickStart","thickEnd")] = bed.tmp[c(2,3)]
        bed.tmp["itemRgb"] = NA
        bed.tmp["blockCount"] = length(tmp)
        bed.tmp["blockSize"] = paste(as.data.frame(tmp)[,4],collapse=",")
        bed.tmp["blockStart"] = paste(as.data.frame(tmp)[,2] - rep(as.numeric(bed.tmp[2]),as.numeric(bed.tmp["blockCount"])),collapse=",")
        bed12 = rbind(bed12,bed.tmp)
      }

    }

    setTxtProgressBar(pb, i)
  }
  return(bed12)
}

