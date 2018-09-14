#' @title qqplot.pvalue
#' @param x can be a vector (p value of one group) or a list of vector (p value of multiple groups).
#' @param pointSize The size of data points
#' @param legendSize The size of points in the legend
#' @export
qqplot.pvalue <- function(x,pointSize = 1,legendSize = 4){
  library(ggplot2)
  if(is.list(x)){
    nn<-sapply(x, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(x))) {
      grp=factor(rep(names(x), nn), levels=names(x))
      names(x)<-NULL
    } else {
      grp=factor(rep(1:length(x), nn))
    }
    pvo<-x
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
      exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
    }
    thin <- unique(data.frame(pvalues = round(pvalues, 3),
                              exp.x = round(exp.x, 3),
                              grp=grp))
    grp = thin$grp
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
    qq.melt <- data.frame(expected=exp.x, Observed = pvalues, label = grp)
    ggplot(data = qq.melt, aes(expected, Observed,colour = label) )+geom_point(size = pointSize)+geom_abline(slope = 1)+xlab("Expected(-log10 p-value)")+ylab("Observed(-log10 p-value)")+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      guides(color = guide_legend(override.aes = list(size=legendSize)))

  }else{
    expected <- -log10(seq(0,1,1/(length(x)-1 ) ))
    p <- -log10(sort(x))

    ## remove some
    tail <- which(p>quantile(x,0.75) )
    head <-seq(1,length(which(p<=quantile(x,0.75) )),10 )
    qq.data <- data.frame(expected=c(expected[head],expected[tail]), Observed=c(p[head],p[tail]) )
    ggplot(data = qq.data, aes(expected, Observed) )+geom_point(size = 1.5,colour = "blue")+geom_abline(slope = 1)+xlab("Expected(-log10 p-value)")+ylab("Observed(-log10 p-value)")+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  }
}
