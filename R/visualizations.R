#' t-SNE visualization of the rare cells
#'
#'Given the fractions of desired rare cells, this function visualize the determined rare cells by t-SNE compared to all cells in data set
#'
#' @param tsne the tsne object returned by Rtsne::Rtsne()
#' @param scores the calculated rareness scores by scRCD or other methods,
#' please remind that the higher rareness scores suggest higher possibility of being rare
#' @param cell.type a one-colunm matrix indicating the cell labels (by clustering or biomarkers) for each cell
#' @param frac.list a vector of desired fractions of rare cells for visualizations
#' @param frac a numeric value of the desired fraction, either frac.list or frac should be input
#'
#' @return A list object containing the x and y of tSNE, and the ggplot objects containing different visualizations depending on input arguments.
#'  If input a vector of fractions, returns the t-SNE visualization plots for each input rare cell fraction
#'  as well as the t-SNE plot for all the cells (cell types labelled in different colors);
#'  if input a single fraction value, the tSNE plot of rare cells at the single fraction,
#'  the visualization of uniformly sampled cells as the input fraction, and the visualization of all the cells.
#' @export
#'
#' @examples
scRCD.visualize<-function(tsne,scores,cell.type,frac.list=NULL,frac=NULL){
  library(ggplot2)
  df<-data.frame(cbind(tsne$Y),scores,cell.type)
  colnames(df)<-c("t_SNE_1","t_SNE_2","scores","cell_ontology_class")
  base_size = 12;
  theme.plt<-theme(legend.key.size = unit(1, "lines"),
                   legend.text = element_text(size = 0.55*base_size,
                                              color = "black"),
                   legend.title = element_text(size = 0.55*base_size,
                                               face = "bold",
                                               hjust = 0,
                                               color = "black"),
                   plot.title = element_text(size = 0.8*base_size,
                                             color = "black"),
                   axis.text.x = element_text(size = base_size*0.6, color = "black",
                                              lineheight = 0.9),
                   axis.text.y = element_text(size = base_size*0.6, color = "black",
                                              lineheight = 0.9))
  if(is.null(frac.list)){
    if(is.null(frac)){
      stop("Should input either frac.list or frac")
    }
    else{
      detect.num = floor(frac*dim(tsne$Y)[1])
      uniform.sample.index<-sample(1:dim(tsne$Y)[1],size = detect.num, replace = F)

      plt1<-ggplot(data = df[order(df[,"scores"],decreasing = T)[1:detect.num],],
                   aes(x = t_SNE_1, y = t_SNE_2))+geom_point(aes(color = cell_ontology_class ), size = 0.35)
      plt1<-plt1+xlim(min(tsne$Y[,1]),max(tsne$Y[,1]))+ylim(min(tsne$Y[,2]),max(tsne$Y[,2]))+ggtitle(paste0(detect.num, " rare cells detected"))+theme.plt

      plt2<-ggplot(data = df,
                   aes(x = t_SNE_1, y = t_SNE_2),)+geom_point(aes(color = cell_ontology_class ), size = 0.35)
      plt2<-plt2+xlim(min(tsne$Y[,1]),max(tsne$Y[,1]))+ylim(min(tsne$Y[,2]),max(tsne$Y[,2]))+ggtitle("Original data")+theme.plt

      plt3<-ggplot(data = df[uniform.sample.index,],
                   aes(x = t_SNE_1, y = t_SNE_2),)+geom_point(aes(color = cell_ontology_class ), size = 0.35)
      plt3<-plt3+xlim(min(tsne$Y[,1]),max(tsne$Y[,1]))+ylim(min(tsne$Y[,2]),max(tsne$Y[,2]))+ggtitle(paste0(detect.num," uniformly sampled cells"))+theme.plt

      return(list(tsne.out = tsne$Y,plt.rare = plt1, plt.all = plt2, plt.unif = plt3))
    }

  }
  else{
    plt.all<-ggplot(data = df,
                    aes(x = t_SNE_1, y = t_SNE_2),)+geom_point(aes(color = cell_ontology_class) , size = 0.35)
    plt.all<-plt.all+xlim(min(tsne$Y[,1]),max(tsne$Y[,1]))+ylim(min(tsne$Y[,2]),max(tsne$Y[,2]))+ggtitle("Original data")+theme.plt

    plt.list<-list()
    for(i in 1:length(frac.list)){
      detect.num = floor(frac.list[i]*dim(tsne$Y)[1])
      plt<-ggplot(data = df[order(df[,"scores"],decreasing = T)[1:detect.num],],
                  aes(x = t_SNE_1, y = t_SNE_2))+geom_point(aes(color = cell_ontology_class) , size = 0.35)
      plt<-plt+xlim(min(tsne$Y[,1]),max(tsne$Y[,1]))+ylim(min(tsne$Y[,2]),max(tsne$Y[,2]))+ggtitle(paste0(frac.list[i]*100, "% rare cells detected"))+theme.plt
      plt.list[[as.character(frac.list[i])]]<-plt
    }
    return(list(tsne.out = tsne$Y , plt.all.cell = plt.all , plt.frac.list = plt.list))
  }

}
