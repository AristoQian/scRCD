#' scRCD main function: to calculate rareness scores of cells for scRNA-seq data
#'
#' Sorting the cells by the rareness scores, which are calculated from a graphical model via this function.
#' First, the anchor cells are randomly selected uniformly from all cells.
#' The bipartite (between anchor cells and all the cells) is constructed.
#' Then the affinity matrix is calculated as the bipartite multiply its transpose.
#' Finally the rareness scores are computed from a spectrum-clustering based optimization problem.
#' This steps are repeated until maximal repetition numbers achieved to avoid the effects of randomly sampling anchor cells.
#' The final scores are the averages of repeatly calculated scores.
#' For each repetition, the scores are calculated by calling rare.cell.scores()
#'
#' @param pc The top PCs selected from PCA or other dimension reduction methods,
#' the columns refer to features and rows refer to samples
#' @param anchor.frac The proportion of anchor cells with respect to all cells
#' @param alpha The percentage of quantile for the distance threshold for bipartite construction.
#' The distance threshold is chosen by the alpha-quantile of the Euclidean distances between cells.
#' The distance threshold is to control the bipartite distance between cells and anchor cells
#' such that only the edges of anchor cells with small distances will remain.
#' While constructing the bipartite and calculating the distance of each cell to anchor cells,
#' the anchor cells with the Euclidean distance (of PCs) from the cell lower than the threshold are selected and the bipartite distances
#' will be calculated by similar to k-NN graphs, the bipartite distance between the cell and the remaining anchor cells will be set to 0.
#' @param t the number of total repetitions to avoid the random effects of uniformly sampling anchor cells
#' @param n.cores the number of cores called for computation, 10 by default
#'
#' @return the rareness scores for all the cell samples, the high rareness score of a cell suggest
#' the "further" this cell is compared to majority cells and a high possibility of this cell to be rare
#' @export
#'
#' @example

scRCD<-function(pc, anchor.frac=0.15, alpha = 0.975, t = 10, n.cores = 10){
  scores.mat<-parallel::mcmapply(function(x){
    scores<-rare.cell.scores(pc = pc, anchor.frac = anchor.frac, alpha = alpha, times = x)
    return(scores)
  },1:t, mc.cores = n.cores)

  scores<-apply(scores.mat, 1, mean)
  return(scores = scores)
}



#' Rare cell detection for single cell RNA-seq data
#'
#' Sorting the cells by the rareness scores, which are calculated from a graphical model via this function.
#' First, the anchor cells are randomly selected uniformly from all cells.
#' The bipartite (between anchor cells and all the cells) is constructed.
#' Then the affinity matrix is calculated as the bipartite multiply its transpose.
#' Finally the rareness scores are computed from a spectrum-clustering based optimization problem.
#'
#' @param pc The top PCs selected from PCA or other dimension reduction methods,
#' the columns refer to features and rows refer to samples
#' @param anchor.frac The proportion of anchor cells with respect to all cells
#' @param alpha The percentage of quantile for the distance threshold for bipartite construction.
#' The distance threshold is chosen by the alpha-quantile of the Euclidean distances between cells.
#' The distance threshold is to control the bipartite distance between cells and anchor cells
#' such that only the edges of anchor cells with small distances will remain.
#' While constructing the bipartite and calculating the distance of each cell to anchor cells,
#' the anchor cells with the Euclidean distance (of PCs) from the cell lower than the threshold are selected and the bipartite distances
#' will be calculated by similar to k-NN graphs, the bipartite distance between the cell and the remaining anchor cells will be set to 0.
#' @param times The current times of repetition
#'
#' @return the rareness scores for all the cell samples, the high rareness score of a cell suggest
#' the "further" this cell is compared to majority cells and a high possibility of this cell to be rare
#' @export
#'
#' @examples

rare.cell.scores<-function(pc, anchor.frac, alpha, times){
  if(mode(pc)=="list"){
    pc<-as.matrix(pc)
  }
  print(paste0(paste0("The ",times),"-th repeat"))
  #the index of uniformly sampled cells
  uniform.sample.index<- sample(1:dim(pc)[1],size = floor(anchor.frac*dim(pc)[1]),replace = F)

  #distance matrix between anchors and other cells
  distance.cell2remain<-function(x,cell.num,red.dat){
    as.matrix(sapply(1:cell.num, function(y){
      if(x!=y){return(norm(as.matrix(red.dat[,x]-red.dat[,y]),"F"))}
      else{return(0)}
    }))
  }

  print("Constructing affinity matrix")
  H<-sapply(uniform.sample.index,function(x){
    distance.cell2remain(x = x, red.dat = t(pc), cell.num = dim(pc)[1])
  })

  #the hausdorff distance between sampled points and retained points
  cell.wise.min.distance<-apply(H,1,function(x){
    sort.x<-sort(x)
    if(sort.x[1]==0){
      return(sort.x[2])
    }
    else{
      return(sort.x[2])
    }
  })

  D.H<-stats::quantile(cell.wise.min.distance,alpha)
  #the bipartite based on Hausdorff distance and K-NN
  B<-H
  B[B>D.H]<-0


  Knn.dist.list<-lapply(1:dim(B)[1],function(x){
    if(length(which(B[x,]>0))==0){
      return(1)
    }
    else{
      dist.anc<-B[x,which(B[x,]>0)]
      sum.dist.anc<-sum(dist.anc)
      len.dist.anc<-length(which(B[x,]>0))
      return(as.vector((D.H-dist.anc)/(len.dist.anc*D.H-sum(dist.anc))))
    }

  })
  for(i in 1:dim(B)[1]){
    if(length(which(B[i,]>0))==0){
      B[i,which(H[i,]==min(H[i,]))]<- 1
    }
    else{
      B[i,which(B[i,]>0)]<-Knn.dist.list[[i]]
    }

  }

  print("Calculating rare cell score for each cell")
  #Laplacian, Degree matrix
  Lambda <- diag(apply(B,2,sum))[which(apply(B,2,sum) != 0),which(apply(B,2,sum) != 0)]
  W <- B[,which(apply(B,2,sum) != 0)]%*%t(B[,which(apply(B,2,sum) != 0)])
  L <- diag(apply(W,1,sum)) - W

  #compute the rareness score of each cell
  scores<-sapply(1:dim(L)[1], function(x){L[x,x]})
  scores<-(scores-min(scores))/(max(scores)-min(scores))
  return(scores)
}
