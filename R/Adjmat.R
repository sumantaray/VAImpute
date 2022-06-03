adjmat <- function (data)
{
  data1=data
  nrow=nrow(data1)
  ncol=ncol(data1)
  nnodes=nrow+ncol
  adjmat=matrix(0,nnodes,nnodes)
  for(i in 1:nrow)
  {
    adjmat[i,(nrow+1):nnodes]=data1[i,]
  }
  
  for(i in 1:ncol)
  {
    adjmat[(nrow+i),1:nrow]=data1[,i]
  }
  
  
  #To generate adjacency matrix with NA, 1 and 0
  ind= which(adjmat>0,arr.ind = TRUE)
  adjmat[ind]=1
  write.csv(adjmat,"~path/data/adjmat_data.csv",row.names = FALSE)
  #To generate test edges
  testind=which(adjmat[1:nrow,(ncol+1):nnodes]==0, arr.ind=TRUE)
  write.csv(testind,"~path/data/testedgepred_data.csv",row.names = FALSE)
  
  #return(norm_log_Data)
}
