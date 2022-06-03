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


dropoutpos= as.matrix(read.table("~path/data/dropoutpos_data.csv", header = FALSE, sep = ','))
for(i in 1:nrow(dropoutpos)){
  adjmat[dropoutpos[i,1],dropoutpos[i,2]]=9999
}
cell_gene_mat_tobeimputed=adjmat[1:nrow,(nrow+1):nnodes]

#To read the NN list
allcell_indx=as.matrix(read.table('~path/data/NNcell_data.csv',header=FALSE,sep=','))
allcell_indx=allcell_indx+1

# Code to find copula similarity between NN cells
library(copula)
datanew<-as.matrix(data1)
n=nrow(datanew)
theta=-0.5
fc=claytonCopula(theta,dim=2)
wcs_allcell=matrix(0,n,20) #number of neighbors
for(i in 1:n)
{
  for (j in 1:10)
  {nn=allcell_indx[i,(j+1)]
  u<-rbind(datanew[i,],datanew[nn,])
  u<-t(u)
  a=pCopula(u,fc)
  wcs_allcell[i,j]=round(mean(a),digits=4)
  }
}

#imputation code
allcell_indx=allcell_indx[,2:ncol(allcell_indx)]
cellind=which(cell_gene_mat_tobeimputed==9999, arr.ind=TRUE)[,1]
geneind=which(cell_gene_mat_tobeimputed==9999, arr.ind=TRUE)[,2]
imputed_val=vector()
cell_gene_mat_imputed=cell_gene_mat_tobeimputed
for(i in 1:length(cellind)){
  cand_cell=as.numeric(allcell_indx[cellind[i],])
  s=0
  for(j in 1:length(cand_cell)){
    if(cell_gene_mat_tobeimputed[cand_cell[j],geneind[i]]!=9999){
      s=s+as.numeric(wcs_allcell[cellind[i],j])*cell_gene_mat_tobeimputed[cand_cell[j],geneind[i]]
    }
  }
  imputed_val[i]=s/sum(wcs_allcell[cellind[i],])
  cell_gene_mat_imputed[cellind[i],geneind[i]]=imputed_val[i]
}

return(cell_gene_mat_imputed)

}