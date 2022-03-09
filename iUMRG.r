######################################################################
##### Multi-layered network-guided propagation modeling (iUMRG) ######
######################################################################
iUMRG <- function(seed,
                 network,
                 r){
  ##RWR function				   
  my_rand_walk <- function(W,p0,r){
    pt <- p0
    delta <- 1
    while(delta > 1e-10){
      pt1 <- (1-r)*W%*%pt+r*p0
      delta <- sum(abs(pt1 - pt))
      pt <- pt1}
    return(pt)
  }		
  ##data deal
  interact <- as.matrix(network)
  nodes <- unique(c(interact[,1],interact[,2]))
  nodes <- cbind(nodes,c(1:length(nodes)))
  rownames(nodes) <- nodes[,1]
  colnames(nodes) <- c("nodes","ID")
  Index1 <- match(interact[,1],nodes[,1])
  Index2 <- match(interact[,2],nodes[,1])
  Net_ID <- cbind(Index1,Index2)	
  
  m <- matrix(0,nrow(nodes),nrow(nodes))
  for(i in 1:nrow(Net_ID)) {
    n1 <- Net_ID[i,1]
    n2 <- Net_ID[i,2]
    m[n1,n2] <- 1
    m[n2,n1] <- 1
  }
  dev <- colSums(m)
  W <- t(t(m) / dev)
  
  seeds <- intersect(seed,rownames(nodes)) #39
  seedsnames <- seeds
  seeds <- nodes[seeds,2]
  seeds <- as.numeric(seeds)
  p0 <- rep(0,nrow(nodes))
  p0[seeds] <- 1/length(seeds) 
  r <- 0.4
  p_final <- my_rand_walk(W,p0,r)
  p_final <- cbind(nodes[,1],p_final)
  p_final <- p_final[-match(seedsnames,p_final[,1]),]
  rank <- order(order(1-as.numeric(p_final[,2])))
  p_final <- cbind(p_final,rank)
  result <- p_final[order(as.numeric(p_final[,3])),]
  return(result)
}


############################################
###          Run iUMRG framework        ####
############################################

HMMN <- read.table("Multi_layered_network.txt",header=T,sep="\t")
uvm_seed <- read.table("seed.txt",header=F,sep="\t")

result <- iUMRG(seed=uvm_seed[,1],HMMN,r=0.4)
write.csv(result[1:50,],"ranking_list_of_SG.csv",quote = F)






