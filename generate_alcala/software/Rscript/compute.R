#######################################################################################################
# Script to Read output of cophylo.out, plot a cophylogeny, the network and compute network metrics   #
#######################################################################################################

require(ape)
require(igraph)

args <- commandArgs(trailingOnly = TRUE)

HT = as.numeric(args[1])
PT = as.numeric(args[2])
infiles  = args[3] 
outfiles = args[4] 
suffix = args[5]

##### USEFUL FONCTIONS #####

# create network from 2 phylogenies and list of association
HPphylo2graph <- function(Hi,Pi,hpedges,index,i=1,distHP=10){
  n=length(Hi$tip.label)
  np=length(Pi$tip.label)
  hpedges2 = hpedges
  hpedges2[,1] = hpedges2[,1]
  
  btH = which( branching.times(Hi)>max(branching.times(Pi)) ) +n
  ghp<-graph.edgelist(rbind(Hi$edge,Pi$edge+2*n-1,hpedges2,hpedges2[,2:1]) )
  V(ghp)$size=c(rep(3,n),rep(1,n-1),rep(3,np),rep(1,np-1)) #tips twice as big
  V(ghp)$type=c(rep(1,2*n-1),rep(2,2*np-1)) #4 types of vertices : (1, black) internal host vertex, (3, green) external host vertex,  (grey) (2, red) idem for parasites
  V(ghp)$type[1:np+2*n-1]=4
  V(ghp)$type[1:n]=3
  V(ghp)$color = c(rep("black",2*n-1),rep("grey",2*np-1))
  E(ghp)$type  = c(rep(1,nrow(Hi$edge)),rep(3,nrow(Pi$edge))     ,rep(2,2*nrow(hpedges)) ) #3 types of edges : (1, black) H-H, (grey) P-P and (2, red) H-P
  E(ghp)$color = c(rep("black",nrow(Hi$edge)),rep("grey",nrow(Pi$edge)),rep("red",2*nrow(hpedges)) )
  E(ghp)$weight = abs( c(Hi$edge.length,Pi$edge.length,rep(1,length(hpedges))) )
  
  ghp2 = delete.vertices(ghp,btH) #remove H nodes that are older than the P root
  layo = matrix(0,nrow=length(V(ghp)),ncol=2)
  layo[1:(2*n-1),1] = node.depth.edgelength(Hi)
  layo[1:(2*n-1),2] = node.height(Hi)/max(node.height(Hi))
  layo[(2*n):(2*n+2*np-2),1] = distHP + max(node.depth.edgelength(Hi)) + max(node.depth.edgelength(Pi)) - node.depth.edgelength(Pi)
  layo[(2*n):(2*n+2*np-2),2] = node.height(Pi)/max(node.height(Pi))
  ghp$layout = layo
  ghp2$layout = ghp$layout[-btH,]
  return(ghp2)
}
  



# compute network metrics
compstats <- function(ghp,nh,np){
  res=c()
  ##### start
  idh = which(V(ghp)$type==3)
  idp = which(V(ghp)$type==4)
  diam = diameter(ghp)
  print("diam") 
  
  evctmp = evcent(ghp)$vector
  evcHm  = mean(evctmp[idh])
  evcPm  = mean(evctmp[idp])
  evcHsd = sd(evctmp[idh])
  evcPsd = sd(evctmp[idp])
  print("evc") #slow
  
  coren    = graph.coreness(ghp)
  corenHm  = mean(coren[idh])
  corenPm  = mean(coren[idp])
  corenHsd = sd(coren[idh])
  corenPsd = sd(coren[idp])
  print("core")
  
  dens   = graph.density(ghp)
  print("dens")
  
  div    = graph.diversity(ghp)
  divHm  = mean(div[idh])
  divPm  = mean(div[idp])
  divHsd = sd(div[idh])
  divPsd = sd(div[idp])
  print("div")
  
  centev = centralization.evcent(ghp)$centralization
  
  degpowH = power.law.fit(degree(ghp)[idh])$alpha
  degpowP = power.law.fit(degree(ghp)[idp])$alpha
  print("degpow")
  
  simjHtmp = similarity.jaccard(ghp,vids=idh)
  simjPtmp = similarity.jaccard(ghp,vids=idp)
  simjHm   = mean( simjHtmp )
  simjPm   = sd( simjHtmp )
  simjHsd  = mean( simjPtmp )
  simjPsd  = sd( simjPtmp )
  print("simj")
  
  # add all new stats in 
  res = c(res,diam, evcHm,evcPm, evcHsd,evcPsd ,corenHm ,corenPm  ,corenHsd ,corenPsd ,dens  ,divHm ,divPm ,divHsd ,divPsd , 
          centev ,degpowH ,degpowP ,simjHm  ,simjPm  ,simjHsd ,simjPsd )
  names(res)=c(names(res), "diameter", "mean_eigenvaluecentrality_Hosts","mean_eigenvaluecentrality_Parasites", "sd_eigenvaluecentrality_Hosts","sd_eigenvaluecentrality_Parasites","mean_coreness_Hosts","mean_coreness_Parasites","sd_coreness_Hosts","sd_coreness_Parasites" ,"density","mean_diversity_Hosts" ,"mean_diversity_Parasites","sd_diversity_Hosts" ,"sd_diversity_Parasites",
               "centralized_eigenvaluecentrality" , "degreepower_Hosts" ,"degreepower_Parasites" ,"mean_similarityjaccard_Hosts"  ,"mean_similarityjaccard_Parasites"  ,"sd_similarityjaccard_Hosts"  ,"sd_similarityjaccard_Parasites"  ) #21 stats
  ## degree
  ddegtip<-function(g,typ="3"){
    tip=which(V(g)$type==typ)
    restmp=degree(g,v=tip)
    return(restmp)
  }
  # tip vertices
  #P
  deg_tp = ddegtip(ghp,typ="4")
  #H
  deg_th = ddegtip(ghp,typ="3")
  print("deg")
  
  res=c(res,mean(deg_th[deg_th>0]),sd(deg_th[deg_th>0]),
        mean(deg_tp[deg_tp>0]),sd(deg_tp[deg_tp>0])
  ) #2*2 stats
  names(res)[ (length(res)-3):length(res)]=c("mean_degree_Hosts","sd_degree_Hosts","mean_degree_Parasites" ,"sd_degree_Parasites") #+ 4 = 25
  ## assortativity
  assort=assortativity.degree(ghp)
  assortn=assortativity.nominal(ghp,types=V(ghp)$type)
  print("assort")
  
  res=c(res,assort,assortn) #2 stats
  names(res)[ (length(res)-1):length(res)]=c("assortativity", "modular_assortativity")
  ## modularity
  modu<-function(g){
    mem = as.numeric(as.factor(V(g)$type) )
    return(modularity(g,mem))
  }
  modu.walk<-function(g){
    mem = walktrap.community(g)
    return(modularity(g,mem$membership))
  }
  coms<-function(g){
    ug= as.undirected(g)
    mem = fastgreedy.community(ug)
    return(table(mem$membership))
  }
  modu=modu(ghp)       #nominal modulatirty
  print("modu")
  
  res=c(res,modu) #1 stats
  names(res)[ (length(res)-2):length(res)]=c("nominal_modularity" ) #+1 = 28
  ## matching index
  match<-function(g){
    typ  = V(g)$type
    typ2 = typ
    typ2[typ%in%c("1","3")]=1
    typ2[typ%in%c("2","4")]=2
    return(maximum.bipartite.matching(g,as.numeric(typ2) )$matching )
  }
  mat = match(ghp) 
  print("mat")
  
  res=c(res, mean(mat,na.rm=T) ) #1 stat
  names(res)[ (length(res)-1):length(res)]=c("mean_maximum_matching") #+1 = 29
  
  ## hubs auth
  limha=0.5
  hatip<-function(g,ha="hub",typ="3"){
    tip=which(V(g)$type==typ)
    if(ha=="hub"){
      res=hub.score(g)$vector[tip]
    }else{
      if(ha=="aut") res=authority.score(g)$vector[tip]
    }
    return(sum(res>limha))
  }
  ## tip vertices
  #P
  hub_tp = hatip(ghp,ha="hub",typ="4")
  #H
  hub_th = hatip(ghp,typ="3")
  #P
  aut_tp = hatip(ghp,ha="aut",typ="4")
  #H
  aut_th = hatip(ghp,ha="aut",typ="3")
  print("hub_aut")
  
  res=c(res, hub_th,aut_th,hub_tp,aut_tp ) #2*2 stats
  names(res)[(length(res)-3):length(res)]=c("number_hubs_Hosts", "number_authorities_Hosts","number_hubs_Parasitess", "number_authorities_Parasites" ) #+4 = 33
  #### end
  return(res) #33 stats
}


######  READ output files #######
#Read Host phylo
H <- list("edge"=NULL,"tip.label"=NULL,"edge.length"=NULL,"Nnode"=NULL)#H phylogeny
H$edge = as.matrix(read.table(paste(infiles,"host_edge.txt",sep="")))
nh = nrow(H$edge)/2 +1 
H$edge.length = scan(paste(infiles,"host_edgelength.txt",sep=""))
H$Nnode       = nh-1
H$tip.label   = as.character(1:(H$Nnode+1) )
class(H)<-"phylo"

#Read Parasite phylo
P<-list("edge"=NULL,"tip.label"=NULL,"edge.length"=NULL,"Nnode"=NULL)#P phylogeny
P$edge        = as.matrix(read.table(paste(outfiles,"edges_0.txt",sep="")) )
P$edge.length = scan(paste(outfiles,"edgelength_0.txt",sep=""))
P$Nnode       = scan(paste(outfiles,"nnode_0.txt",sep=""))
P$tip.label   = as.character(1:(P$Nnode+1) )
class(P)<-"phylo"

#Read H-P assoc
HP = matrix(scan(paste(outfiles,"hpassoc_0.txt",sep="")),ncol=2,byrow=T)

#remove extinct Parasite lineages
maxP = HP[nrow(HP),1] - nh*2+1
nPtips = length(P$tip.label)
if(nPtips>maxP) P = drop.tip( P , (maxP+1):nPtips ) 


#Create Network
ghp = HPphylo2graph(H,P,HP)

#Compute Network metrics
stats = compstats(ghp)

##### PLOT #####

# plot cophylo
pdf(paste(suffix,"tanglegram.pdf",sep=""),h=3.5,w=3.5*3)
par( mai = c(1, 0, 0.5, 0) )#par(mfrow=c(1,3))
layout(matrix(c(rep(1,3),2,rep(3,3)), nrow = 1, byrow = TRUE))
plot(H,col=1, show.tip.label = TRUE,x.lim=c(0,HT),main="Host phylogeny")
print( HT )
axis(1,at=seq(0,HT,length.out=6),labels=rev(seq(0,HT,length.out=6)) )
plot(-1,-1,xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="",main="associations")
segments(rep(0,nrow(HP)) , (HP[,2]-1)/H$Nnode, rep(1,nrow(HP)) , (HP[,1]-(nh-1)*2-2)/P$Nnode,col=rgb(1,0,0,1))
plot(P,edge.color="gray",tip.color="gray",direction="leftwards",cex=1,show.tip.label = TRUE,x.lim=c(0,HT),main="Parasite Phylogeny")
axis(1,at=seq(0,PT,length.out=6)+HT-PT,labels=seq(0,PT,length.out=6) )
dev.off()

# plot network
lay = layout.kamada.kawai(ghp)
pdf(paste(suffix,"network.pdf",sep=""),h=7,w=7)
plot(ghp,layout=lay,vertex.label=NA,edge.arrow.mode="-")
dev.off()

# write stats
write( rbind(names(stats),stats), file=paste(suffix,"stats.txt",sep=""),ncolumns=2)


