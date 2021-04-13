library(ape)
library(phytools)
library(castor)
library(phangorn)
library(caper)

#Process data:
mv <- read.table("Mean_value.csv", sep="\t",head=TRUE,check.names=FALSE)
species_name <- mv[,1]
bodysize <- mv[,5]
pmax <- mv[,6]
tr <- c(read.tree("tree1.tre"),read.tree("tree2.tre"))

char_subset <- c("P_stride/PL","P_gauge/PL","P_minus_M_gauge/PL",
                 "GAD/PL","P_orientation","M_minus_P_orientation")
all_char <- c("P_stride/PL", "P_pace/PL","P_gauge/PL","P_minus_M_gauge/PL",
              "GAD/PL", "PM_distance/PL","P_pace_ang","M_pace_ang",
              "P_orientation","M_orientation","M_minus_P_orientation")
log_char <- c("PL (pes length)","pIV length")
specimen_file <- "tracks_specimens.csv"
treefile <- "trees.csv"
specimens <- read.table(specimen_file, sep="\t",head=TRUE,check.names=FALSE)
sp <- split(specimens,specimens$`No Species`)
species <- do.call(rbind,lapply(sp,function(x)colMeans(x[,colnames(specimens)%in%c(all_char,log_char)],na.rm=TRUE)))
trees <- list()
for(i in 1:length(tr)){
  trees[[i]] <- multi2di(tr[[i]],tips=TRUE)
  trees[[i]]$edge.length[trees[[i]]$edge.length==0]<-0.1
}
dat <- list(specimens=specimens,species=species,trees=trees)
n_char <- length(char_subset)
n_totalchar <- length(all_char)
n_tax <- length(unique(specimens[,1]))
n_tree <- length(tr)

# Mean Tracks
iA <- list()
for(j in 1:n_tax){
  all <- specimens[specimens[,1]==j,colnames(specimens)%in%all_char]
  if(nrow(all)>1){
    crv <- sapply(1:n_char,function(x){val <- all[,colnames(all)%in%char_subset[x]];(val-mean(val,na.rm=T))/sd(val,na.rm=T)})
    iCRV <- crv[!rowSums(is.na(all)),]
    inst <- all[!rowSums(is.na(all)),]
    if(nrow(inst)>1){
      iA[[j]] <- inst[which.min(apply(iCRV,1,function(x)sqrt(sum(x^2)))),]
    }else{
      iA[[j]] <- inst
    }
  }else{
    iA[[j]] <- all
  }
}
meantracks <- do.call(rbind,iA)
meantracks$taxa <- 1:n_tax

# ACE based on mean tracks
res <- list()
res[[1]]<-list()
for(i in 1:n_totalchar){
  res[[1]][[i]] <- list()
  for(j in 1:n_tree){
    res[[1]][[i]][[j]] <- reconstruct(sapply(trees[[j]]$tip.label,function(x)meantracks[meantracks$taxa==x,colnames(meantracks)==all_char[i]]),trees[[j]],method="ML")
  }
}
amt <- list(ACE=res, MEAN=meantracks)

save(dat,amt,file="savepoint0.Rdata")

#Plot ACE of last 11 characters
TR <- dat$trees
dir.create("with_meantrack")
setwd("with_meantrack")
for(i in 1:n_totalchar){
    tip_values <- amt$MEAN[,colnames(amt$MEAN)%in%all_char[i]]
    for(j in 1:n_tree){
      all_ace <- lapply(amt$ACE,function(x)sapply(x,function(y)y[[j]]$ace))
      node_values <- apply(do.call(cbind,lapply(all_ace, function(x)x[,i])),1,mean)
      TR[[j]]$tip.label <- species_name
      names(tip_values) <- species_name
      pdf(sprintf("%s_tree%i.pdf",gsub("/","-",all_char[i]),j),bg="transparent",width=16,height=16)
      par(cex=1)
      contMap(TR[[j]],tip_values,method="user",anc.states=node_values,direction="upwards")
      dev.off()
    }
}

mt <- do.call(rbind,lapply(rownames(amt$MEAN),function(x)dat$specimens[dat$specimens[,3]==x,]))

#Make ACE for the log characters
pmax_ace <- reconstruct(sapply(dat$trees[[1]]$tip.label,function(x)log(pmax[mv[,2]==x],10)),dat$trees[[1]],method="ML")
pmax_ace$ace <- 10^pmax_ace$ace
pmax_ace$CI95 <- 10^pmax_ace$CI95
pmax_ace2 <- reconstruct(sapply(dat$trees[[2]]$tip.label,function(x)log(pmax[mv[,2]==x],10)),dat$trees[[2]],method="ML")
pmax_ace2$ace <- 10^pmax_ace2$ace
pmax_ace2$CI95 <- 10^pmax_ace2$CI95
pmax_ace <- list(pmax_ace, pmax_ace2)

pl_ace <- reconstruct(sapply(dat$trees[[1]]$tip.label,function(x)log(mt[mt[,1]==x,"PL (pes length)"],10)),dat$trees[[1]],method="ML")
pl_ace$ace <- 10^pl_ace$ace
pl_ace$CI95 <- 10^pl_ace$CI95
pl_ace2 <- reconstruct(sapply(dat$trees[[2]]$tip.label,function(x)log(mt[mt[,1]==x,"PL (pes length)"],10)),dat$trees[[2]],method="ML")
pl_ace2$ace <- 10^pl_ace2$ace
pl_ace2$CI95 <- 10^pl_ace2$CI95
pl_ace <- list(pl_ace, pl_ace2)

piv_ace <- reconstruct(sapply(dat$trees[[1]]$tip.label,function(x)log(mt[mt[,1]==x,"pIV length"],10)),dat$trees[[1]],method="ML")
piv_ace$ace <- 10^piv_ace$ace
piv_ace$CI95 <- 10^piv_ace$CI95
piv_ace2 <- reconstruct(sapply(dat$trees[[2]]$tip.label,function(x)log(mt[mt[,1]==x,"pIV length"],10)),dat$trees[[2]],method="ML")
piv_ace2$ace <- 10^piv_ace2$ace
piv_ace2$CI95 <- 10^piv_ace2$CI95
piv_ace <- list(piv_ace, piv_ace2)

#Save all ACE results
t1A<-sapply(amt$ACE[[1]],function(x)x[[1]]$ace)
t1u<-sapply(amt$ACE[[1]],function(x)x[[1]]$CI95[,2])
t1l<-sapply(amt$ACE[[1]],function(x)x[[1]]$CI95[,1])
M1<-matrix(sprintf("%.3f (%.3f-%.3f)",t1A,t1l,t1u),nrow=nrow(t1A),ncol=ncol(t1A))

t2A<-sapply(amt$ACE[[1]],function(x)x[[2]]$ace)
t2u<-sapply(amt$ACE[[1]],function(x)x[[2]]$CI95[,2])
t2l<-sapply(amt$ACE[[1]],function(x)x[[2]]$CI95[,1])
M2<-matrix(sprintf("%.3f (%.3f-%.3f)",t2A,t2l,t2u),nrow=nrow(t2A),ncol=ncol(t2A))

M1 <- as.data.frame(M1)
M2 <- as.data.frame(M2)
colnames(M1)<-colnames(M2)<-all_char
M1$piv <- sprintf("%.3f (%.3f-%.3f)", piv_ace[[1]]$ace,piv_ace[[1]]$CI95[,1],piv_ace[[1]]$CI95[,2])
M2$piv <- sprintf("%.3f (%.3f-%.3f)", piv_ace[[2]]$ace,piv_ace[[2]]$CI95[,1],piv_ace[[2]]$CI95[,2])
M1$pl <- sprintf("%.3f (%.3f-%.3f)", pl_ace[[1]]$ace,pl_ace[[1]]$CI95[,1],pl_ace[[1]]$CI95[,2])
M2$pl <- sprintf("%.3f (%.3f-%.3f)", pl_ace[[2]]$ace,pl_ace[[2]]$CI95[,1],pl_ace[[2]]$CI95[,2])
M1$pmax <- sprintf("%.3f (%.3f-%.3f)", pmax_ace[[1]]$ace,pmax_ace[[1]]$CI95[,1],pmax_ace[[1]]$CI95[,2])
M2$pmax <- sprintf("%.3f (%.3f-%.3f)", pmax_ace[[2]]$ace,pmax_ace[[2]]$CI95[,1],pmax_ace[[2]]$CI95[,2])
M1$node<-sapply(Descendants(dat$trees[[1]],23:43),function(x)paste(dat$trees[[1]]$tip.label[x],collapse="+"))
M2$node<-sapply(Descendants(dat$trees[[2]],23:43),function(x)paste(dat$trees[[2]]$tip.label[x],collapse="+"))

#write.table(M1,"ACE1.csv",sep="\t",row.names=FALSE)
#write.table(M2,"ACE2.csv",sep="\t",row.names=FALSE)

# Make plots for the log characters
for(i in 1:2){
  tip.values1 <- sapply(dat$trees[[i]]$tip.label,function(x)mt[mt[,1]==x,"pIV length"])
  tip.values2 <- sapply(dat$trees[[i]]$tip.label,function(x)mt[mt[,1]==x,"PL (pes length)"])
  tip.values3 <- sapply(dat$trees[[i]]$tip.label,function(x)pmax[mv[,2]==x])
  TR <- dat$trees[[i]]
  TR$tip.label <- species_name
  names(tip.values1) <- species_name
  names(tip.values2) <- species_name
  names(tip.values3) <- species_name
  pdf(sprintf("pIV_tree%i.pdf",i),bg="transparent",width=16,height=16)
  par(cex=1)
  contMap(TR,tip.values1,method="user",anc.states=piv_ace[[i]]$ace,direction="upwards")
  dev.off()
  pdf(sprintf("pl_tree%i.pdf",i),bg="transparent",width=16,height=16)
  par(cex=1)
  contMap(TR,tip.values2,method="user",anc.states=pl_ace[[i]]$ace,direction="upwards")
  dev.off()
  pdf(sprintf("pmax_tree%i.pdf",i),bg="transparent",width=16,height=16)
  par(cex=1)
  contMap(TR,tip.values3,method="user",anc.states=pmax_ace[[i]]$ace,direction="upwards")
  dev.off()
}

#Body Mark attempt
bs_ace <- hsp_independent_contrasts(dat$trees[[1]],sapply(dat$trees[[1]]$tip.label,function(x)bodysize[mv[,2]==x]))$states
bs_ace2 <- hsp_independent_contrasts(dat$trees[[2]],sapply(dat$trees[[2]]$tip.label,function(x)bodysize[mv[,2]==x]))$states
write.table(cbind(bs_ace,c(1:22,sapply(Descendants(dat$trees[[1]],23:43),function(x)paste(dat$trees[[1]]$tip.label[x],collapse="+")))),file="bs_ace.csv",sep="\t")
write.table(cbind(bs_ace2,c(1:22,sapply(Descendants(dat$trees[[2]],23:43),function(x)paste(dat$trees[[2]]$tip.label[x],collapse="+")))),file="bs_ace2.csv",sep="\t")

bi_ace <- hsp_mk_model(dat$trees[[1]],sapply(dat$trees[[1]]$tip.label,function(x)1+as.integer(bodysize[mv[,2]==x]>0)),Nstates=2)
bi_ace2 <- hsp_mk_model(dat$trees[[2]],sapply(dat$trees[[2]]$tip.label,function(x)1+as.integer(bodysize[mv[,2]==x]>0)),Nstates=2)

t1 <- dat$trees[[1]]
t1$tip.label <- sapply(t1$tip.label,function(x)species_name[as.integer(x)])
pdf("bodymarks_binary.pdf",h=16,w=16)
par(cex=1)
plot(t1,label.offset=1)
tiplabels(pie=bi_ace$likelihood[1:22,],cex=0.5,piecol=c("black","white"))
nodelabels(pie=bi_ace$likelihood[23:43,],cex=0.5,piecol=c("black","white"))
dev.off()

t2 <- dat$trees[[2]]
t2$tip.label <- sapply(t2$tip.label,function(x)species_name[as.integer(x)])
pdf("bodymarks_binary2.pdf",h=16,w=16)
par(cex=1)
plot(t2,label.offset=1)
tiplabels(pie=bi_ace2$likelihood[1:22,],cex=0.5,piecol=c("black","white"))
nodelabels(pie=bi_ace2$likelihood[23:43,],cex=0.5,piecol=c("black","white"))
dev.off()

bs <- cut(bodysize,c(0,.25,.5,.75,1),right=TRUE,include.lowest=TRUE)
bmp1_ace1 <- hsp_mk_model(dat$trees[[1]],
                          sapply(dat$trees[[1]]$tip.label,function(x)as.integer(bs[mv[,2]==x])),
                          rate_model="SUEDE",Nstates=4)
bmp1_ace2 <- hsp_mk_model(dat$trees[[2]],
                          sapply(dat$trees[[2]]$tip.label,function(x)as.integer(bs[mv[,2]==x])),
                          rate_model="SUEDE",Nstates=4)

bs <- as.integer(cut(bodysize,c(0,1e-5,1/3,1,2),right=FALSE))
bmp2_ace1 <- hsp_mk_model(dat$trees[[1]],
                          sapply(dat$trees[[1]]$tip.label,function(x)as.integer(bs[mv[,2]==x])),
                          rate_model="SUEDE",Nstates=4)
bmp2_ace2 <- hsp_mk_model(dat$trees[[2]],
                          sapply(dat$trees[[2]]$tip.label,function(x)as.integer(bs[mv[,2]==x])),
                          rate_model="SUEDE",Nstates=4)
pdf("bodymarks_none_few_all1.pdf",h=16,w=16)
par(cex=1)
plot(t1,label.offset=1)
tiplabels(pie=bmp2_ace1$likelihood[1:22,],cex=0.5,piecol=c("black","grey50","grey80","white"))
nodelabels(pie=bmp2_ace1$likelihood[23:43,],cex=0.5,piecol=c("black","grey50","grey80","white"))
dev.off()

pdf("bodymarks_none_few_all2.pdf",h=16,w=16)
par(cex=1)
plot(t2,label.offset=1)
tiplabels(pie=bmp2_ace2$likelihood[1:22,],cex=0.5,piecol=c("black","grey50","grey80","white"))
nodelabels(pie=bmp2_ace2$likelihood[23:43,],cex=0.5,piecol=c("black","grey50","grey80","white"))
dev.off()


#Comparisons with Node K
nodeK <- as.numeric(gsub("^([0-9.-]+) .+$","\\1",M1[M1[,15]=="15+16+17+18+19+20+21+22",1:14]))
names(nodeK)<-colnames(M1)[1:14]
nodeK <- nodeK[1:11]
sp <- dat$specimens[,names(nodeK)[1:11]]
#sp <- sp[rowSums(is.na(sp))<1,]
write.table(data.frame(species=dat$specimens$`Species`,specimen=dat$specimens$`Specimen name`,dist_to_K=sqrt(colSums(((t(sp)-nodeK)/apply(sp,2,var,na.rm=T))^2))),file="dist2k_specimens.csv",sep="\t")

nodeK <- as.numeric(gsub("^([0-9.-]+) .+$","\\1",M2[M2[,15]=="15+17+18+19+20+21+22+16",1:14]))
names(nodeK)<-colnames(M2)[1:14]
nodeK <- nodeK[1:11]
sp <- dat$specimens[,names(nodeK)[1:11]]
write.table(data.frame(species=dat$specimens$`Species`,specimen=dat$specimens$`Specimen name`,dist_to_K=sqrt(colSums(((t(sp)-nodeK)/apply(sp,2,var,na.rm=T))^2))),file="dist2k_specimens__2.csv",sep="\t")

#Bootstrap based on polymorphism
n_trials <- 1e4
n_char <- length(all_char)
library(doSNOW)
cl <- parallel::makeCluster(2)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_trials, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
set.seed(20180822)
aas <- foreach(i = seq_len(n_trials), .options.snow = opts) %dopar% {
    R<-list()
    for(j in 1:n_char){
      R[[j]]<-list()
      sp <- dat$specimens[!is.na(dat$specimens[,colnames(dat$specimens)==all_char[j]]),]
      specimen_set <- sapply(split(sp$`No Specimen`,sp$`No Species`),function(x)ifelse(length(x)>1,sample(x,1),x))
      values <- sapply(specimen_set,function(x)sp[sp$`No Specimen`==x,colnames(sp)==all_char[j]])
      for(k in 1:n_tree){
        R[[j]][[k]]<-ape::reconstruct(values,dat$trees[[k]],method="ML")
      }
    }
    R
  }
close(pb)
stopCluster(cl)

d <- dat$specimens
d$`PL (pes length)` <- log(d$PL,10)
d$`pIV length` <- log(d$pIV,10)
cl <- parallel::makeCluster(2)
registerDoSNOW(cl)
pb <- txtProgressBar(max = n_trials, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
set.seed(20180822)
n_char <- length(log_char)
aas_log <- foreach(i = seq_len(n_trials), .options.snow = opts) %dopar% {
  R<-list()
  for(j in 1:n_char){
    R[[j]]<-list()
    sp <- d[!is.na(d[,colnames(d)==log_char[j]]),]
    specimen_set <- sapply(split(sp$`No Specimen`,sp$`No Species`),function(x)ifelse(length(x)>1,sample(x,1),x))
    values <- sapply(specimen_set,function(x)sp[sp$`No Specimen`==x,colnames(sp)==log_char[j]])
    for(k in 1:n_tree){
      R[[j]][[k]]<-ape::reconstruct(values,dat$trees[[k]],method="ML")
    }
  }
  R
}
close(pb)
stopCluster(cl)

#Making plots
n_char <- length(all_char)
for(j in 1:n_tree){
  all_ace <- lapply(aas,function(x)sapply(x,function(y)y[[j]]$ace))
  for(i in 1:n_char){
    char <- do.call(cbind,lapply(all_ace, function(x)x[,i]))
    rownames(char) <- dat$trees[[j]]$node.label
    char <- char[rownames(char)!="",]
    n_node - nrow(char)
    pdf((sprintf("%s_tree%i_distrib.pdf",gsub("/","-",all_char[i]),j)),h=8,w=8)
    layout(matrix(n_node:1,ncol=1),height=c(rep(1,n_node-1),1.5))
    for(k in 1:n_node){
      if(k==1){
        par(mar=c(3,15,0,1))
      }else if (k==n_node){
        par(mar=c(0,15,1,1))
      }else{par(mar=c(0,15,0,1))}
      de <- density(char[k,])
      plot(de,yaxs="i",
           xaxs="i",xlim=range(char),ann=FALSE,axes=FALSE)
      q <- quantile(char[k,],c(.25,.75))
      px <- c(q[1],de$x[de$x>q[1]&de$x<q[2]],q[2])
      py <- c(0,de$y[de$x>q[1]&de$x<q[2]],0)
      polygon(cbind(px,py),col="grey50")
      abline(v=mean(char[k,]),col="red")
      box(lwd=2)
      if(k==1){axis(1)}
      mtext(rownames(char)[k],2,1,las=2,cex=0.6)
    }
    dev.off()
    
    pdf((sprintf("%s_tree%i_boxplot.pdf",gsub("/","-",all_char[i]),j)),h=4,w=8)
    par(mar=c(3,15,1,1))
    plot(NA,ann=FALSE,axes=FALSE,ylim=c(0,n_node+1),xlim=range(char),xaxs="i",yaxs="i")
    boxplot(t(char),pars=list(ann=FALSE,axes=FALSE),pch=19,cex=0.7,add=TRUE,horizontal=TRUE)
    box(lwd=2)
    axis(1)
    axis(2,las=2,labels=rownames(char),at=1:n_node)
    dev.off()
  }
}

n_char <- length(log_char)
for(j in 1:n_tree){
  all_ace <- lapply(aas_log,function(x)sapply(x,function(y)y[[j]]$ace))
  for(i in 1:n_char){
    char <- do.call(cbind,lapply(all_ace, function(x)x[,i]))
    rownames(char) <- dat$trees[[j]]$node.label
    char <- char[rownames(char)!="",]
    n_node <- nrow(char)
    pdf((sprintf("%s_tree%i_distrib.pdf",gsub("/","-",log_char[i]),j)),h=8,w=8)
    layout(matrix(n_node:1,ncol=1),height=c(rep(1,n_node-1),1.5))
    for(k in 1:n_node){
      if(k==1){
        par(mar=c(3,15,0,1))
      }else if (k==n_node){
        par(mar=c(0,15,1,1))
      }else{par(mar=c(0,15,0,1))}
      de <- density(char[k,])
      plot(de,yaxs="i",xaxs="i",xlim=range(char),ann=FALSE,axes=FALSE)
      q <- quantile(char[k,],c(.25,.75))
      px <- c(q[1],de$x[de$x>q[1]&de$x<q[2]],q[2])
      py <- c(0,de$y[de$x>q[1]&de$x<q[2]],0)
      polygon(cbind(px,py),col="grey50")
      abline(v=mean(char[k,]),col="red")
      box(lwd=2)
      if(k==1){axis(1, at=log(c(10,20,50,100,200,500,1000),10),label=c(10,20,50,100,200,500,1000))}
      mtext(rownames(char)[k],2,1,las=2,cex=0.6)
    }
    dev.off()
    
    pdf((sprintf("%s_tree%i_boxplot.pdf",gsub("/","-",log_char[i]),j)),h=4,w=8)
    par(mar=c(3,15,1,1))
    plot(NA,ann=FALSE,axes=FALSE,ylim=c(0,n_node+1),xlim=range(char),xaxs="i",yaxs="i")
    boxplot(t(char),pars=list(ann=FALSE,axes=FALSE),pch=19,cex=0.7,add=TRUE,horizontal=TRUE)
    box(lwd=2)
    axis(1, at=log(c(10,20,50,100,200,500,1000),10),label=c(10,20,50,100,200,500,1000))
    axis(2,las=2,labels=rownames(char),at=1:n_node)
    dev.off()
  }
}

#Correlation between characters
crunchyA <- list()
AMT <- amt$MEAN
colnames(AMT) <- gsub("[/() ]",".",colnames(AMT))
for(i in 1:2){
  dat$trees[[i]]$node.label[dat$trees[[i]]$node.label==""] <- paste("poly",1:3,sep="")
  tracks <- comparative.data(dat$trees[[i]],AMT,"taxa")
  crunchyA[[i]]<-list()
  for(j in 1:11){
    crunchyA[[i]][[j]]<-lapply(colnames(AMT)[-c(j,ncol(AMT))], function(x)crunch(reformulate(x,response=colnames(AMT)[j],intercept=FALSE),data=tracks))
  }
}
r2<-do.call(rbind,lapply(crunchyA[[1]],function(x)sapply(x,function(y)summary(y)$adj.r.squared)))
pv<-do.call(rbind,lapply(crunchyA[[1]],function(x)sapply(x,function(y)summary(y)$coef[4])))
write.table(pv,file="crunchy1.csv",sep="\t")
write.table(r2,file="crunchy1_r2.csv",sep="\t")

r2<-do.call(rbind,lapply(crunchyA[[2]],function(x)sapply(x,function(y)summary(y)$adj.r.squared)))
pv<-do.call(rbind,lapply(crunchyA[[2]],function(x)sapply(x,function(y)summary(y)$coef[4])))
write.table(pv,file="crunchy2.csv",sep="\t")
write.table(r2,file="crunchy2_r2.csv",sep="\t")


#Akaike
dev.BM <- function(p) {
  nb.tip <- Ntip(phy)
  nb.node <- Nnode(phy)
  tip <- phy$edge[, 2] <= nb.tip
  if (p[1] < 0) return(1e+100)
  x1 <- p[-1][phy$edge[, 1] - nb.tip]
  x2 <- numeric(length(x1))
  x2[tip] <- x[phy$edge[tip, 2]]
  x2[!tip] <- p[-1][phy$edge[!tip, 2] - nb.tip]
  -2 * (-sum((x1 - x2)^2/phy$edge.length)/(2 * p[1]) - nb.node * log(p[1]))
}
ll <- sig <- list()
for(j in 1:n_tree){
  phy <- dat$trees[[j]]
  a1 <- list()
  for(i in seq_along(all_char)){
    if(j==1)  pp <- t1A[,i]
    if(j==2)  pp <- t2A[,i]
    x <- amt$MEAN[,all_char[i]]
    a1[[i]] <- nlm(dev.BM,c(1,pp),hessian=TRUE)
  }
  ll[[j]]<-sapply(a1,function(x)x$minimum)
}
aic <- sapply(ll,function(x)-2*sum(x)+2)
daic <- aic-min(aic)
