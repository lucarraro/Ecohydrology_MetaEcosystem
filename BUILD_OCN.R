rm(list=ls())
library(OCNet)
library(R.matlab)
library(spam)
library(terra)
library(rivnet)

plotFigures <- FALSE # set to TRUE to produce figures

if (!file.exists("utilities/OCN_480.rda")){
  OCN_480 <- NULL
  save(OCN_480,file="utilities/OCN_480.rda")
set.seed(1)
OCN_480 <- create_OCN(480, 480, displayUpdates = 2, nUpdates = 200,
                      cellsize=100, outletPos=round(480/10),
                      typeInitialState = "V",
                      initialNoCoolingPhase = 0.1, coolingRate = 0.5)
OCN_480 <- landscape_OCN(OCN_480, zMin=200, slope0 = 0.005)
OCN_480 <- aggregate_OCN(OCN_480, thrA=5e6, maxReachLength = 4000, equalizeLength=TRUE)
OCN_480 <- paths_OCN(OCN_480, level="AG", includePaths = TRUE)
OCN_480 <- rivergeometry_OCN(OCN_480, widthMax = 19.2, velocityMax = 1.25, depthMax = 2.4)
save(OCN_480,file="utilities/OCN_480.rda")
} else {load('utilities/OCN_480.rda')}

hw <- which(colSums(OCN_480$AG$W)==0)
pathLengths <- OCN_480$AG$downstreamPathLength[hw,OCN_480$AG$outlet]

nNodesPath <- numeric(length(hw))
for (i in 1:length(hw)){
nNodesPath[i] <- length(OCN_480$AG$downstreamPath[[hw[i]]][[OCN_480$AG$outlet]])}
ss <- sort(pathLengths, index.return=T)

ind_mpl <- 191 
pathLong <- OCN_480$AG$downstreamPath[[which.max(OCN_480$AG$downstreamPathLength[,OCN_480$AG$outlet])]][[OCN_480$AG$outlet]]
pathMedian <- OCN_480$AG$downstreamPath[[ind_mpl]][[OCN_480$AG$outlet]]

vec480 <- numeric(OCN_480$AG$nNodes)
vec480[pathMedian] <- vec480[pathMedian] + 2
vec480[pathLong] <- vec480[pathLong] + 3

if(!file.exists('utilities/OCN_144.mat')){
writeMat('utilities/OCN_480.mat',W=as.matrix(OCN_480$AG$W), downNode=OCN_480$AG$downNode, A=OCN_480$AG$A, Z=OCN_480$AG$Z, 
         width=OCN_480$AG$width, depth=OCN_480$AG$depth, velocity=OCN_480$AG$velocity, slope=OCN_480$AG$slope,
         As=OCN_480$SC$ALocal, L=OCN_480$AG$leng, streamOrder=OCN_480$AG$streamOrder,
         distToOutlet=OCN_480$AG$downstreamPathLength[,OCN_480$AG$outlet],
         pathLong=pathLong, pathMedian=pathMedian)
}



if (!file.exists("utilities/OCN_144.rda")){
  OCN_144 <- NULL
  save(OCN_144,file="utilities/OCN_144.rda")
  set.seed(11)
  OCN_144 <- create_OCN(144, 1600, displayUpdates = 2, nUpdates = 200,
                        typeInitialState = "V",
                        cellsize=100, outletPos=round(144/10),
                        initialNoCoolingPhase = 0.1, coolingRate = 0.5)
  OCN_144 <- landscape_OCN(OCN_144, zMin=200, slope0 = 0.005)
  OCN_144 <- aggregate_OCN(OCN_144, thrA=5e6, maxReachLength = 4000, equalizeLength=TRUE)
  OCN_144 <- paths_OCN(OCN_144, level="AG", includePaths = TRUE)
  OCN_144 <- rivergeometry_OCN(OCN_144, widthMax = 19.2, velocityMax = 1.25, depthMax = 2.4)
  save(OCN_144,file="utilities/OCN_144.rda")
} else {load('utilities/OCN_144.rda')}

  ind_mpl <- 330
  pathLong=OCN_144$AG$downstreamPath[[which.max(OCN_144$AG$downstreamPathLength[,OCN_144$AG$outlet])]][[OCN_144$AG$outlet]]
  pathMedian=OCN_144$AG$downstreamPath[[ind_mpl]][[OCN_144$AG$outlet]]
  vec144 <- numeric(OCN_144$AG$nNodes)
  vec144[pathMedian] <- vec144[pathMedian] + 2
  vec144[pathLong] <- vec144[pathLong] + 3

  if(!file.exists('utilities/OCN_144.mat')){
  writeMat('utilities/OCN_144.mat',W=as.matrix(OCN_144$AG$W), downNode=OCN_144$AG$downNode, A=OCN_144$AG$A, Z=OCN_144$AG$Z, 
           width=OCN_144$AG$width, depth=OCN_144$AG$depth, velocity=OCN_144$AG$velocity, slope=OCN_144$AG$slope,
           As=OCN_144$SC$ALocal, L=OCN_144$AG$leng, streamOrder=OCN_144$AG$streamOrder,
           distToOutlet=OCN_144$AG$downstreamPathLength[,OCN_144$AG$outlet],
           pathLong=pathLong, pathMedian=pathMedian)
  }

# Plot Figures 1, S1a, S7 ####  
if(plotFigures){
# Fig. S1  
pdf(file="Fig1_draft.pdf",width=24/2.54,height=18/2.54)
par(mfrow=c(1,2))
plot(vec480,OCN_480,colPalette=hcl.colors(1000,'Rocket',rev=T), addLegend=F)
sbar(d=1000,xy=c(min(OCN_480$FD$X), max(OCN_480$FD$Y)))
plot(vec144,OCN_144,colPalette=hcl.colors(1000,'Rocket',rev=T), addLegend=F)
sbar(d=1000,xy=c(min(OCN_144$FD$X), max(OCN_144$FD$Y)))
dev.off()

# Fig. S1a
bioDF <- readMat('utilities/bioDF_144_N5D5_def.mat')
bioDF <- bioDF$bioDF
medBioDF <- apply(bioDF,1,median)

ff<- which(vec144==3 | vec144==5)
gg<- sort(OCN_144$AG$A[ff],index.return=T)
ff<- ff[gg$ix]

pdf(file="FigS1a_draft.pdf",width=24/2.54,height=12/2.54)
plot(vec144,OCN_144,colPalette=hcl.colors(1000,'Rocket',rev=T), 
     addLegend=F,ylim=c(146300,160000),axes=FALSE)
points_colorscale(OCN_144$AG$X,
                  OCN_144$AG$Y,
                  medBioDF,
                  legend.lab="Y_DF")
dev.off()

# Fig. S7
nn=which(OCN_144$AG$Y>145000 & OCN_144$AG$nUpstream==1)
pp=vector("list",length(nn))
for (i in 1:length(nn)){
  pp[[i]] <- OCN_144$AG$downstreamPath[[nn[i]]][[66]]
}
names(pp) <- paste0('n',1:(length(nn)-1))
writeMat('utilities/zoomin144.mat', pp=pp) # for Fig. S1bc (in Matlab)

dA <- sort(c(unique(OCN_144$FD$A, OCN_480$FD$A)))
PA_144 <- PA_480 <- numeric(length(dA))
for (i in 1:length(PA_144)){
  PA_144[i] <- sum(OCN_144$FD$A>=dA[i])/length(OCN_144$FD$A)
  PA_480[i] <- sum(OCN_480$FD$A>=dA[i])/length(OCN_480$FD$A)
}

pdf(file="FigS7_draft.pdf",width=24/2.54,height=18/2.54)
plot(dA,PA_144,pch=19,col=rgb(0.66, 0.21, 0.66), log="xy")
points(dA,PA_480,pch=19,col=rgb(0.21, 0.59, 0.21))
dev.off()

}

