rm(list=ls())
library(R.matlab)

# Niche model
# bodymass: a numeric vector of species body masses
# IF.Prop: a 0-1 value or "free": the proportion of nutrient feeders, which include basal species (producers) and consumers feeding also on nutrient. "free" means no constraint and take whatever the model generated.
# Connectance: a 0-1 value for the desired connectance of the generated web (not counting links to the nutrient block).
# Cannibalism: TRUE/FALSE, whether cannibalistic links are allowed in the generated web.
# Mutual.Feeding: TRUE/FALSE, whether bi-directional consumption is allowed in the generated web.
# Precise.Control: TRUE/FALSE, whether to force the generated connectance to be closer to the desired value by removing one layer of stochasticity.
# Beta.Scale: a scaling value that influence the beta distribution from with the foraging centre and radius are drawn. In general, a large value avoids super small radius.
# Checkpoint: TRUE/FALSE, whether to active a check point to disallow isolated nodes.

Nic = function(bodymass,
               IF.Prop = "free",  # when "free" usually it's the smallest one(s) being the NF.
               Connectance = 0.25,
               Cannibalism = TRUE,
               Mutual.Feeding = TRUE,
               Precise.Control = FALSE,
               Beta.Scale = 100,
               Checkpoint = FALSE,
               nD = 5){

  richness = length(bodymass)
  if(IF.Prop != "free") no.NF = round(richness * IF.Prop)
  # maybe we have inorganic feeders in general, then we partition them into nutrient and detritus feeders by flipping a coin
  # (larger more likely to eat detritus, smaller nutrients)
  # is it possible to have species that feed on both detritus and nutrients? maybe not, for simplicity
  # let's say 5% of nutrient feeders, 5% of detritus feeders. If the quotas are not attained, then we pick the smallest species and assign them to NF or DF
  # maybe also N and D have same size (just for simplicity). Then we flip a 50/50 coin
  # the 50/50 coin could be either N or D feeder; or one coin for N/D, then one for D given N and one for N given D (so that species feed on at least one type)
  # this gives 50% D+N; 25% D; 25% N
  # alternatively, a three sided dice (so we have less D+N)
  iter = 0
  unacceptable = FALSE
  repeat{ # The checkpoint repeat loop.
    inorganic.feeder = rep("N", richness) # a vector to later store identities of nutrient feeders.
    iter = iter + 1
    # Assign niche values
    if(Precise.Control){
      n = rank(bodymass)/richness
    } else {
      n = runif(richness)
      n = n[order(n)][rank(bodymass)] # regenerated, uniform distributed niche values according to bodymass ranking
    }

    # An empty diet matrix to be filled
    D = matrix(0, richness, richness)

    # Corrected connectance depending on the selected criteria
    if(Cannibalism){
      corrected.connectance = Connectance
    } else {
      corrected.connectance = (Connectance * richness^2)/(richness * (richness - 1))
    }

    # Assign diets
    alpha = (2 * corrected.connectance) * Beta.Scale
    beta = (1 - 2 * corrected.connectance) * Beta.Scale
    radius = n * rbeta(richness, alpha, beta)
    centre = unlist(lapply(c(1:richness), function(x){runif(1, radius[x]/2, n[x])})) # due to the setting here, none of the diets will cover n=0 or lower.
    for(i in 1:richness){
      diet = which(n <= centre[i] + radius[i]/2 & n >= centre[i] - radius[i]/2)
      D[diet, i] = 1
      if(length(setdiff(diet, i)) == 0) inorganic.feeder[i] = "Y"  # Those eat nothing (other than itself) are assigned to feed on nutrient.
      # These should be assigned to either nutrients or detritus (or both?)
      if(Cannibalism == F) diag(D) = 0 # just wire-up as no constraint then remove cannibalistic links, as in Cannibalism == F case connectances is already corrected.
    } # End of i for-loop.
    no.diet = which(inorganic.feeder == "Y") # those without feeding on any others, and already be assigned as NF.

    if(IF.Prop != "free"){
      # IF TOO MANY, DROP; IF TOO FEW, ADD THE LIGHTEST
      if(length(which(inorganic.feeder == "Y")) > no.NF){
        unacceptable = TRUE
      } else if(length(which(inorganic.feeder == "Y")) < no.NF){
        no.new.NF = no.NF - length(which(inorganic.feeder == "Y"))
        inorganic.feeder[which(inorganic.feeder == "N")[(length(which(inorganic.feeder == "N"))-no.new.NF+1):length(which(inorganic.feeder == "N"))]] <- "Y"
      }
    }
    if(iter > 100) { message("The given community and niche model settings are unlikely to generate a valid food web."); D = NA; break}
    if(Mutual.Feeding == FALSE) {if(any((D + t(D))[upper.tri(D)] > 1)) next} # if mutual feeding i.e., two-node loop detected, regenerate a new web.
    if (unacceptable) next
    if(Checkpoint == FALSE) {break} # If the checkpoint isn't activated, accept whatever web generated with the first iteration.
    if(Checkpoint == TRUE & prod(apply(D, 2, sum)[-no.diet]) > 0 & prod(apply(D, 1, sum)[no.diet]) > 0) {break} # If no isolated nodes, pass. Otherwise regenerate a new web.

  } # End of the checkpoint repeat loop.

  # now partition nutrient feeders into D, N, D+N
  # here I have a fixed minimum IF, but no control on D and N
  IF <- which(inorganic.feeder == "Y")
  DF <- sample(IF, nD)
  NF <- setdiff(IF, DF)

  return(list(D = D, DF = DF, NF = NF))
}

# Generate food web realizations
nSp <- 98 # makes 100 with nutrients & detritus
MassMean <- 10^-2
MassSD <-  10
a0 <- 1

nFW <- 30
bodymass_list <- dietMatrix_list <- detritusFeeders_list <- nutrientFeeders_list <- r_list <- A_list <- vector("list",nFW)
nms <- NULL
for (i in 1:nFW){
  nms <- c(nms, paste0('w',i))
}
names(bodymass_list) <- names(dietMatrix_list) <- names(detritusFeeders_list) <- names(nutrientFeeders_list) <- nms
names(r_list) <- names(A_list) <- nms

for (ind in 1:nFW){
  cat(sprintf('ind: %d \n',ind))
  set.seed(ind)
  nDfeed <- 2+ 3*floor((ind-1)/10) # first 10 have 2 Df; then 5 Df; then 8 Df
  flag <- NA; k <- 1
  while(isTRUE(is.na(flag))){
    bodymass <- rlnorm(nSp, log(MassMean), log(MassSD))
    bodymass <- sort(bodymass,decreasing = T)
    set.seed(ind+100*k)
    Web <-  Nic(bodymass, IF.Prop=0.1, Connectance=0.1, Cannibalism=F, Mutual.Feeding=F, nD=nDfeed, Checkpoint = T) # no mutual feeding, no cannibalism allows to calculate energy influxes from D, N
    flag <- Web$D
    k <- k+1
  }
  
  D <- Web$D # the generated food web as a diet matrix, cols eat rows.
  DF <- c(Web$DF, Web$DNF)
  NF <- c(Web$NF, Web$DNF)

  r_list[[ind]] <- -4.15e-8*bodymass^-0.25
  D_N <- cbind(rbind(D, numeric(nSp)),numeric(nSp+1))     # Diet matrix with detritus added
  D_N <- cbind(rbind(D_N, numeric(nSp+1)),numeric(nSp+2)) # Diet matrix with nutrients added
  D_N[nSp+1,DF] <- 1
  D_N[nSp+2,NF] <- 1
  bodymass_N <- c(bodymass,4e-5,2e-5) # detritus mass set to 4e-5, nutrient mass set to 2e-5
  bodymass_list[[ind]] <- bodymass
  dietMatrix_list[[ind]] <- D
  detritusFeeders_list[[ind]] <- DF
  nutrientFeeders_list[[ind]] <- NF

  A <- matrix(0,nSp+2,nSp+2)
  for (i in 1:(nSp+2)){
    for (j in 1:(nSp+2)){
      A[i,j] <- -2.72*D_N[i,j]*bodymass_N[i]^0.63*bodymass_N[j]^0.42 +
        0.5*2.72*D_N[j,i]*bodymass_N[i]^(-0.37)*bodymass_N[j]^1.42
    }
  }
  diag(A) <- diag(A) + c(-a0*bodymass^0.5, 0, 0)
  A_list[[ind]] <- A
}

writeMat('utilities/30FW_DF_2_5_8.mat',bodymass_list=bodymass_list, dietMatrix_list=dietMatrix_list,
         detritusFeeders_list=detritusFeeders_list,
         nutrientFeeders_list=nutrientFeeders_list, r_list=r_list, A_list=A_list)

