## Optimization problem to solve for partitioning parameters for soybean

multiyear_BioCro_optim_obj <- function(optim_params, biocro.fun, ExpData, num_rows, weights, wts2, RootVals){

  cost.avg <- 0

  for (i in 1:length(biocro.fun)) {

    gro.opt <- match.fun(biocro.fun[[i]])

    result <- gro.opt(optim_params)

#add a penalty on kLeaf, kStem, kShell & kGrain
    doy = result$doy
    kLeaf   = result$kLeaf
    kStem   = result$kStem
    kShell  = result$kShell
    kGrain  = result$kGrain
#get the first doys when these terms become non-zero
    doy_leaf  = doy[kLeaf>0.01][1]
    doy_stem  = doy[kStem>0.01][1]
    doy_shell = doy[kShell>0.01][1]
    doy_grain = doy[kGrain>0.01][1]
#add a penalty if 
#1. doy_leaf and doy_stem are separated by more than 5 days
#2. doy_shell or doy_grain occurs before doy_leaf
    if(is.na(doy_leaf)|is.na(doy_stem)|is.na(doy_shell)|is.na(doy_grain)){
       penalty = 9999
    }else if(abs(doy_leaf - doy_stem) > 5){
       penalty = 9999  #now i'm simply not allowing such case 
#    }else if(doy_shell < doy_leaf | doy_grain < doy_leaf){
#       penalty = 9999
#    }else if((doy_leaf-doy[1]) < 7 | (doy_stem-doy[1]) < 7){ #make sure leaf does not start very early
#       penalty = 9999
    }else{
       penalty = 0
    }

    if (nrow(result) < num_rows[i]) {
      # if simulation does not complete, assign a very high cost and exit
      cost.avg <- 1e10
      break

    } else{

      TrueValues <- ExpData[[i]]

      # Predicted values at DOY equal the DOYs in experimental data
      Pred <- data.frame("DOY"=TrueValues$DOY)
      doy_inds <- which(result$time %in% TrueValues$DOY)
      Pred$Stem <- result$Stem[doy_inds]
      Pred$Leaf <- result$Leaf[doy_inds]
      Pred$Shell  <- result$Shell[doy_inds]
      Pred$Seed <- result$Grain[doy_inds]

      doy_inds.Root <- which(result$time %in% RootVals[[i]]$DOY)
      Pred.Root <- result$Root[doy_inds.Root]


      # factor to scale experimental and simulated results between 0 and ~1 for all components
      scale.leaf <- max(TrueValues$Leaf)
      scale.stem <- max(TrueValues$Stem)
      scale.shell  <- max(TrueValues$Shell)
      scale.seed <- max(TrueValues$Seed)
      scale.Root <- max(RootVals[[i]]$Root)

      # weights
      wts <- weights[[i]]

      # weighted rmses
      err.stem <- sum(wts$Stem*(((Pred$Stem-TrueValues$Stem)/scale.stem)^2))/length(Pred$Stem)
      err.leaf <- sum(wts$Leaf*(((Pred$Leaf-TrueValues$Leaf)/scale.leaf)^2))/length(Pred$Leaf)
      err.shell <- sum(wts$Shell*(((Pred$Shell-TrueValues$Shell)/scale.shell)^2))/length(Pred$Shell)
      err.seed <- sum(wts$Seed*(((Pred$Seed-TrueValues$Seed)/scale.seed)^2))/length(Pred$Seed)
      err.Root <- sum(((Pred.Root-RootVals[[i]]$Root)/scale.Root)^2)/length(Pred.Root)

      cost <- (wts2$Stem*err.stem + wts2$Shell*err.shell + 
               wts2$Leaf*err.leaf + wts2$Root*err.Root +
               wts2$Seed*err.seed)

      cost <- round(100 * cost,2) / length(biocro.fun)
      cost.avg <- cost.avg + cost + penalty

    }

  }

  return(cost.avg)

}

