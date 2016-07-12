#=======================================================================================================
# program: function.R
# note: code submission for paper "Promoting similarity of model sparsity structures 
#                                  in integrative analysis of cancer genetic data".
#
# In this file: 
#    -  int.spBoost defines the function for the integrative sparse boosting.
#    -  AFT.getWeightedData applies the AFT weighting scheme as provided in the Appendix.
#=======================================================================================================
#########################################################################################################
#------------------------------------------------------------------------------------------------------//
# Function: int.spBoost                                                                                //
# Input:                                                                                               //
#           xList           as a list of length M that xList[[m]] is a sample size(n) by dim(p) matrix //
#           yList           as a list of length M that yList[[m]] is a vector                          //
#           lambda_df       as the switch for sparse boosting(=1) and boosting(=0)                     //
#           lambda_percent  as the tuning parameter for the similarity penalty                         //
#           n.boosting      as the number of boosting steps                                            //
#           increase        as the step size used in boosting                                          //
# Output:                                                                                              //
#    A list L of length two: estList and critList.                                                     //
#       - estList:  a list of length M that estList[[m]] is a n.boosting by p matrix                   //
#       - critList: a list of length 2 that both critList[["hdbic"]] and critList[["hdbic.perc"]]      //
#                   are M by n.boosting matrices                                                       //
#------------------------------------------------------------------------------------------------------//
#########################################################################################################

int.spBoost = function(xList,yList,lambda_percent,lambda_df=1,n.boosting=150,increase=0.1){
  
  # set the constants
  
  M             = length(xList)
  n.rv          = ncol(xList[[1]])
  n.sample      = sapply(yList,length)
  
  Identy_M_List = inxtx_vector_List = hatmatrix_List = Nothatmatrix_List = vector("list", M)
  
  sd_yList      = sapply(yList,sd)
  yList         = lapply(yList,function(x){ as.matrix((x-mean(x))/sd(x),ncol=1)}) 

  for (m in 1:M){  
    Identy_M_List[[m]]     = diag(1, n.sample[m]) 
    xList[[m]]             = apply(xList[[m]],2,function(x){(x-mean(x))}) 
    work.x.list            = split(xList[[m]], col(xList[[m]]))
    inxtx_vector_List[[m]] = sapply(work.x.list, function(x) {solve(t(x)%*%  x)}  )
    Nothatmatrix_List[[m]] = lapply(work.x.list, function(x) {solve(t(x)%*%  x)%*%t(x)}  )
    hatmatrix_List[[m]]    = lapply(1:n.rv, function(i)  work.x.list[[i]] %*% Nothatmatrix_List[[m]][[i]]) 
  }

  #  initialize containers in the main loop
  
  U_List  = F_List  =  act_sel_List = B_List = est_beta_List = vector("list", M)
  
  for (m in 1:M){
    U_List[[m]]          = yList[[m]]
    F_List[[m]]          = 0*yList[[m]]
    B_List[[m]]          = diag(0, n.sample[m])
    est_beta_List[[m]]   = matrix(0, nrow=(n.boosting+1), ncol= n.rv)
  }
  
  crit_List        = vector("list", 2)
  names(crit_List) = c("hdbic","hdbic.perc")
  
  for (i in 1:length(crit_List)){
    crit_List[[i]] = matrix(NA,nrow=M,ncol=n.boosting)
  }
  
  # (m,i.rv,1) = 1 if  i.rv is selected by the m set. 
  # (m,i.rv,2) = 1 if  i.rv is not selcted by the m set, but is selcted by other set.
  # (m.i.rv,3) = 1 if  i.rv is not selected by any of the sets.
  
  sel_tracker_array      = array(0,c(M,n.rv,3))  
  sel_tracker_array[,,3] = 1
  
  # main loop
  
  for (i in 1:n.boosting){                         
    for (m in 1:M){
      
      # get minimum sel
      
      xtB.use         = t(xList[[m]])%*% B_List[[m]]
      trace_vec.use   = sum(diag(B_List[[m]])) + 1 - inxtx_vector_List[[m]] * sapply(1:n.rv, function(i)  sum(xtB.use[i,]*xList[[m]][,i])) 
      RSS_vec.use     = sapply(1:n.rv, function(i.rv) sum((hatmatrix_List[[m]][[i.rv]]%*% U_List[[m]]-U_List[[m]])^2))
      
      if (i==1) { 
        percent.use   = 0  # datasets are not affected in the first loop.
      } else{
        percent.use   = percent * ((sel_tracker_array[m,,1]-rep(1,n.rv))==0) +
          (length(act)+1)/(3*length(unique(act))) *  ((sel_tracker_array[m,,2]-rep(1,n.rv))==0) +
          (length(act)+1)/(3*(length(unique(act))+1)) * ((sel_tracker_array[m,,3]-rep(1,n.rv))==0)
      }
      
      criterion.use = log(RSS_vec.use /n.sample[m]) + lambda_df* trace_vec.use *log(n.sample[m])/n.sample[m] - lambda_percent * percent.use
      sel.use       = which.min(criterion.use)
      
      sel_tracker_array[m ,sel.use,1]  = 1
      sel_tracker_array[m ,sel.use,2]  = 0
      sel_tracker_array[  ,sel.use,3]  = 0 
      sel_tracker_array[-m,sel.use,2]  = (sel_tracker_array[-m,sel.use,1]==0)*1
      
      # update percentage
      act_sel_List[[m]] = unique(c(sel.use, act_sel_List[[m]]))
      act               = unlist(act_sel_List)
      percent           = length(act)/(3*length(unique(act)))
      
      ## update beta
      
      beta.change = increase*Nothatmatrix_List[[m]][[sel.use]] %*% U_List[[m]]
      f.change    = increase*hatmatrix_List[[m]][[sel.use]]%*% U_List[[m]]
      
      est_beta_List[[m]][(i+1),]        = est_beta_List[[m]][i,]            
      est_beta_List[[m]][(i+1),sel.use] = est_beta_List[[m]][(i+1),sel.use] + beta.change
      
      F_List[[m]]     = F_List[[m]] + f.change
      U_List[[m]]     = yList[[m]]  - F_List[[m]]
      B_List[[m]]     = B_List[[m]] + increase *  hatmatrix_List[[m]][[sel.use]] %*% (Identy_M_List[[m]]- B_List[[m]])  
      
      part.df.m   = sum(diag(B_List[[m]]))*log(n.sample[m])/n.sample[m]
      avg.rss.m   = sum(U_List[[m]]^2)/n.sample[m]

      crit_List[["hdbic"]][m,i]      = log(avg.rss.m) + part.df.m*log(n.rv)
      crit_List[["hdbic.perc"]][m,i] = log(avg.rss.m) + part.df.m*log(n.rv) - lambda_percent * percent 
    }}
  
  for (m in 1:M){
    est_beta_List[[m]] = Matrix((est_beta_List[[m]][-1,]* sd_yList[m]),sparse = TRUE)
  }
  return(list(estList=est_beta_List,critList=crit_List))
}


## The following code applies AFT weighting scheme as provided in the Appendix.

AFT.getWeightedData <- function(x,obstime,delta.ifobserve){
  
  order.obstime             <- order(obstime)
  sort.x                    <- x[order.obstime,]
  sort.obstime              <- obstime[order.obstime]
  sort.delta                <- delta.ifobserve[order.obstime]
  
  n                         <- nrow(x)
  w.vec                     <- rep(0,n)
  w.vec[1]                  <- sort.delta[1]/n
  w.vec[-1]                 <- sort.delta[-1]*cumprod( (((n-1):1)/(n:2))^sort.delta[-n] )/((n-1):1)
  
  weightedmean.sort.obstime <- sum(sort.obstime*w.vec)/sum(w.vec)
  weightedmean.sort.x       <- apply(sort.x,2,function(y,w){sum(y*w)/sum(w)},w=w.vec)
  
  xstar                     <- sort.x*0
  ystar                     <- sort.obstime*0
  
  for (i in which(w.vec!=0)){
    xstar[i,] <- sqrt(w.vec[i])*(sort.x[i,]-weightedmean.sort.x)
    ystar[i]  <- sqrt(w.vec[i])*(sort.obstime[i]-weightedmean.sort.obstime)
  }
  
  return(list(x.use=xstar,y.use=ystar))
}
############################################# END of file  #################################################################

