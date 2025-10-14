#######################################################
# Oct 14, 2025
# All required functions needed in parallel_master.r
######################################################
num_age = 9

update_outputdate = function(shortDate=FALSE) {
  ifelse(shortDate, substr(gsub("-", "", Sys.Date()),3,8), gsub("-", "", Sys.Date()))
}

fun_ageing <- function(y){
  y_new <- y
  num_y <- length(y)
  y_change <- y[1:(num_y-1)]/5
  y_new[1:(num_y-1)] <- y_new[1:(num_y-1)] - y_change
  y_new[2:num_y] <- y_new[2:num_y] + y_change
  return(y_new)
} # end fun_ageing

epimodel_fast <- function(initialstates, parameters, configuration, 
                          intervention, adjContactMatrix)
{
  P = parameters
  S = initialstates
  i = 0 # i relates to time step
  # Calculate number of deaths for both total population and advanced stages
  deaths = initialstates
  deaths <- lapply(deaths, function(x) {
    if (is.matrix(x) || is.data.frame(x)) {
      x[1, ] <- 0  # set first row to be 0
    }
    return(x)
  })
  ntreated = rep(0,configuration$nmonths)
  contactMat = P$contactMat
  alpha = P$alpha

  if(missing(alpha)){alpha = P$alpha}
  if (!missing(adjContactMatrix)){
    if (length(adjContactMatrix) == nrow(contactMat)){
      adjContactMatrix = exp(adjContactMatrix) # make sure all elements are non-negative
      contactMat = t(sapply(seq(length(adjContactMatrix)), function(x) adjContactMatrix[x]*contactMat[x,]))
    } 
    else{
      print( "length(adjContactMatrix)!=nrow(contactMat). Please check.")
      stop()
    }
  }
  
  # Some trigger parameters
  ageing = TRUE
  
  # Some pre-calculated results
  ## For D, J and X
  term_1 = diag(c(-1,P$kappa,1-P$kappa,rep(0,5)))
  term_2 = list('first'=diag(c(rep(1,7),0)),'second'=diag(P$pi * P$rho), 'third'=diag(P$rho))
  term_3 = diag(alpha)
  term_4 = matrix(c(0,P$kappa/P$iota,rep(0,6),
                    0,-P$kappa/P$iota,rep(0,6),
                    rep(0,2),-P$lambda_1,rep(0,5),
                    rep(0,2),P$lambda_1,-P$lambda_2,rep(0,4),
                    rep(0,3),P$lambda_2,-(P$lambda_3+P$lambda_4),rep(0,3),
                    rep(0,4),P$lambda_3,-(P$lambda_4+P$lambda_5+P$mu_4),rep(0,2),
                    rep(0,4),P$lambda_4,P$lambda_4,-(P$lambda_5+P$mu_5),rep(0,1),
                    rep(0,5),P$lambda_5,P$lambda_5,-P$mu_6),
                  nrow = 8,ncol = 8,byrow=T)
  term_5 = list('first'=diag(c(rep(1,5),rep(0,3))), 'second'=diag(P$mu))
  term_6 = matrix(c(P$beta,rep(0,71)),nrow=8,ncol=9)
  
  # For mortality
  term_7 = matrix(c(rep(P$mu,5),rep(c(P$mu_4,P$mu_5,P$mu_6),each=9)),nrow=8,ncol=9,byrow=T)
  
  ## Initial D, J and X
  current_D = rbind(S$Du[1,],S$D0[1,],S$D1[1,],S$D2[1,],
                    S$D3[1,],S$D4[1,],S$D5[1,],S$D6[1,])
  current_J = rbind(S$Ju[1,],S$J0[1,],S$J1[1,],S$J2[1,],
                    S$J3[1,],S$J4[1,],S$J5[1,],S$J6[1,])
  current_X = rbind(S$Xu[1,],S$X0[1,],S$X1[1,],S$X2[1,],
                    S$X3[1,],S$X4[1,],S$X5[1,],S$X6[1,])
  
  ## Inidial D, J and X for mortality
  death_D = rbind(deaths$Du[1,],deaths$D0[1,],deaths$D1[1,],deaths$D2[1,],
                  deaths$D3[1,],deaths$D4[1,],deaths$D5[1,],deaths$D6[1,])
  death_J = rbind(deaths$Ju[1,],deaths$J0[1,],deaths$J1[1,],deaths$J2[1,],
                  deaths$J3[1,],deaths$J4[1,],deaths$J5[1,],deaths$J6[1,])
  death_X = rbind(deaths$Xu[1,],deaths$X0[1,],deaths$X1[1,],deaths$X2[1,],
                  deaths$X3[1,],deaths$X4[1,],deaths$X5[1,],deaths$X6[1,])
  
  # Run the simulation for each day
  for(month in 1:configuration$nmonths){
    for(day in 1:configuration$nsteps_permonth){
      # Move on to next time point
      i = i+1
      
      Di = S$D0[i,]+S$D1[i,]+S$D2[i,]+S$D3[i,]+S$D4[i,]+S$D5[i,]+S$D6[i,] # Di = number of infected drug users by age group
      Ni = Di + S$Du[i,] # Ni = number of drug users by age groups
      S$total4[i+1,] = S$total4[i,]
      S$total5[i+1,] = S$total5[i,]
      S$total6[i+1,] = S$total6[i,]
      
      #1. Drug Users (new transmission)
      new_infected     = P$tau * S$Du[i,] * as.vector(colSums(P$delta * t(contactMat) * (Di/Ni)))
      new_infected_all = do.call("rbind", rep(list(new_infected), 8))
      
      rates_D = term_1 %*% new_infected_all + term_2$first %*% (current_J %*% term_2$second - current_D %*% term_3) + term_4 %*% current_D - term_5$first %*% current_D %*% term_5$second + term_6
      
      #2. Prisoners (no drugs)
      rates_J = term_2$first %*% (-current_J %*% term_2$third + current_D %*% term_3) + term_4 %*% current_J - term_5$first %*% current_J %*% term_5$second
      
      #3. Ex-Prisoners (no drugs)
      rates_X = term_2$first %*% current_J %*% (term_2$third - term_2$second) + term_4 %*% current_X - term_5$first %*% current_X %*% term_5$second
      
      # Update deaths
      death_D = death_D + current_D * term_7
      death_J = death_J + current_J * term_7
      death_X = death_X + current_X * term_7
      
      deaths$Du[i+1,] = death_D[1,];deaths$Ju[i+1,] = death_J[1,];deaths$Xu[i+1,] = death_X[1,]
      deaths$D0[i+1,] = death_D[2,];deaths$J0[i+1,] = death_J[2,];deaths$X0[i+1,] = death_X[2,]
      deaths$D1[i+1,] = death_D[3,];deaths$J1[i+1,] = death_J[3,];deaths$X1[i+1,] = death_X[3,]
      deaths$D2[i+1,] = death_D[4,];deaths$J2[i+1,] = death_J[4,];deaths$X2[i+1,] = death_X[4,]
      deaths$D3[i+1,] = death_D[5,];deaths$J3[i+1,] = death_J[5,];deaths$X3[i+1,] = death_X[5,]
      deaths$D4[i+1,] = death_D[6,];deaths$J4[i+1,] = death_J[6,];deaths$X4[i+1,] = death_X[6,]
      deaths$D5[i+1,] = death_D[7,];deaths$J5[i+1,] = death_J[7,];deaths$X5[i+1,] = death_X[7,]
      deaths$D6[i+1,] = death_D[8,];deaths$J6[i+1,] = death_J[8,];deaths$X6[i+1,] = death_X[8,]
      
      # Update current D, J and X
      current_D = current_D + rates_D
      current_J = current_J + rates_J
      current_X = current_X + rates_X
      
      S$Du[i+1,] = current_D[1,];S$Ju[i+1,] = current_J[1,];S$Xu[i+1,] = current_X[1,]
      S$D0[i+1,] = current_D[2,];S$J0[i+1,] = current_J[2,];S$X0[i+1,] = current_X[2,]
      S$D1[i+1,] = current_D[3,];S$J1[i+1,] = current_J[3,];S$X1[i+1,] = current_X[3,]
      S$D2[i+1,] = current_D[4,];S$J2[i+1,] = current_J[4,];S$X2[i+1,] = current_X[4,]
      S$D3[i+1,] = current_D[5,];S$J3[i+1,] = current_J[5,];S$X3[i+1,] = current_X[5,]
      S$D4[i+1,] = current_D[6,];S$J4[i+1,] = current_J[6,];S$X4[i+1,] = current_X[6,]
      S$D5[i+1,] = current_D[7,];S$J5[i+1,] = current_J[7,];S$X5[i+1,] = current_X[7,]
      S$D6[i+1,] = current_D[8,];S$J6[i+1,] = current_J[8,];S$X6[i+1,] = current_X[8,]
      
      S$total4[i+1,] = S$total4[i+1,] + P$lambda_3 * (S$X3[i,] + S$J3[i,] + S$D3[i,]) # sic: cumulate
      S$total5[i+1,] = S$total5[i+1,] + P$lambda_4 * (S$X4[i,] + S$J4[i,] + S$D4[i,]) 
      S$total6[i+1,] = S$total6[i+1,] + P$lambda_5 * (S$X5[i,] + S$J5[i,] + S$D5[i,]) 
      
      # Update time status
      S$t[i+1] = S$t[i] + 1
      deaths$t[i+1] = deaths$t[i] + 1
    } # end day for-loop
    
    # process for every month
    # treatment schedule - treat those in stages 1-3
    # percentages within groups:
    ntreated[month] = 0
    
    tx = rep(0,9)
    cures = rep(0,9)
    temp = S$D1[i+1,]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D1[i+1,] = S$D1[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$D2[i+1,]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D2[i+1,] = S$D2[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$D3[i+1,]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_2; S$D3[i+1,] = S$D3[i+1,] - temp + (temp*(1-parameters$theta_2))
    S$Du[i+1,] = S$Du[i+1,] + cures
    
    cures = rep(0,9)
    temp = S$J1[i+1,]*intervention$pJ1; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J1[i+1,] = S$J1[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$J2[i+1,]*intervention$pJ2; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J2[i+1,] = S$J2[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$J3[i+1,]*intervention$pJ3; tx = tx + temp; cures = cures + temp*parameters$theta_2; S$J3[i+1,] = S$J3[i+1,] - temp + (temp*(1-parameters$theta_2))
    S$Ju[i+1,] = S$Ju[i+1,] + cures
    
    cures = rep(0,9)
    temp = S$X1[i+1,]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X1[i+1,] = S$X1[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$X2[i+1,]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X2[i+1,] = S$X2[i+1,] - temp + (temp*(1-parameters$theta_1))
    temp = S$X3[i+1,]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_2; S$X3[i+1,] = S$X3[i+1,] - temp + (temp*(1-parameters$theta_2))
    S$Xu[i+1,] = S$Xu[i+1,] + cures 
    
    # Update current_D, J, X if there is intervention
    current_D[1,] = S$Du[i+1,]; current_D[3,] = S$D1[i+1,]; current_D[4,] = S$D2[i+1,]; current_D[5,] = S$D3[i+1,]
    current_J[1,] = S$Ju[i+1,]; current_J[3,] = S$J1[i+1,]; current_J[4,] = S$J2[i+1,]; current_J[5,] = S$J3[i+1,]
    current_X[1,] = S$Xu[i+1,]; current_X[3,] = S$X1[i+1,]; current_X[4,] = S$X2[i+1,]; current_X[5,] = S$X3[i+1,]
    
    ntreated[month] = ntreated[month] + sum(tx)
    
    if (ageing){
      # ageing process
      if ((month%%12)==0){
        current_D[1,] = S$Du[i+1,] = fun_ageing(S$Du[i+1,]);current_J[1,] = S$Ju[i+1,] = fun_ageing(S$Ju[i+1,]);current_X[1,] = S$Xu[i+1,] = fun_ageing(S$Xu[i+1,])
        current_D[2,] = S$D0[i+1,] = fun_ageing(S$D0[i+1,]);current_J[2,] = S$J0[i+1,] = fun_ageing(S$J0[i+1,]);current_X[2,] = S$X0[i+1,] = fun_ageing(S$X0[i+1,])
        current_D[3,] = S$D1[i+1,] = fun_ageing(S$D1[i+1,]);current_J[3,] = S$J1[i+1,] = fun_ageing(S$J1[i+1,]);current_X[3,] = S$X1[i+1,] = fun_ageing(S$X1[i+1,])
        current_D[4,] = S$D2[i+1,] = fun_ageing(S$D2[i+1,]);current_J[4,] = S$J2[i+1,] = fun_ageing(S$J2[i+1,]);current_X[4,] = S$X2[i+1,] = fun_ageing(S$X2[i+1,])
        current_D[5,] = S$D3[i+1,] = fun_ageing(S$D3[i+1,]);current_J[5,] = S$J3[i+1,] = fun_ageing(S$J3[i+1,]);current_X[5,] = S$X3[i+1,] = fun_ageing(S$X3[i+1,])
        current_D[6,] = S$D4[i+1,] = fun_ageing(S$D4[i+1,]);current_J[6,] = S$J4[i+1,] = fun_ageing(S$J4[i+1,]);current_X[6,] = S$X4[i+1,] = fun_ageing(S$X4[i+1,])
        current_D[7,] = S$D5[i+1,] = fun_ageing(S$D5[i+1,]);current_J[7,] = S$J5[i+1,] = fun_ageing(S$J5[i+1,]);current_X[7,] = S$X5[i+1,] = fun_ageing(S$X5[i+1,])
        current_D[8,] = S$D6[i+1,] = fun_ageing(S$D6[i+1,]);current_J[8,] = S$J6[i+1,] = fun_ageing(S$J6[i+1,]);current_X[8,] = S$X6[i+1,] = fun_ageing(S$X6[i+1,])
      } # end ageing process
    } # if false
    
  } # end month for-loop
  
  S$D=list(S$Du,S$D0,S$D1,S$D2,S$D3,S$D4,S$D5,S$D6) # list of matrix of (ntime, num_age)
  S$J=list(S$Ju,S$J0,S$J1,S$J2,S$J3,S$J4,S$J5,S$J6)
  S$X=list(S$Xu,S$X0,S$X1,S$X2,S$X3,S$X4,S$X5,S$X6)
  
  deaths$D = list(deaths$Du,deaths$D0,deaths$D1,deaths$D2,deaths$D3,deaths$D4,deaths$D5,deaths$D6)
  deaths$J = list(deaths$Ju,deaths$J0,deaths$J1,deaths$J2,deaths$J3,deaths$J4,deaths$J5,deaths$J6)
  deaths$X = list(deaths$Xu,deaths$X0,deaths$X1,deaths$X2,deaths$X3,deaths$X4,deaths$X5,deaths$X6)
  
  return(list(states = S, treats = ntreated, deaths = deaths))
}

create_configuration = function(ny=50, nsteps_permonth=30)
{
  configuration = within(list(), {
    nyears = ny
    nmonths = nyears*12
    nsteps_permonth = nsteps_permonth
    nsteps = nmonths*nsteps_permonth
  })
  return(configuration)
}


set_parameters = function(nsteps_permonth=30){
  nsteps_peryr = 12 * nsteps_permonth 
  parameters=list(
    #======= parameters from prison data ===========
    alpha = c(0.067,
              0.21,
              0.21,
              0.15,
              0.20,
              0.31,
              0.14,
              0.4,
              0.025)/nsteps_peryr, #0.5374/360,   # arrest rate per current IDU per year
    rho   = c(0.5472816,
              0.6104789,
              0.5816654,
              0.5151572,
              0.4358284,
              0.4178032,
              0.4031776,
              0.4029714,
              0.4315638)/nsteps_peryr,   # release rate per jailed former user per year
    delta = 0.1,
    beta  = 1861/nsteps_peryr  ,   # initiation (or "birth") rate per year
    tau   = 0.0057, #0.0017, #(with the assumption that all IDU events were actually associated with sharing) #0.0057(per-event probability of HCV infection after IDU following a sharing event), #0.0063,  #0.0000492/360  
    contactMat = rbind(c(7, 4, 1, 1, 0, 1, 1, 0, 0),#4.74
                       c(11, 34, 21, 11, 6, 2, 1, 1, 0.1*1),#2.14
                       c(7, 30, 80, 62, 30, 10, 2, 3, 0.1*3), #0.86
                       c(2, 10, 60, 121, 65, 38, 15, 4, 0.1*4), #0.47
                       c(1, 11, 22, 67, 107, 41, 18, 5, 0.1*5), #0.51
                       c(0, 4, 6, 22, 32, 31, 10, 4, 0.1*4),#1.34
                       c(0, 1, 1, 8, 10, 15, 11, 1, 0.1*1), #2.37
                       c(0, 10, 10, 11, 13, 13, 7, 6, 0.1*6),#1.36
                       c(0, 10, 10, 11, 13, 13, 7, 6, 0.1*6))/(3 * nsteps_permonth), #dividing by 90 to convert from 3-month rate to daily rate; # matrix(c(...), ncol=9, byrow=TRUE) or directly use rbind()
    # age-specific; updated
    pi=c(0.6561922, 0.6596096, 0.645692, 0.6578092, 0.6687155, 0.6996991, 0.6669786, 0.6351017, 0.5369245), #recidivism probability
    #mu=c(0.001267,0.0003,0.0003,0.0004,0.0005,0.0007,0.0014,0.0023,0.0356)/nsteps_peryr,
    mu=c(0.001267,0.0003,0.0003,0.0004,0.0005,0.0007,0.0014,0.0023,0.0161)/nsteps_peryr,
    
    #======= additional parameters ======
    kappa = 0.26, # proportion of infections leading to spontaneous clearance (Uniform (0.22, 0.29)) (NICE)
    iota  = 0.5*nsteps_peryr, # duration in which spontaneous clearers are infectious
    
    # transition probabilities per year (UK NICE)
    lambda_1 = 0.025/nsteps_peryr, #rbeta(1, 38.0859, 1485.3516)/360   mild to moderate
    lambda_2 = 0.037/nsteps_peryr, #rbeta(1, 26.905,700.2582)/360      
    lambda_3 = 0.039/nsteps_peryr, #rbeta(1, 14.6168,360.1732)/360     
    lambda_4 = 0.014/nsteps_peryr, #rbeta(1, 1.9326,136.1074)/360      
    lambda_5 = 0.03/nsteps_peryr , #rbeta(1, 6.5256,210.9945)/360      HCC to LT
    
    mu_4 = 0.13/nsteps_peryr,      #rbeta(1, 147.03,983.97)/360   decomp to death
    mu_5 = 0.43/nsteps_peryr,      #rbeta(1, 117.1033,155.23)/360 HCC to death
    mu_6 = 0.21/nsteps_peryr,      #rbeta(1, 16.2762,61.2294)/360 LT to death
    
    theta_1 = 0.97,       #cure rate (SVR) for genotype 1, mild/moderate for Harvoni  ##0.45 (NICE) runif(1,0.4,0.5)                         
    theta_2 = 0.94        #cure rate (SVR) for genotype 1,        ##0.25 (NICE)   
    
  )
  return(parameters)
}

make_initial_states=function(configuration,equilibrium=FALSE,parameters)
{
  N=configuration$nsteps + 1
  output = list()
  output$Du = output$D0 = output$D1 = output$D2 = output$D3 = output$D4 = output$D5 = output$D6 = matrix(0, ncol=num_age, nrow=N) # rep(0,N)
  output$Ju = output$J0 = output$J1 = output$J2 = output$J3 = output$J4 = output$J5 = output$J6 = matrix(0, ncol=num_age, nrow=N) # rep(0,N)
  output$Xu = output$X0 = output$X1 = output$X2 = output$X3 = output$X4 = output$X5 = output$X6 = matrix(0, ncol=num_age, nrow=N) # rep(0,N)
  output$t = rep(0,N)
  output$total6 = output$total5 = output$total4 = matrix(0, ncol=num_age, nrow=N) # rep(0,N)
  
  if(!equilibrium)
  {
    output$Du[1,] = c(504,1304,1292,938,772,686,518,498,466) #c(2292,8500,11762,11432,11470,12632,11826,9780,9010) #2*number of people in J assuming arrest prob of 50% #c(400,700,1100,1000,1000,1200,1200,1000,1200) # output$Du[1] = 8464
    output$D0[1,] = c(110, 290, 366, 328, 424, 598, 444, 380, 266)#c(802.2, 2975.0, 4116.7, 4001.2, 4014.5, 4421.2, 4139.1, 3423.0, 3153.5) #assuming prev of 35% #c(rep(100,9)) #c(10,60,100,100,100,200,150,140,140) #1000 
    output$Ju[1,] = c(252, 652, 646, 469, 386, 343, 259, 249, 233)# based on prison data
    output$J0[1,] = c(55,145,183,164,212,299,222,190,133) # based on prison data
    output$Xu[1,] = c(900,7200,12000,11000,11000,13000,10000,10000,10000) #90000#136943
  }
  
  if(equilibrium)#haven't done this yet
  {
    load('output/equilibrium_240807.rdata')
    output$Du[1,] = equilibrium$Du; output$Ju[1,] = equilibrium$Ju; output$Xu[1,] = equilibrium$Xu
    output$D0[1,] = equilibrium$D0; output$J0[1,] = equilibrium$J0; output$X0[1,] = equilibrium$X0
    output$D1[1,] = equilibrium$D1; output$J1[1,] = equilibrium$J1; output$X1[1,] = equilibrium$X1
    output$D2[1,] = equilibrium$D2; output$J2[1,] = equilibrium$J2; output$X2[1,] = equilibrium$X2
    output$D3[1,] = equilibrium$D3; output$J3[1,] = equilibrium$J3; output$X3[1,] = equilibrium$X3
    output$D4[1,] = equilibrium$D4; output$J4[1,] = equilibrium$J4; output$X4[1,] = equilibrium$X4
    output$D5[1,] = equilibrium$D5; output$J5[1,] = equilibrium$J5; output$X5[1,] = equilibrium$X5
    output$D6[1,] = equilibrium$D6; output$J6[1,] = equilibrium$J6; output$X6[1,] = equilibrium$X6
  }
  output$total4[1,] = output$D4[1,] + output$J4[1,] + output$X4[1,]
  output$total5[1,] = output$D5[1,] + output$J5[1,] + output$X5[1,]
  output$total6[1,] = output$D6[1,] + output$J6[1,] + output$X6[1,]
  return(output)
}

final_prevalence = function(X) # among jailed people
{
  i = nrow(X$states$D6); #as a matrix # length(X$states$D6) as a vector
  num_agegp = ncol(X$states$D6);
  # single prevalence (non age-specific)
  Ji_sumage <- sapply(X$states$J,function(x) sum(x[i,])) # sum each J at time i across all age groups; output is a vector of Ju, J0,J1,...,J6 at time i
  infected = sum(Ji_sumage[-1])
  total = sum(Ji_sumage)
  prev = infected/total
  prevalence_byage <- rep(0, num_agegp)
  for (iage in 1:num_agegp){
    Ji_sumage_age <- sapply(X$states$J,function(x) x[i,iage]) # sum each J at time i across all age groups; output is a vector of Ju, J0,J1,...,J6 at time i
    infected_age = sum(Ji_sumage_age[-1])
    total_age = sum(Ji_sumage_age)
    prevalence_byage[iage] <- infected_age/total_age}
  return(list(prevalence_overall=prev, prevalence_byage=prevalence_byage))}


store_states=function(X, outfile='output/equilibrium_240807.rdata')
{
  i = nrow(X$states$Du) # length(X$states$Du)
  equilibrium = list(
    Du = X$states$Du[i,], Ju = X$states$Ju[i,], Xu = X$states$Xu[i,],
    D0 = X$states$D0[i,], J0 = X$states$J0[i,], X0 = X$states$X0[i,],
    D1 = X$states$D1[i,], J1 = X$states$J1[i,], X1 = X$states$X1[i,],
    D2 = X$states$D2[i,], J2 = X$states$J2[i,], X2 = X$states$X2[i,],
    D3 = X$states$D3[i,], J3 = X$states$J3[i,], X3 = X$states$X3[i,],
    D4 = X$states$D4[i,], J4 = X$states$J4[i,], X4 = X$states$X4[i,],
    D5 = X$states$D5[i,], J5 = X$states$J5[i,], X5 = X$states$X5[i,],
    D6 = X$states$D6[i,], J6 = X$states$J6[i,], X6 = X$states$X6[i,]
  )
  save(equilibrium,file=outfile)
}


tabulator = function(X,intervention,file='output/table.rds',header=FALSE,resolution='LOW')
{
  if(header)
  {
    cat('Strategy,t,Du1,Du2,Du3,Du4,Du5,Du6,Du7,Du8,Du9,D01,D02,D03,D04,D05,D06,D07,D08,D09,D11,D12,D13,D14,D15,D16,D17,D18,D19,D21,D22,D23,D24,D25,D26,D27,D28,D29,D31,D32,D33,D34,D35,D36,D37,D38,D39,D41,D42,D43,D44,D45,D46,D47,D48,D49,D51,D52,D53,D54,D55,D56,D57,D58,D59,D61,D62,D63,D64,D65,D66,D67,D68,D69,Ju1,Ju2,Ju3,Ju4,Ju5,Ju6,Ju7,Ju8,Ju9,J01,J02,J03,J04,J05,J06,J07,J08,J09,J11,J12,J13,J14,J15,J16,J17,J18,J19,J21,J22,J23,J24,J25,J26,J27,J28,J29,J31,J32,J33,J34,J35,J36,J37,J38,J39,J41,J42,J43,J44,J45,J46,J47,J48,J49,J51,J52,J53,J54,J55,J56,J57,J58,J59,J61,J62,J63,J64,J65,J66,J67,J68,J69,Xu1,Xu2,Xu3,Xu4,Xu5,Xu6,Xu7,Xu8,Xu9,X01,X02,X03,X04,X05,X06,X07,X08,X09,X11,X12,X13,X14,X15,X16,X17,X18,X19,X21,X22,X23,X24,X25,X26,X27,X28,X29,X31,X32,X33,X34,X35,X36,X37,X38,X39,X41,X42,X43,X44,X45,X46,X47,X48,X49,X51,X52,X53,X54,X55,X56,X57,X58,X59,X61,X62,X63,X64,X65,X66,X67,X68,X69,ctx,C41,C42,C43,C44,C45,C46,C47,C48,C49,C51,C52,C53,C54,C55,C56,C57,C58,C59,C61,C62,C63,C64,C65,C66,C67,C68,C69\n',file=file)
    return(NULL)
  }
  if(!header)
  {
    if(resolution=='HIGH')i = seq(1,length(X$states$t),30)
    if(resolution=='HIGH')j = seq(1,length(X$treats),1)
    if(resolution=='LOW')i = seq(1,length(X$states$t),360)
    if(resolution=='LOW')j = seq(12,length(X$treats),12)
    #output = cbind(X$states$t[i],X$states$D[i,],X$states$J[i,],X$states$X[i,])
    combD <- as.data.frame(do.call(cbind, X$states$D))
    combJ <- as.data.frame(do.call(cbind, X$states$J))
    combX <- as.data.frame(do.call(cbind, X$states$X))
    output = cbind(X$states$t[i],combD[i,],combJ[i,],combX[i,])
    temp1 = rep(intervention$strategy,length(i))
    temp2 = c(0,cumsum(X$treats)[j])
    temp3 = as.data.frame(cbind(X$states$total4,X$states$total5,X$states$total6))[i,]
    
    # For deaths
    combD_d = as.data.frame(do.call(cbind, X$deaths$D))
    combJ_d = as.data.frame(do.call(cbind, X$deaths$J))
    combX_d = as.data.frame(do.call(cbind, X$deaths$X))
    output_d = cbind(X$states$t[i],combD_d[i,],combJ_d[i,],combX_d[i,])
    
    output_f = cbind(temp1,output,temp2,temp3,output_d)
    write.table(output_f,file = file,sep=',',col.names = FALSE,row.names = FALSE,append = TRUE)
  }
}

rdstabulator = function(X,intervention)
{
  i = seq(1,length(X$states$t),30)
  j = seq(1,length(X$treats),1)
  #output = cbind(X$states$t[i],X$states$D[i,],X$states$J[i,],X$states$X[i,])
  combD <- as.data.frame(do.call(cbind, X$states$D))
  combJ <- as.data.frame(do.call(cbind, X$states$J))
  combX <- as.data.frame(do.call(cbind, X$states$X))
  output = cbind(X$states$t[i],combD[i,],combJ[i,],combX[i,])
  temp1 = rep(intervention$strategy,length(i))
  temp2 = c(0,cumsum(X$treats)[j])
  temp3 = as.data.frame(cbind(X$states$total4,X$states$total5,X$states$total6))[i,]
  output_f = cbind(temp1,output,temp2,temp3)
  return(output_f)
}

#================================= MCMC ======================================#
fun_final_prev = function(xparam){
  X = epimodel_fast(initialstates,parameters,configuration,intervention, 
                    adjContactMat=xparam)
  return( final_prevalence(X) )
}

fun_loglikelihood = function(xparam, obs_data = obs_data){
  # assign parameters (xparam)
  X = epimodel_fast(initialstates,parameters,configuration,intervention, 
                    adjContactMat=xparam)
  final_prev = final_prevalence(X)
  final_prev_byage = final_prev$prevalence_byage
  
  logLL = with(obs_data, pos *log(final_prev_byage) + (tot-pos) *log(1-final_prev_byage))	#binomial 
  return(sum(logLL))
}  # end fun_loglikelihood

fun_LBUB = function(x){ # to indicate the boundary
  return( all(LB<=x & x<=UB) ) # single &/| for comparing all elements in the vector
} # fun_LBUB

fun_logLL_MCMC = function(x){
  # set LB / UB within the function
  if (!fun_LBUB(x)){
    return( -1e7 )
  }
  return( fun_loglikelihood(xparam=x, obs_data) )
} # fun_logLL_MCMC

fun_plotprev = function(prev, obs_data, arrows_length=0.05){
  num_agegp = length(prev)
  plot(NA, xlim=c(1,num_agegp)+c(-0.25,0.25), ylim=c(0,1), xlab="Age group",ylab="Age-specific prevalence at the final step (%)", axes=FALSE)
  axis(1, at=1:num_agegp, labels=labs)
  axis(2, at=seq(0,1,by=0.2), labels=100*seq(0,1,by=0.2), las=1)
  obs_prev = with(obs_data, 
                  sapply(seq(num_agegp),
                         function(ii) c(obsprev=pos[ii]/tot[ii], lowlim=binom.test(pos[ii],tot[ii])$conf.int[1], upplim=binom.test(pos[ii],tot[ii])$conf.int[2])))
  arrows(x0=seq(num_agegp), y0=obs_prev["lowlim",],y1=obs_prev["upplim",], code=3,angle=90, length=arrows_length)
  points(x=seq(num_agegp), y=obs_prev["obsprev",], pch=1)
  points(x=seq(num_agegp), y=prev, pch=4)
  legend("topright", legend=c("Fitted model","Empirical data"), pch=c(4,1), bty="n")
} # end fun_plotprev


