strategies = read.csv("strategies_new.csv");
strategies$nD = strategies$nD/12;   # Original is annual number so convert it to monthly
strategies$nJ = strategies$nJ/12;   # "
strategies$nX = strategies$nX/12;   # "
strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
strategies$pJ1 = strategies$pJ1/12; # "
strategies$pJ2 = strategies$pJ2/12; # "
strategies$pJ3 = strategies$pJ3/12; # "
strategies$pX = strategies$pX/12;   # "



epimodel = function(initialstates,parameters,configuration,intervention)
{
  P=parameters
  S=initialstates
  i=0
  ntreated = rep(0,configuration$nmonths)
  for(month in 1:configuration$nmonths)
  {
    for(day in 1:30)
    {
      i = i+1
      S$total4[i+1] = S$total4[i]
      S$total5[i+1] = S$total5[i]
      S$total6[i+1] = S$total6[i]
      #Drug Users (in the community)
      Di = (S$D0[i]+S$D1[i]+S$D2[i]+S$D3[i]+S$D4[i]+S$D5[i]+S$D6[i])
      rates_Du = 0 - P$tau*S$Du[i]*Di         + P$pi*P$rho*S$Ju[i]  - P$alpha*S$Du[i] + P$kappa*S$D0[i]/P$iota    - P$mu*S$Du[i]    + P$beta           #Uninfected
      rates_D0 = P$kappa*P$tau*S$Du[i]*Di     + P$pi*P$rho*S$J0[i]  - P$alpha*S$D0[i] - P$kappa*S$D0[i]/P$iota    - P$mu*S$D0[i]                     #acute infection
      rates_D1 = (1-P$kappa)*P$tau*S$Du[i]*Di + P$pi*P$rho*S$J1[i]  - P$alpha*S$D1[i] - P$lambda_1*S$D1[i] - P$mu*S$D1[i]                     #mild
      rates_D2 = P$lambda_1*S$D1[i]           + P$pi*P$rho*S$J2[i]  - P$alpha*S$D2[i] - P$lambda_2*S$D2[i] - P$mu*S$D2[i]                     #moderate
      rates_D3 = P$lambda_2*S$D2[i]           + P$pi*P$rho*S$J3[i]  - P$alpha*S$D3[i] - P$lambda_3*S$D3[i] - P$mu*S$D3[i]    - P$lambda_4*S$D3[i] #comp cirrhosis
      rates_D4 = P$lambda_3*S$D3[i]           + P$pi*P$rho*S$J4[i]  - P$alpha*S$D4[i] - P$lambda_4*S$D4[i] - P$mu_4*S$D4[i]  - P$lambda_5*S$D4[i] #decomp cirrhosis
      rates_D5 = P$lambda_4*S$D4[i]           + P$pi*P$rho*S$J5[i]  - P$alpha*S$D5[i] - P$lambda_5*S$D5[i] - P$mu_5*S$D5[i]  + P$lambda_4*S$D3[i] #HCC
      rates_D6 = P$lambda_5*S$D5[i]                                                                        - P$mu_6*S$D6[i]  + P$lambda_5*S$D4[i] #LT - no further transitions to J and X for LT stage
      
      S$Du[i+1] <- S$Du[i] + rates_Du
      S$D0[i+1] <- S$D0[i] + rates_D0
      S$D1[i+1] <- S$D1[i] + rates_D1
      S$D2[i+1] <- S$D2[i] + rates_D2
      S$D3[i+1] <- S$D3[i] + rates_D3
      S$D4[i+1] <- S$D4[i] + rates_D4
      S$D5[i+1] <- S$D5[i] + rates_D5
      S$D6[i+1] <- S$D6[i] + rates_D6
      S$total4[i+1] <- S$total4[i+1] + P$lambda_3*S$D3[i] # sic: cumulate
      S$total5[i+1] <- S$total5[i+1] + P$lambda_4*S$D4[i] 
      S$total6[i+1] <- S$total6[i+1] + P$lambda_5*S$D5[i] 
      
      #IN Prison: drugs X
      rates_Ju = P$alpha*S$Du[i] - P$pi*P$rho*S$Ju[i] - (1-P$pi)*P$rho*S$Ju[i] + P$kappa*S$J0[i]/P$iota                    - P$mu*S$Ju[i]                    #Uninfected
      rates_J0 = P$alpha*S$D0[i] - P$pi*P$rho*S$J0[i] - (1-P$pi)*P$rho*S$J0[i] - P$kappa*S$J0[i]/P$iota                    - P$mu*S$J0[i]                    #acute infection
      rates_J1 = P$alpha*S$D1[i] - P$pi*P$rho*S$J1[i] - (1-P$pi)*P$rho*S$J1[i]                      - P$lambda_1*S$J1[i] - P$mu*S$J1[i]                    #mild
      rates_J2 = P$alpha*S$D2[i] - P$pi*P$rho*S$J2[i] - (1-P$pi)*P$rho*S$J2[i] + P$lambda_1*S$J1[i] - P$lambda_2*S$J2[i] - P$mu*S$J2[i]                    #moderate
      rates_J3 = P$alpha*S$D3[i] - P$pi*P$rho*S$J3[i] - (1-P$pi)*P$rho*S$J3[i] + P$lambda_2*S$J2[i] - P$lambda_3*S$J3[i] - P$mu*S$J3[i]   - P$lambda_4*S$J3[i] #comp cirrhosis
      rates_J4 = P$alpha*S$D4[i] - P$pi*P$rho*S$J4[i] - (1-P$pi)*P$rho*S$J4[i] + P$lambda_3*S$J3[i] - P$lambda_4*S$J4[i] - P$mu_4*S$J4[i] - P$lambda_5*S$J4[i] #decomp cirrhosis
      rates_J5 = P$alpha*S$D5[i] - P$pi*P$rho*S$J5[i] - (1-P$pi)*P$rho*S$J5[i] + P$lambda_4*S$J4[i] - P$lambda_5*S$J5[i] - P$mu_5*S$J5[i] + P$lambda_4*S$J3[i] #HCC
      rates_J6 =                                                                 P$lambda_5*S$J5[i]                      - P$mu_6*S$J6[i] + P$lambda_5*S$J4[i] #LT - no transitions from D6 (assume they don't drugs and get arrested at this stage)
      S$Ju[i+1] <- S$Ju[i] + rates_Ju
      S$J0[i+1] <- S$J0[i] + rates_J0
      S$J1[i+1] <- S$J1[i] + rates_J1
      S$J2[i+1] <- S$J2[i] + rates_J2
      S$J3[i+1] <- S$J3[i] + rates_J3
      S$J4[i+1] <- S$J4[i] + rates_J4
      S$J5[i+1] <- S$J5[i] + rates_J5
      S$J6[i+1] <- S$J6[i] + rates_J6
      S$total4[i+1] <- S$total4[i+1] + P$lambda_3*S$J3[i] # sic: cumulate
      S$total5[i+1] <- S$total5[i+1] + P$lambda_4*S$J4[i] 
      S$total6[i+1] <- S$total6[i+1] + P$lambda_5*S$J5[i] 
      
      #OUT of prison: drugs X
      rates_Xu = (1-P$pi)*P$rho*S$Ju[i] + P$kappa*S$X0[i]/P$iota                     - P$mu*S$Xu[i]                    #Uninfected
      rates_X0 = (1-P$pi)*P$rho*S$J0[i] - P$kappa*S$X0[i]/P$iota                     - P$mu*S$X0[i]                    #acute infection
      rates_X1 = (1-P$pi)*P$rho*S$J1[i]                      - P$lambda_1*S$X1[i] - P$mu*S$X1[i]                    #mild
      rates_X2 = (1-P$pi)*P$rho*S$J2[i] + P$lambda_1*S$X1[i] - P$lambda_2*S$X2[i] - P$mu*S$X2[i]                    #moderate
      rates_X3 = (1-P$pi)*P$rho*S$J3[i] + P$lambda_2*S$X2[i] - P$lambda_3*S$X3[i] - P$mu*S$X3[i]   - P$lambda_4*S$X3[i] #comp cirrhosis
      rates_X4 = (1-P$pi)*P$rho*S$J4[i] + P$lambda_3*S$X3[i] - P$lambda_4*S$X4[i] - P$mu_4*S$X4[i] - P$lambda_5*S$X4[i] #decomp cirrhosis
      rates_X5 = (1-P$pi)*P$rho*S$J5[i] + P$lambda_4*S$X4[i] - P$lambda_5*S$X5[i] - P$mu_5*S$X5[i] + P$lambda_4*S$X3[i] #HCC
      rates_X6 =                          P$lambda_5*S$X5[i]                      - P$mu_6*S$X6[i] + P$lambda_5*S$X4[i] #LT - no transitions from D6 (assume they don't drugs and get arrested at this stage)
      S$Xu[i+1] <- S$Xu[i] + rates_Xu
      S$X0[i+1] <- S$X0[i] + rates_X0
      S$X1[i+1] <- S$X1[i] + rates_X1
      S$X2[i+1] <- S$X2[i] + rates_X2
      S$X3[i+1] <- S$X3[i] + rates_X3
      S$X4[i+1] <- S$X4[i] + rates_X4
      S$X5[i+1] <- S$X5[i] + rates_X5
      S$X6[i+1] <- S$X6[i] + rates_X6
      S$total4[i+1] <- S$total4[i+1] + P$lambda_3*S$X3[i] # sic: cumulate
      S$total5[i+1] <- S$total5[i+1] + P$lambda_4*S$X4[i] 
      S$total6[i+1] <- S$total6[i+1] + P$lambda_5*S$X5[i] 
      
      
      S$t[i+1] <- S$t[i] + 1
    }
    # treatment schedule:
    # percentages within groups:
    ntreated[month] = 0
    tx = 0
    cures = 0
    temp = S$D0[i+1]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D0[i+1] = S$D0[i+1] - temp
    temp = S$D1[i+1]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D1[i+1] = S$D1[i+1] - temp
    temp = S$D2[i+1]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D2[i+1] = S$D2[i+1] - temp
    temp = S$D3[i+1]*intervention$pD; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$D3[i+1] = S$D3[i+1] - temp
    S$Du[i+1] = S$Du[i+1] + cures
    
    cures = 0
    temp = S$J0[i+1]*intervention$pJ1; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J0[i+1] = S$J0[i+1] - temp # sic. States cannot be distinguished
    temp = S$J1[i+1]*intervention$pJ1; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J1[i+1] = S$J1[i+1] - temp
    temp = S$J2[i+1]*intervention$pJ2; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J2[i+1] = S$J2[i+1] - temp
    temp = S$J3[i+1]*intervention$pJ3; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$J3[i+1] = S$J3[i+1] - temp
    S$Ju[i+1] = S$Ju[i+1] + cures
    
    cures = 0
    temp = S$X0[i+1]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X0[i+1] = S$X0[i+1] - temp
    temp = S$X1[i+1]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X1[i+1] = S$X1[i+1] - temp
    temp = S$X2[i+1]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X2[i+1] = S$X2[i+1] - temp
    temp = S$X3[i+1]*intervention$pX; tx = tx + temp; cures = cures + temp*parameters$theta_1; S$X3[i+1] = S$X3[i+1] - temp
    S$Xu[i+1] = S$Xu[i+1] + cures
    
    # numbers
    if(intervention$nJ>0)
    {
      cures = 0
      nleftover = intervention$nJ
      treatgroups = c(3)
      if(intervention$nJ>S$J3[i+1])treatgroups = c(3,2)
      if(intervention$nJ>(S$J3[i+1]+S$J2[i+1]))treatgroups = c(3,2,1,0)
      if(sum(treatgroups==3)==1)
      {
        temp = nleftover # use up all the treatments on them
        if(S$J3[i+1]<temp)temp = S$J3[i+1] # still got some treatment left after treating all the stage 3s!
        nleftover = nleftover - temp
        tx = tx + temp
        cures = cures + temp*parameters$theta_1
        S$J3[i+1] = S$J3[i+1] - cures
        S$Ju[i+1] = S$Ju[i+1] + cures
      }
      if(sum(treatgroups==2)==1)
      {
        temp = nleftover # use up all the treatments on them
        if(S$J2[i+1]<temp)temp = S$J2[i+1] # still got some treatment left after treating all the stage 2s!
        nleftover = nleftover - temp
        tx = tx + temp
        cures = cures + temp*parameters$theta_1
        S$J2[i+1] = S$J2[i+1] - cures
        S$Ju[i+1] = S$Ju[i+1] + cures
      }
      if(sum(treatgroups==1)==1) # then treat stage 0s and 1s 'randomly'
      {
        temp = nleftover # use up all the treatments on them
        if((S$J1[i+1]+S$J0[i+1])<=nleftover){temp1 = S$J1[i+1]; temp0 = S$J0[i+1]} # still got some treatment left after treating all the stage 0s and 1s
        if((S$J1[i+1]+S$J0[i+1])>nleftover)
        {
          temp1 = nleftover*S$J1[i+1]/(S$J1[i+1]+S$J0[i+1]);
          temp0 = nleftover*S$J0[i+1]/(S$J1[i+1]+S$J0[i+1])
        }
        nleftover = nleftover - temp0 - temp1
        tx = tx + temp0 + temp1
        cures = cures + (temp0+temp1)*parameters$theta_1
        S$J0[i+1] = S$J0[i+1] - temp0*parameters$theta_1
        S$J1[i+1] = S$J1[i+1] - temp1*parameters$theta_1
        S$Ju[i+1] = S$Ju[i+1] + cures
      }
    }
    ntreated[month] = ntreated[month] +tx
  }
  R = list(
    D=c(rates_Du,rates_D0,rates_D1,rates_D2,rates_D3,rates_D4,rates_D5,rates_D6),
    J=c(rates_Ju,rates_J0,rates_J1,rates_J2,rates_J3,rates_J4,rates_J5,rates_J6),
    X=c(rates_Xu,rates_X0,rates_X1,rates_X2,rates_X3,rates_X4,rates_X5,rates_X6)
  )
  S$D=cbind(S$Du,S$D0,S$D1,S$D2,S$D3,S$D4,S$D5,S$D6)
  S$J=cbind(S$Ju,S$J0,S$J1,S$J2,S$J3,S$J4,S$J5,S$J6)
  S$X=cbind(S$Xu,S$X0,S$X1,S$X2,S$X3,S$X4,S$X5,S$X6)
  
  return(list(states=S,rates=R,treats = ntreated))
}

create_configuration = function(ny=50)
{
  configuration=list(
    nyears = ny,
    nmonths = ny*12,
    nsteps = ny*360
  )
  return(configuration)
}



parameters = list(
  #======= parameters from prison data ===========
  alpha = 0.41/360,     #0.5374/360,   # arrest rate per current IDU per year
  pi    = 0.6408    ,   # recidivism probability per jailed former user per incarceration
  rho   = 0.4714/360,   # release rate per jailed former user per year
  mu    = 0.0120/360,   # mortality rate per person per year
  beta  = 1861/360  ,   # initiation (or "birth") rate per year
  tau   = 0.0000492/360*0.61,  # transmission rate, per infected-uninfected pair of current IDU per year 0.0000331
  
  #======= additional parameters ======
  kappa = 0.26, # proportion of infections leading to spontaneous clearance (Uniform (0.22, 0.29)) (NICE)
  iota  = 360/2, # duration in which spontaneous clearers are infectious
  # transition probabilities per year (UK NICE)
  lambda_1 = 0.025/360, #rbeta(1, 38.0859, 1485.3516)/360   mild to moderate
  lambda_2 = 0.037/360, #rbeta(1, 26.905,700.2582)/360      moderate to comp cirrhosis
  lambda_3 = 0.039/360, #rbeta(1, 14.6168,360.1732)/360     comp cirrhosis to decomp cirrhosis
  lambda_4 = 0.014/360, #rbeta(1, 1.9326,136.1074)/360      decomp cirrhosis to HCC
  lambda_5 = 0.03/360 , #rbeta(1, 6.5256,210.9945)/360      HCC to LT
  
  mu_4 = 0.13/360,      #rbeta(1, 147.03,983.97)/360   decomp to death
  mu_5 = 0.43/360,      #rbeta(1, 117.1033,155.23)/360 HCC to death
  mu_6 = 0.21/360,      #rbeta(1, 16.2762,61.2294)/360 LT to death
  
  theta_1 = 0.97,       #cure rate (SVR) for genotype 1, mild/moderate for Harvoni  ##0.45 (NICE) runif(1,0.4,0.5)                         
  theta_2 = 0.94        #cure rate (SVR) for genotype 1, cirrhosis for Harvoni      ##0.25 (NICE)    
)

make_initial_states=function(configuration,equilibrium=TRUE,parameters)
{
  N=configuration$nsteps
  output = list()
  output$Du = output$D0 = output$D1 = output$D2 = output$D3 = output$D4 = output$D5 = output$D6 = rep(0,N)
  output$Ju = output$J0 = output$J1 = output$J2 = output$J3 = output$J4 = output$J5 = output$J6 = rep(0,N)
  output$Xu = output$X0 = output$X1 = output$X2 = output$X3 = output$X4 = output$X5 = output$X6 = rep(0,N)
  output$t = output$total6 = output$total5 = output$total4 = rep(0,N)
  
  if(!equilibrium)
  {
    output$Du[1] = 8464
    output$D0[1] = 1000 
    output$Ju[1] = 9705
    output$Xu[1] = 90000#136943
  }
  
  if(equilibrium)#haven't done this yet
  {
    load('equilibrium_12Feb24.rdata')
    output$Du[1] = equilibrium$Du; output$Ju[1] = equilibrium$Ju; output$Xu[1] = equilibrium$Xu
    output$D0[1] = equilibrium$D0; output$J0[1] = equilibrium$J0; output$X0[1] = equilibrium$X0
    output$D1[1] = equilibrium$D1; output$J1[1] = equilibrium$J1; output$X1[1] = equilibrium$X1
    output$D2[1] = equilibrium$D2; output$J2[1] = equilibrium$J2; output$X2[1] = equilibrium$X2
    output$D3[1] = equilibrium$D3; output$J3[1] = equilibrium$J3; output$X3[1] = equilibrium$X3
    output$D4[1] = equilibrium$D4; output$J4[1] = equilibrium$J4; output$X4[1] = equilibrium$X4
    output$D5[1] = equilibrium$D5; output$J5[1] = equilibrium$J5; output$X5[1] = equilibrium$X5
    output$D6[1] = equilibrium$D6; output$J6[1] = equilibrium$J6; output$X6[1] = equilibrium$X6
  }
  output$total4[1] = output$D4[1] + output$J4[1] + output$X4[1]
  output$total5[1] = output$D5[1] + output$J5[1] + output$X5[1]
  output$total6[1] = output$D6[1] + output$J6[1] + output$X6[1]
  return(output)
}

final_prevalence = function(X) # among jailed people
{
  i = length(X$states$D6)
  infected = sum(X$states$J[i,-1])
  total = sum(X$states$J[i,])
  return(infected/total)
}


store_states=function(X,outfile = 'equilibrium_Bo.rdata')
{
  i = length(X$states$Du)
  equilibrium = list(
    Du = X$states$Du[i], Ju = X$states$Ju[i], Xu = X$states$Xu[i],
    D0 = X$states$D0[i], J0 = X$states$J0[i], X0 = X$states$X0[i],
    D1 = X$states$D1[i], J1 = X$states$J1[i], X1 = X$states$X1[i],
    D2 = X$states$D2[i], J2 = X$states$J2[i], X2 = X$states$X2[i],
    D3 = X$states$D3[i], J3 = X$states$J3[i], X3 = X$states$X3[i],
    D4 = X$states$D4[i], J4 = X$states$J4[i], X4 = X$states$X4[i],
    D5 = X$states$D5[i], J5 = X$states$J5[i], X5 = X$states$X5[i],
    D6 = X$states$D6[i], J6 = X$states$J6[i], X6 = X$states$X6[i]
  )
  save(equilibrium,file=outfile)
}


tabulator = function(X,intervention,file='output/table_new.csv',header=FALSE,resolution='LOW')
{
  if(header)
  {
    cat('Strategy,t,Du,D0,D1,D2,D3,D4,D5,D6,Ju,J0,J1,J2,J3,J4,J5,J6,Xu,X0,X1,X2,X3,X4,X5,X6,ctx,C4,C5,C6\n',file=file)
    return(NULL)
  }
  if(!header)
  {
    if(resolution=='HIGH')i = seq(1,length(X$states$t),30)
    if(resolution=='HIGH')j = seq(1,length(X$treats),1)
    if(resolution=='LOW')i = seq(1,length(X$states$t),360)
    if(resolution=='LOW')j = seq(12,length(X$treats),12)
    output = cbind(X$states$t[i],X$states$D[i,],X$states$J[i,],X$states$X[i,])
    temp1 = rep(intervention$strategy,length(i))
    temp2 = c(0,cumsum(X$treats)[j])
    temp3 = cbind(X$states$total4,X$states$total5,X$states$total6)[i,]
    output = cbind(temp1,output,temp2,temp3)
    write.table(output,file = file,sep=',',col.names = FALSE,row.names = FALSE,append = TRUE)
  }
}




