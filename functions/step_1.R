####################################################################################
# May 27, 2020
# Revision for manuscript
# Figure for number of infections: 37 treatment scenarios in total
# Figure including combined strategies (total: 5 columns)
####################################################################################
library(grid)

# Boolean flag for different table versions
Bo <- F

if (Bo){
  TAB = read.csv('settings/tableHIGH_1_wFINALstrategies_bo.csv')
} else {
  TAB = read.csv('settings/tableHIGH_1_wFINALstrategies_yichen.csv')
}

TAB$t = TAB$t/360 # convert to 'years'

strategies = read.csv("settings/Final_strategies.csv")
 # "
strategies$pD = strategies$pD/12;   # Original is *loosely the annual proportion
strategies$pJ1 = strategies$pJ1/12; # "
strategies$pJ2 = strategies$pJ2/12; # "
strategies$pJ3 = strategies$pJ3/12; # "
strategies$pX = strategies$pX/12;   # "

#GC = read.csv('input/graph_config.csv',as.is=TRUE)

liner = function(str,yvec,tab=TAB)
{
  SC=strategies$strategy_class
  n = sum(SC==SC[str])
  x = sum(SC[1:str]==SC[str])
  if(SC[str]==0)cl = 'black'
  if(SC[str]==2)cl = colorRampPalette(c('purple','plum','violet'))(n)[x] #put NA to make it invisible
  if(SC[str]==3)cl = colorRampPalette(c('yellowgreen','darkgreen'))(n)[x]
  if(SC[str]==4)cl = colorRampPalette(c('skyblue','navyblue'))(n)[x]
  if(SC[str]==5)cl = colorRampPalette(c('gold','orange','red','darkred'))(n)[x]
  #if(SC[str]==6)cl = colorRampPalette(c('lightgrey','darkgrey'))(n)[x]
  
  rows = which(tab$Strategy==str)
  grid.lines(c(-1,tab$t[rows]),c(yvec[rows][1],yvec[rows]),default.units = 'native',gp=gpar(col=cl))
}
D <- cbind(TAB[,1:2], TAB[ , grepl( "D" , names( TAB ) ) ]); dim(D) #check
J <- cbind(TAB[,1:2], TAB[ , grepl( "J" , names( TAB ) ) ]); dim(J) #check
X <- cbind(TAB[,1:2], TAB[ , grepl( "X" , names( TAB ) ) ]); dim(X) #check

paneller = function(column=1,row=1)
{
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  
  pushViewport(plotViewport(c(1,1,1,1),
                            xscale=c(-1,50),yscale=c(0,c(12.8,3.9,31,47)[row])))

  if(row==1)yvec <- a <- rowSums(D[,12:74])/1000
  if(row==2)yvec <- b <- rowSums(J[,12:74])/1000
  if(row==3)yvec <- c <- rowSums(X[,12:74])/1000
  if(row==4)yvec <- d <- (rowSums(D[,12:74])+rowSums(J[,12:74])+rowSums(X[,12:74]))/1000
  
  strids = which(strategies$strategy_class==0)
  for(str in strids)liner(str,yvec)
  
  if(column==1)strids = which(strategies$strategy_class==3)
  if(column==2)strids = which(strategies$strategy_class==2)
  if(column==3)strids = which(strategies$strategy_class==4)
  if(column==4)strids = which(strategies$strategy_class==5)
  
  for(str in strids)liner(str,yvec)

  grid.text(paste('(',letters[(row-1)*4+column],')',sep=''),unit(0.5,'lines'),unit(1,'npc')+unit(-1,'lines'),hjust=0)
  if(column==1)grid.text(c('Infections (000s):\nCommunity PWID','Infections (000s):\nIncarcerated ex-PWID','Infections (000s):\nCommunity ex-PWID','Infections (000s):\nTotal')[row],x=unit(-3.5,'lines'),rot=90)
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  
  if(row==4)grid.xaxis(at=seq(0,50,10))
  if(row!=4)grid.xaxis(at=seq(0,50,10),label=FALSE)
  
  if(row!=3)grid.yaxis(label=c(TRUE,FALSE,FALSE,FALSE)[column])
  if(row==3)grid.yaxis(at=c(0,5,10,15,20,25),label=c(TRUE,FALSE,FALSE,FALSE,FALSE)[column])
  #if(row==1)grid.yaxis(at=c(0,3,6,9),label=c(TRUE,FALSE,FALSE,FALSE)[column])
  
  if(row==4)grid.text('Time (y)',y=unit(-3,'lines'))
  if(row==1)grid.text(c('Current PWID (%)','Prisoners (all)','Former PWID (%)', 'Combined')[column],y=unit(1,'npc')+unit(1,'lines'))
  #grid.text(GC$ylab[GC$state==state],x=unit(-4,'lines'),rot=90)
  popViewport()
  popViewport()
}

legend = function(column)
{
  row = c(4,4,4,4)[column]
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(c(1,1,1,1)))
  if(column==1)strids = which(strategies$strategy_class==3)
  if(column==2)strids = which(strategies$strategy_class==2)
  if(column==3)strids = which(strategies$strategy_class==4)
  if(column==4)strids = which(strategies$strategy_class==5)
  
  SC=strategies$strategy_class
  for(str in strids)
  {
    n = sum(SC==SC[str])
    x = sum(SC[1:str]==SC[str])
    if(SC[str]==0)cl = 'black'
    if(SC[str]==2)cl = colorRampPalette(c('purple','plum','violet'))(n)[x]
    if(SC[str]==3)cl = colorRampPalette(c('yellowgreen','darkgreen'))(n)[x]
    if(SC[str]==4)cl = colorRampPalette(c('skyblue','navyblue'))(n)[x]
    if(SC[str]==5)cl = colorRampPalette(c('gold','orange','red','darkred'))(n)[x]
    
    if(column!=2)
    {
      grid.lines(c(0.1,0.1),0.3-c(x-1,x)/n*0.25,gp=gpar(col=cl))
      if(x==1)grid.lines(c(0.1,0.12),0.3-c(x-1,x-1)/n*0.25,gp=gpar(col=cl))
      if(x==n)grid.lines(c(0.1,0.12),0.3-c(x,x)/n*0.25,gp=gpar(col=cl))
      if(x==n)grid.text(c('15%/y','','15%/y','13%/y')[column],unit(0.12,'npc')+unit(0.5,'lines'),0.05,hjust=0,gp=gpar(cex=0.75))
      if(x==1)grid.text(c('1%/y','','1%/y','10%/y')[column],unit(0.12,'npc')+unit(0.5,'lines'),0.3,hjust=0,gp=gpar(cex=0.75))
    }
    if(column==2)
    {
      grid.lines(c(0.075,0.15),unit(rep(c(0.5,1.5,2.5)[x],2),'lines'),gp=gpar(col=cl))
      grid.text(c('Up to C.C.','Moderate/C.C.','C.C.')[x],unit(0.15,'npc')+unit(0.5,'lines'),unit((x-0.5)/0.75,'lines'),hjust=0,gp=gpar(cex=0.75))
    }
    
  }
  
  
  popViewport()
  popViewport()
}

png('output/fig/Figure_1.png',height=15,width=15,units='cm',res=900,pointsize=8)
#pdf('output/figure_rev.pdf',height = 11, width = 8.5, paper = "letter")
pushViewport(plotViewport(c(3,5,2,0)))
pushViewport(viewport(layout=grid.layout(nrow=4,ncol=4)))

for(i in 1:4)
{
  for(j in 1:4)paneller(column=i,row=j)
  legend(i)
}

popViewport()
popViewport()
#dev.print(pdf, "Figure_rev.pdf")

dev.off()
  
