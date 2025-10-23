####################################################################################
# May 27, 2020
# Revision for manuscript
# Figure for number of complications: 37 treatment scenarios in total
####################################################################################
library(grid)

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
C4 <- cbind(TAB[,1:2], TAB[ , grepl( "C4" , names( TAB ) ) ]); dim(C4) #check
C5 <- cbind(TAB[,1:2], TAB[ , grepl( "C5" , names( TAB ) ) ]); dim(C5) #check
C6 <- cbind(TAB[,1:2], TAB[ , grepl( "C6" , names( TAB ) ) ]); dim(C6) #check

paneller = function(column=1,row=1)
{
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  
  pushViewport(plotViewport(c(1,1,1,1),
                            xscale=c(-1,50),yscale=c(0,c(13.8,1190,45)[row])))
  
  
  if(row==1)yvec <- aa<-rowSums(C4[,3:11])/1000
  if(row==2)yvec <- bb<-rowSums(C5[,3:11])/1
  #if(row==4)yvec <- cc<-rowSums(C6[,3:11])/1
  if(row==3)yvec <- dd<-(TAB$ctx)/1000
  strids = which(strategies$strategy_class==0)
  for(str in strids)liner(str,yvec)
  
  if(column==1)strids = which(strategies$strategy_class==3)
  #if(column==2)strids = which(strategies$strategy_class==1)
  if(column==2)strids = which(strategies$strategy_class==2)
  if(column==3)strids = which(strategies$strategy_class==4)
  if(column==4)strids = which(strategies$strategy_class==5) 
  #if(column==5)strids = which(strategies$strategy_class==6)
  
  for(str in strids)liner(str,yvec)
  
  grid.text(paste('(',letters[(row-1)*4+column],')',sep=''),unit(0.5,'lines'),unit(1,'npc')+unit(-1,'lines'),hjust=0)
  if(column==1)grid.text(c('Decompensated\ncirrhosis cases (000s)','Hepatocellular\ncarcinomas','Number of\ntreatments (000s)')[row],x=unit(-4,'lines'),rot=90)
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  if(row==3)grid.xaxis(at=seq(0,50,10))
  if(row!=3)grid.xaxis(at=seq(0,50,10),label=FALSE)
  
  grid.yaxis(label=c(TRUE,FALSE,FALSE,FALSE)[column])
  
  if(row==3)grid.text('Time (y)',y=unit(-3,'lines'))
  if(row==1)grid.text(c('Current PWID (%)','Prisoners (all)','Former PWID (%)','Combined')[column],y=unit(1,'npc')+unit(1,'lines'))
  #grid.text(GC$ylab[GC$state==state],x=unit(-4,'lines'),rot=90)
  popViewport()
  popViewport()
}

legend = function(column)
{
  row = c(2,2,2,2)[column]
  pushViewport(viewport(layout.pos.col=column,layout.pos.row=row))
  pushViewport(plotViewport(c(1,1,1,1)))
  if(column==1)strids = which(strategies$strategy_class==3)
  if(column==2)strids = which(strategies$strategy_class==2)
  if(column==3)strids = which(strategies$strategy_class==4)
  if(column==4)strids = which(strategies$strategy_class==5) 
  #if(column==5)strids = which(strategies$strategy_class==6)
  
  SC=strategies$strategy_class
  for(str in strids)
  {
    z=0.3
    n = sum(SC==SC[str])
    x = sum(SC[1:str]==SC[str])
    if(SC[str]==0)cl = 'black'
    if(SC[str]==2)cl = colorRampPalette(c('purple','plum','violet'))(n)[x]
    if(SC[str]==3)cl = colorRampPalette(c('yellowgreen','darkgreen'))(n)[x]
    if(SC[str]==4)cl = colorRampPalette(c('skyblue','navyblue'))(n)[x]
    if(SC[str]==5)cl = colorRampPalette(c('gold','orange','red','darkred'))(n)[x]
    #if(SC[str]==6)cl = colorRampPalette(c('lightgrey','darkgrey'))(n)[x]
    
    if(column!=2)
    {
      grid.lines(z+c(0.1,0.1),0.3-c(x-1,x)/n*0.25,gp=gpar(col=cl))
      if(x==1)grid.lines(z+c(0.1,0.12),0.3-c(x-1,x-1)/n*0.25,gp=gpar(col=cl))
      if(x==n)grid.lines(z+c(0.1,0.12),0.3-c(x,x)/n*0.25,gp=gpar(col=cl))
      if(x==n)grid.text(c('15%/y','','15%/y','13%/y')[column],unit(z+0.12,'npc')+unit(0.5,'lines'),0.05,hjust=0,gp=gpar(cex=0.75))
      if(x==1)grid.text(c('1%/y','','1%/y','10%/y')[column],unit(z+0.12,'npc')+unit(0.5,'lines'),0.3,hjust=0,gp=gpar(cex=0.75))
    }
    if(column==2)
    {
      grid.lines(z+c(0.075,0.15),unit(rep(c(0.5,1.5,2.5)[x],2),'lines'),gp=gpar(col=cl))
      grid.text(c('Up to C.C.','Moderate/C.C.','C.C.')[x],unit(z+0.15,'npc')+unit(0.5,'lines'),unit((x-0.5)/0.75,'lines'),hjust=0,gp=gpar(cex=0.75))
      
    }
  }
  
  
  popViewport()
  popViewport()
}

png('output/fig/Figure_2.png',height=11,width=15,units='cm',res=900,pointsize=8)
pushViewport(plotViewport(c(3,5,2,0)))
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=4)))

for(i in 1:4)
{
  for(j in 1:3)paneller(column=i,row=j)
  legend(i)
}

popViewport()
popViewport()
dev.off()

