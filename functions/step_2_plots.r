library(grid)

TAB = read.csv('settings/tableHIGH_1_wFINALstrategies_yichen.csv')
TAB$t = TAB$t/360 # convert to 'years'

GC = read.csv('settings/graph_config.csv',as.is=TRUE)

liner = function(str,yvec,tab=TAB)
{
  SC=strategies$strategy_class
  n = sum(SC==SC[str])
  x = sum(SC[1:str]==SC[str])
  if(SC[str]==0)cl = 'black'
  if(SC[str]==1)cl = colorRampPalette(c('gold','orange','red','darkred'))(n)[x]
  if(SC[str]==2)cl = colorRampPalette(c('violet','plum','purple'))(n)[x]
  if(SC[str]==3)cl = colorRampPalette(c('yellowgreen','darkgreen'))(n)[x]
  if(SC[str]==4)cl = colorRampPalette(c('skyblue','navyblue'))(n)[x]
  
  rows = which(tab$Strategy==str)
  grid.lines(c(-1,tab$t[rows]),c(yvec[rows][1],yvec[rows]),default.units = 'native',gp=gpar(col=cl))
}

for(state in 1:dim(GC)[1])
{
  png(paste('output/figure_1',letters[state],'.png',sep=''),height=8,width=8,units='cm',res=300,pointsize=10)
  pushViewport(plotViewport(c(5,5,1,1),
                            xscale=c(-1,30),yscale=c(0,GC$ymax[GC$state==state])))
  cat('yvec <- ',GC$columns[GC$state==state],sep='',file='temp.r')
  source('temp.r')
  for(str in 1:dim(strategies)[1])liner(str,yvec)
  
  
  file.remove('temp.r')
  
  grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
  grid.xaxis()
  grid.yaxis()
  grid.text('Time (y)',y=unit(-3,'lines'))
  grid.text(GC$ylab[GC$state==state],x=unit(-4,'lines'),rot=90)
  popViewport()
  dev.off()
  
}
