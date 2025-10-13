library(grid)

cl=c('steelblue','yellow2','orange','orangered','darkred','purple','darkblue','black')

png('output/figure_1_D.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(4,4,1,1),
                          xscale=c(0,max(X$states$t)/360),yscale=c(0,1.05*max(X$states$Du)/1000)))
grid.lines(X$states$t/360,X$states$Du/1000,default.units = 'native',gp=gpar(col=cl[1]))
grid.lines(X$states$t/360,X$states$D0/1000,default.units = 'native',gp=gpar(col=cl[2]))
grid.lines(X$states$t/360,X$states$D1/1000,default.units = 'native',gp=gpar(col=cl[3]))
grid.lines(X$states$t/360,X$states$D2/1000,default.units = 'native',gp=gpar(col=cl[4]))
grid.lines(X$states$t/360,X$states$D3/1000,default.units = 'native',gp=gpar(col=cl[5]))
grid.lines(X$states$t/360,X$states$D4/1000,default.units = 'native',gp=gpar(col=cl[6]))
grid.lines(X$states$t/360,X$states$D5/1000,default.units = 'native',gp=gpar(col=cl[7]))
grid.lines(X$states$t/360,X$states$D6/1000,default.units = 'native',gp=gpar(col=cl[8]))

grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
grid.xaxis()
grid.yaxis()
grid.text('Time (y)',y=unit(-3,'lines'))
grid.text('Drug users out of jail (000s)',x=unit(-3,'lines'),rot=90)
popViewport()
dev.off()


png('output/figure_1_J.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(4,4,1,1),
                          xscale=c(0,max(X$states$t)/360),yscale=c(0,1.05*max(X$states$Ju)/1000)))
grid.lines(X$states$t/360,X$states$Ju/1000,default.units = 'native',gp=gpar(col=cl[1]))
grid.lines(X$states$t/360,X$states$J0/1000,default.units = 'native',gp=gpar(col=cl[2]))
grid.lines(X$states$t/360,X$states$J1/1000,default.units = 'native',gp=gpar(col=cl[3]))
grid.lines(X$states$t/360,X$states$J2/1000,default.units = 'native',gp=gpar(col=cl[4]))
grid.lines(X$states$t/360,X$states$J3/1000,default.units = 'native',gp=gpar(col=cl[5]))
grid.lines(X$states$t/360,X$states$J4/1000,default.units = 'native',gp=gpar(col=cl[6]))
grid.lines(X$states$t/360,X$states$J5/1000,default.units = 'native',gp=gpar(col=cl[7]))
grid.lines(X$states$t/360,X$states$J6/1000,default.units = 'native',gp=gpar(col=cl[8]))

grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
grid.xaxis()
grid.yaxis()
grid.text('Time (y)',y=unit(-3,'lines'))
grid.text('Jailed drug users (000s)',x=unit(-3,'lines'),rot=90)
popViewport()
dev.off()



png('output/figure_1_X.png',height=8,width=8,units='cm',res=300,pointsize=10)
pushViewport(plotViewport(c(4,4,1,1),
                          xscale=c(0,max(X$states$t)/360),yscale=c(0,1.05*max(X$states$Xu)/1000)))
grid.lines(X$states$t/360,X$states$Xu/1000,default.units = 'native',gp=gpar(col=cl[1]))
grid.lines(X$states$t/360,X$states$X0/1000,default.units = 'native',gp=gpar(col=cl[2]))
grid.lines(X$states$t/360,X$states$X1/1000,default.units = 'native',gp=gpar(col=cl[3]))
grid.lines(X$states$t/360,X$states$X2/1000,default.units = 'native',gp=gpar(col=cl[4]))
grid.lines(X$states$t/360,X$states$X3/1000,default.units = 'native',gp=gpar(col=cl[5]))
grid.lines(X$states$t/360,X$states$X4/1000,default.units = 'native',gp=gpar(col=cl[6]))
grid.lines(X$states$t/360,X$states$X5/1000,default.units = 'native',gp=gpar(col=cl[7]))
grid.lines(X$states$t/360,X$states$X6/1000,default.units = 'native',gp=gpar(col=cl[8]))

grid.lines(c(0,1,1,0,0),c(0,0,1,1,0))
grid.xaxis()
grid.yaxis()
grid.text('Time (y)',y=unit(-3,'lines'))
grid.text('Ex-drug users (000s)',x=unit(-3,'lines'),rot=90)
popViewport()
dev.off()

rm(cl)
