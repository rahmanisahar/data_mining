m31 <- read.csv('~/Dropbox/for_SOM/primary_data/m31_primery_data_flx_per_pc.csv',head=TRUE,sep=',')
m31[is.na(m31)]=0
set.seed(5)
dara.sc<-scale(m31[4:35])
m31.som <- som(data =dara.sc, grid=somgrid(2,3,"hexagonal"))
plot(m31.som)
