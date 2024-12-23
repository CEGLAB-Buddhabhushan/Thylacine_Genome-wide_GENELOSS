library(ape)
a<-read.tree("VWA7.nwk")
b<-unroot(a)
write.tree(b,"VWA7.nwk.tree")
