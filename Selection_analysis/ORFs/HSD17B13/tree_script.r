library(ape)
a<-read.tree("HSD17B13.nwk")
b<-unroot(a)
write.tree(b,"HSD17B13.nwk.tree")
