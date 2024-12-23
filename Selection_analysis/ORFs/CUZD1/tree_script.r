library(ape)
a<-read.tree("CUZD1.nwk")
b<-unroot(a)
write.tree(b,"CUZD1.nwk.tree")
