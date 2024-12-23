library(ape)
a<-read.tree("SAMD9.nwk")
b<-unroot(a)
write.tree(b,"SAMD9.nwk.tree")
