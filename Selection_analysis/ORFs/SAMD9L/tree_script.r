library(ape)
a<-read.tree("SAMD9L.nwk")
b<-unroot(a)
write.tree(b,"SAMD9L.nwk.tree")
