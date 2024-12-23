library(ape)
a<-read.tree("19spCUZD1.nwk")
b<-unroot(a)
write.tree(b,"19spCUZD1.nwk.tree")
