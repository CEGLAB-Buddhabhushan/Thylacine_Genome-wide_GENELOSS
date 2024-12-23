library(ape)
a<-read.tree("Species.no_internodes.nwk")
b<-unroot(a)
write.tree(b,"Species.no_internodes.unroot.nwk")
