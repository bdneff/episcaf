seqs=read.csv("/Users/jalt0715/Dropbox/Work-TGen/Projects/structureModeling/PepSeqLibraries/DP3/content/scaffoldedEpitopeControls/scaffolded403.csv")
source("/Users/jalt0715/Dropbox/Work-TGen/Projects/structureModeling/PepSeqLibraries/DP3/content/scaffoldedEpitopeControls/functions/mutate6.R")
source("/Users/jalt0715/Dropbox/Work-TGen/Projects/structureModeling/PepSeqLibraries/DP3/content/scaffoldedEpitopeControls/functions/mutateToA.R")

originals=seqs$scaffoldEPITOPE
mut1s=unlist(lapply(seqs$scaffoldEPITOPE,function(X) mutate6(X,"PPDDGG",1,max_tries=1000)))
mut4s=unlist(lapply(seqs$scaffoldEPITOPE[1:403],function(X) mutate6(X,"PPDDGG",4,max_tries=1000)))
Ala1=unlist(lapply(seqs$scaffoldEPITOPE,function(X) mutateToA(X,1)))
Ala2=unlist(lapply(seqs$scaffoldEPITOPE,function(X) try(mutateToA(X,2))))

out=rbind(
cbind(paste0(seqs$Design_ID,""),seqs$Target,originals),
cbind(paste0(seqs$Design_ID,"_scaffoldMutX1"),seqs$Target,mut1s),
cbind(paste0(seqs$Design_ID,"_scaffoldMutX4"),seqs$Target,mut4s),
cbind(paste0(seqs$Design_ID,"_epitopeIsland1Mut>A"),seqs$Target,Ala1),
cbind(paste0(seqs$Design_ID,"__epitopeIsland2Mut>A"),seqs$Target,Ala2))

out=cbind(out,toupper(out[,3]),substring(toupper(out[,3]),1,103))
colnames(out)=c("design_ID","target","scaffoldEPITOPE","designedSequence","sequence")
out=out[which(nchar(out[,"scaffoldEPITOPE"])!=0),]

out=cbind(out[,"sequence"],"scaffoldedAbEpitope","RFD",out[,"designedSequence"],nchar(out[,"designedSequence"]),out[,"design_ID"],out[,"target"],out[,"scaffoldEPITOPE"])
colnames(out)=c("sequence","category","model","designedSequence","designedSequenceLength","design_ID","target","scaffoldEPITOPE")


write.csv(out,file="/Users/jalt0715/Desktop/scaffoldedControls.csv",quote=F,row.names=F)