#### DERIVE OBSERVED MUTSPEC IN CANCER
rm(list=ls(all=TRUE))
setwd("~/signatures/")

library('dplyr')
library('gtools')

##### read observed and merge it with annotation in order to annotate
ObsAll = read.table("./data/external/CancerDataFromCampbell/mtDNA_snv_Oct2016.txt", sep = '\t', header = TRUE)
ObsAll=select(ObsAll,position,ref,var,tissue,Annot,X,X.1,is_nonsense)

ExtAnn = read.table("./data/external/CancerDataFromCampbell/HumanMtDnaRefSeq.fasta.ExtensiveMut.vcf.ann", sep = '\t', header = FALSE)
ExtAnn = ExtAnn[c(2,4,5,8)]
names(ExtAnn)=c('position','ref','var','annotation')
After = ExtAnn$ref[c(4:nrow(ExtAnn))]; After = c(After,'Na','Na','Na'); length(After)
Before = ExtAnn$ref[c(1:(nrow(ExtAnn)-3))]; Before = c('Na','Na','Na',Before); length(Before)
ExtAnn$Before = Before; ExtAnn$After = After;  
ExtAnn$Context = paste(ExtAnn$Before,ExtAnn$ref,ExtAnn$After,sep='')
ExtAnn$Subst = paste(ExtAnn$ref,ExtAnn$var,sep='')
ExtAnn$names = paste(ExtAnn$Subst, ExtAnn$Context , sep=': ') 
ExtAnn = ExtAnn[c(4:(nrow(ExtAnn)-3)),]

##### MERGE ObsAll & ExtAnn (all.x = TRUE)
dim(ObsAll)
ObsAll = merge(ObsAll,ExtAnn, by = c('position','ref','var'))
dim(ObsAll)
ObsAll = ObsAll[order(ObsAll$position),]  
extractGene<-function(x) {unlist(strsplit(x,"\\|"))[4];}; 
ObsAll$ProtCodGene = apply(as.matrix(ObsAll$annotation),1,FUN=extractGene); table(ObsAll$ProtCodGene)
extractSubstType<-function(x) {unlist(strsplit(x,"\\|"))[2];} # synonymous_variant
ObsAll$SubstType = apply(as.matrix(ObsAll$annotation),1,FUN=extractSubstType); table(ObsAll$SubstType)
ObsAll = select(ObsAll,Subst,names,ProtCodGene,SubstType)

##### DEFINE MAIN SUBSETS OF MUTATIONS:
SubsetGenes = c('withoutNd6')  # c('withNd6','withoutNd6','Nd6Only')
SubsetSubstitutions = c('synonymous','AllExceptStopGains') #  SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons')

### PREPARE ARTIFICIAL MATRIX WITH 192-COMPONENT CONTEXT:
Context = data.frame(permutations(4, 3, repeats.allowed = TRUE)); names(Context) = c('Before','Middle','After')
Subst = data.frame(permutations(4, 2, repeats.allowed = TRUE)); names(Subst) = c('ref','var')
Context = merge(Context,Subst)
Context[] <- lapply(Context, function(x) (gsub("1", "A", x))); Context[] <- lapply(Context, function(x) (gsub("2", "T", x))); 
Context[] <- lapply(Context, function(x) (gsub("3", "G", x))); Context[] <- lapply(Context, function(x) (gsub("4", "C", x)))
Context = Context[Context$ref != Context$var,];
Context = Context[Context$ref == Context$Middle,]; nrow(Context) # 192
Context$names = paste(Context$ref,Context$var,': ',Context$Before,Context$Middle,Context$After,sep = '')
Context$NumbObs = 0; Artificial192=select(Context,names,NumbObs)
Context = select(Context,names,NumbObs)

##### RUN ALL FILTERS FOR SUBSETS:
PcBeforeFilters = ObsAll

MutSpec12Fin = data.frame()
MutSpec192Fin = data.frame()

for (genes in SubsetGenes) # SubsetGenes = c('withNd6','WithoutNd6','Nd6Only')
{
for (subst in SubsetSubstitutions) # SubsetSubstitutions = c('synonymous','non-synonymous','SynAndNons')
{
  # genes = 'withoutNd6'; subst = 'AllExceptStopGains'
  # if (genes == 'withNd6')  {}
  if (genes == 'withoutNd6' & subst == 'synonymous')          {
    Pc = PcBeforeFilters[PcBeforeFilters$ProtCodGene != 'ND6',];
    Pc = Pc[Pc$SubstType == 'synonymous_variant',]            } 
  if (genes == 'withoutNd6' & subst == 'AllExceptStopGains')  {
    Pc = PcBeforeFilters[PcBeforeFilters$ProtCodGene != 'ND6',];
    Pc = Pc[Pc$SubstType != 'stop_gained',]                   } 
  # if (genes == 'Nd6Only')  {Pc = PcBeforeFilters[PcBeforeFilters$ProtCodGene == 'ND6',];}    
    
######### MutSpec12
MutSpec12 = as.data.frame(table(Pc$Subst))
names(MutSpec12)=c('Subst','NumbObs')
MutSpec12 = MutSpec12[order(MutSpec12$Subst),]
MutSpec12Dummy = data.frame(c('AC','AG','AT','CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG'),c(rep(0,12)))
names(MutSpec12Dummy)=c('Subst','NumbObs')
MutSpec12 = rbind(MutSpec12,MutSpec12Dummy); MutSpec12 = aggregate(MutSpec12$NumbObs, by = list(MutSpec12$Subst), FUN = sum)
names(MutSpec12)=c('Subst','NumbObs')
MutSpec12$Genes = genes;
MutSpec12$SubstType = subst;  # MutSpec12$Subst
MutSpec12Fin = rbind(MutSpec12Fin,MutSpec12)
      
######### MutSpec192
MutSpec192 = as.data.frame(table(Pc$names)); 
MutSpec192 = MutSpec192[MutSpec192$Freq > 0,] 
names(MutSpec192)=c('names','NumbObs')
MutSpec192 = rbind(MutSpec192,Context)  # add unobserved mutations (zeros)
MutSpec192 = aggregate(MutSpec192$NumbObs, by = list(MutSpec192$names),FUN = sum)
names(MutSpec192)=c('names','NumbObs')
MutSpec192 = MutSpec192[order(MutSpec192$names),]
MutSpec192$Genes = genes; MutSpec192$SubstType = subst; 
MutSpec192Fin = rbind(MutSpec192Fin,MutSpec192)
      
}}

write.table(MutSpec12Fin,'data/processed/04.MutSpecObserved.Cancers.MutSpecObs12.txt')
write.table(MutSpec192Fin,'data/processed/04.MutSpecObserved.Cancers.MutSpecObs192.txt')
