setwd("/mnt/raid5/ChIPseqArchive/20220106_CutandRun/")
require(Rsamtools);require(plyr);library(reshape2); library(dplyr); library(ggplot2); library(ggsignif); library(ggpubr);library(Rmisc);library(tidyr);
library(InPAS); library(remotes); library(GenomicRanges);library(Gviz);library(rtracklayer);library(Rsamtools);library(GenomeInfoDb)
options(ucscChromosomeNames=FALSE)

install.packages("./genome/BSGenome.Pf3D7.PlasmoDB.51_1.0.tar.gz", repos = NULL, type="source")

library(BSGenome.Pf3D7.PlasmoDB.51)
Pf3D7.seqinfo<-data.frame(cbind(Chr=seqlevels(BSGenome.Pf3D7.PlasmoDB.51),chr=c(paste0("Chr",substring(names(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))[1:14],7,8)),"Mito","Apico"), data.frame(seqinfo(BSGenome.Pf3D7.PlasmoDB.51))))
genome=Pf3D7

renameChrs<-function(gr){  
  gr<-data.frame(gr)
  gr$seqnames<-as.character(gr$seqnames)
  if(grepl("Pf3D7_",gr$seqnames[1],fixed=T)){
    gr$seqnames<-Pf3D7.seqinfo$chr[match(gr$seqnames,Pf3D7.seqinfo$Chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$chr)]
    genome(gr)<-"51"
  }else{
    gr$seqnames<-Pf3D7.seqinfo$Chr[match(gr$seqnames,Pf3D7.seqinfo$chr)]
    gr<-makeGRangesFromDataFrame(data.frame(gr),keep.extra.columns = T)
    seqlengths(gr)<-Pf3D7.seqinfo$seqlengths[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    isCircular(gr)<-Pf3D7.seqinfo$isCircular[match(names(seqlengths(gr)),Pf3D7.seqinfo$Chr)]
    genome(gr)<-"51"
  }
  return(gr)
}



GFF<-sort(import.gff3("./genome/PlasmoDB-51_Pfalciparum3D7.gff"))
seqlengths(GFF)[1:16]<-seqlengths(Pf3D7)[c(1:14,16,15)]
isCircular(GFF)<-c(rep(F,14),T,T)
genome(GFF)<-genome(Pf3D7)
seqinfo(GFF)

# modify GFF for compatibility with later functions
gff<-renameChrs(GFF)
gff$feature<-gff$type
names(gff)<-gff$ID

gff$gene<-substr(gff$ID,regexpr("PF3D7_",gff$ID),regexpr("PF3D7_",gff$ID)+12)
gff[gff@seqnames %in% c("Apico","Mito")]$gene<-substr(gff[gff@seqnames %in% c("Apico","Mito")]$ID,regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID),regexpr("PF3D7_",gff[gff@seqnames %in% c("Apico","Mito")]$ID)+13)
gff$transcript<-gff$ID
gff[gff$type %in% c("protein_coding_gene","pseudogene","ncRNA_gene")]$transcript<-NA
gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$transcript<-gff[gff$type %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")]$Parent
gff$exon<-NA; gff$exon[gff$feature=="exon"]<-gff$ID[gff$feature=="exon"]
gff$source<-NULL; gff$phase<-NULL; gff$Note<-NULL; gff$type<-NULL; gff$protein_source_id<-NULL; gff$score<-NULL

#read in FE bedgraphs
H3K9me3_CutRun_rep1<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S9.X1.NF54.MC.K9.E_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K9me3_CutRun_rep2<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S18.X1.NF54.MC.K9.F_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K9me3_CutRun_rep3<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S26.X1.NF54.MC.K9.G_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K9me3_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/H3K9me3_ChIP_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
HP1_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/HP1_ChIP_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_CutRun_rep1<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S10.X1.NF54.MC.K4.E_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_CutRun_rep2<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S19.X1.NF54.MC.K4.F_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_CutRun_rep3<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/X1.S27.X1.NF54.MC.K4.G_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/H3K4me3_ChIP_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

S26_0.25<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S26_K9_0.25_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
S26_0.0625<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S26_K9_0.0625_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
S26_0.015625<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S26_K9_0.015625_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
S27_0.25<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S27_K4_0.25_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
S27_0.0625<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S27_K4_0.0625_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
S27_0.015625<-renameChrs(keepSeqlevels(import.bedGraph("./singleSamplePeaks/SPMR/S27_K4_0.015625_SPMR.comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))


#
H3K9me3_CutRun_rep1_norm<-H3K9me3_CutRun_rep1
H3K9me3_CutRun_rep2_norm<-H3K9me3_CutRun_rep2
H3K9me3_CutRun_rep3_norm<-H3K9me3_CutRun_rep3
H3K9me3_ChIP_norm<-H3K9me3_ChIP
HP1_ChIP_norm<-HP1_ChIP
H3K4me3_CutRun_rep1_norm<-H3K4me3_CutRun_rep1
H3K4me3_CutRun_rep2_norm<-H3K4me3_CutRun_rep2
H3K4me3_CutRun_rep3_norm<-H3K4me3_CutRun_rep3
H3K4me3_ChIP_norm<-H3K4me3_ChIP

S26_0.25_norm<-S26_0.25
S26_0.0625_norm<-S26_0.0625
S26_0.015625_norm<-S26_0.015625
S27_0.25_norm<-S27_0.25
S27_0.0625_norm<-S27_0.0625
S27_0.015625_norm<-S27_0.015625

#lazy scaling to 1ish
H3K9me3_CutRun_rep1_norm$score<-H3K9me3_CutRun_rep1_norm$score/9.5
H3K9me3_CutRun_rep2_norm$score<-H3K9me3_CutRun_rep2_norm$score/7
H3K9me3_CutRun_rep3_norm$score<-H3K9me3_CutRun_rep3_norm$score/14.2
H3K9me3_ChIP_norm$score<-H3K9me3_ChIP_norm$score/4.8
HP1_ChIP_norm$score<-HP1_ChIP_norm$score/8
H3K4me3_CutRun_rep1_norm$score<-H3K4me3_CutRun_rep1_norm$score/2.6
H3K4me3_CutRun_rep2_norm$score<-H3K4me3_CutRun_rep2_norm$score/2.4
H3K4me3_CutRun_rep3_norm$score<-H3K4me3_CutRun_rep3_norm$score/2.5
H3K4me3_ChIP_norm$score<-H3K4me3_ChIP_norm$score/2.7
S26_0.25_norm$score<-S26_0.25_norm$score/15
S26_0.0625_norm$score<-S26_0.0625_norm$score/15
S26_0.015625_norm$score<-S26_0.015625_norm$score/14
S27_0.25_norm$score<-S27_0.25_norm$score/2.3
S27_0.0625_norm$score<-S27_0.0625_norm$score/2.3
S27_0.015625_norm$score<-S27_0.015625_norm$score/2.8

plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],H3K9me3_CutRun_rep1_norm)))
plotme$STRAND<-"plus"
plotme$STRAND[plotme$strand=="-"]<-"minus"

CDS.track<-GeneRegionTrack(
  plotme,
  chr="Chr08",
  feature=plotme$STRAND,
  col="black",
  name="genes",
  transcriptAnnotation = "gene",
  plus="blue",
  minus="red",
  cex.group=0.75,
  cex.title=1,
  just.group="right",
  fontsize.group=10,
  fontcolor.group="black"
)

genomeAxis.track <- GenomeAxisTrack(
                      name="chr8",littleTicks=F,col="black",showTitle=T,rotation.title=0,
                      fill.range="transparent",fontcolor="black",labelPos="below",cex.title=0.70,col.border.title="white",font.size=14,cex.title=1
                      ) 

H3K9me3_CutRun_rep1.track <- DataTrack(
                                      H3K9me3_CutRun_rep1_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                      type="histogram",col.histogram="green4",window=-1,windowSize = 1000,
                                      name="H3K9me3\nCUT&RUN", ylim=c(-.1,1.1), 
                                      ) 
H3K9me3_CutRun_rep2.track <- DataTrack(
                                      H3K9me3_CutRun_rep2_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                      type="histogram",col.histogram="green4",window=-1,windowSize = 1000,
                                      name="H3K9me3\nCUT&RUN", ylim=c(-.1,1.1), 
                                      ) 
H3K9me3_CutRun_rep3.track <- DataTrack(
                                      H3K9me3_CutRun_rep3_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                      type="histogram",col.histogram="green4",window=-1,windowSize = 1000,
                                      name="H3K9me3\nCUT&RUN", ylim=c(-.1,1.1), 
                                      ) 
H3K9me3_ChIP.track<- DataTrack(
                              H3K9me3_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                              type="histogram",col.histogram="green2",window=-1, windowSize = 1000,
                              name="H3K9me3\nChIP",fontsize.legend=22,ylim=c(-.1,1.1), 
                              )
HP1_ChIP.track<- DataTrack(
                          HP1_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                          type="histogram",col.histogram="springgreen2",window=-1, windowSize = 1000,
                          name="HP1\nChIP",fontsize.legend=22,ylim=c(-.1,1.1),
                          ) 
H3K4me3_CutRun_rep1.track<- DataTrack(
                                H3K4me3_CutRun_rep1_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
                                name="H3K4me3\nCUT&RUN",fontsize.legend=22, ylim=c(-.1,1.1),
                                )
H3K4me3_CutRun_rep2.track<- DataTrack(
                                H3K4me3_CutRun_rep2_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
                                name="H3K4me3\nCUT&RUN",fontsize.legend=22, ylim=c(-.1,1.1),
                              )
H3K4me3_CutRun_rep3.track<- DataTrack(
                                    H3K4me3_CutRun_rep3_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                    type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
                                    name="H3K4me3\nCUT&RUN",fontsize.legend=22, ylim=c(-.1,1.1),
                                  )
H3K4me3_ChIP.track <- DataTrack(
                              H3K4me3_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                              type="histogram",col.histogram="purple",window=-1, windowSize = 1000,
                              name="H3K4me3\nChIP",fontsize.legend=22,ylim=c(-.1,1.1), 
                              )

S26_0.25_norm.track <- DataTrack(
  S26_0.25_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="green4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)
S26_0.0625_norm.track <- DataTrack(
  S26_0.0625_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="green4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)
S26_0.015625_norm.track <- DataTrack(
  S26_0.015625_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="green4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)

S27_0.25.track <- DataTrack(
  S27_0.25_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)
S27_0.0625_norm.track <- DataTrack(
  S27_0.0625_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)
S27_0.015625_norm.track <- DataTrack(
  S27_0.015625_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
  type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
  name="H3K4me3\nChIP",fontsize.legend=22,#ylim=c(-.1,1.1), 
)


# Figure 2
plotTracks(c(H3K9me3_CutRun_rep1.track,H3K9me3_CutRun_rep2.track,H3K9me3_CutRun_rep3.track,H3K4me3_CutRun_rep1.track,H3K4me3_CutRun_rep2.track,H3K4me3_CutRun_rep3.track,CDS.track,genomeAxis.track),
           from=398114,to=484205,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)

# Figure 3
plotTracks(c(H3K9me3_CutRun_rep3.track,S26_0.25_norm.track,S26_0.0625_norm.track,S26_0.015625_norm.track,H3K4me3_CutRun_rep3.track,S27_0.25.track,S27_0.0625_norm.track,S27_0.015625_norm.track,CDS.track,genomeAxis.track),
           from=398114,to=484205,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,0.5)
)

# Figure 4
plotTracks(c(H3K9me3_CutRun_rep3.track,H3K9me3_ChIP.track,HP1_ChIP.track,H3K4me3_CutRun_rep2.track,H3K4me3_ChIP.track,CDS.track,genomeAxis.track),
           from=398114,to=484205,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,2,0.5)
)



