setwd("/mnt/raid5/ChIPseqArchive/20220106_CutandRun/MethodsPaper")
require(Rsamtools);require(plyr);library(reshape2); library(dplyr); library(ggplot2); library(ggsignif); library(ggpubr);library(Rmisc);library(tidyr);
library(InPAS); library(remotes); library(GenomicRanges);library(Gviz);library(rtracklayer);library(Rsamtools);library(GenomeInfoDb)

options(ucscChromosomeNames=FALSE)
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

H3K9me3_CutRun<-renameChrs(keepSeqlevels(import.bedGraph("X1.S2.X1.NF54.MC.K9.A_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K9me3_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("H3K9me3_ChIP_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
HP1_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("Voss_GSM2743114_Pf-3D7-Schizont-HP1_ChipSeq_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_CutRun<-renameChrs(keepSeqlevels(import.bedGraph("X1.S27.X1.NF54.MC.K4.G_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))
H3K4me3_ChIP<-renameChrs(keepSeqlevels(import.bedGraph("Stunnenberg_GSM588522_Pf-3D7-H3K4me3_40hpi_ChipSeq_SPMR_FE_comp.bdg", genome="Pf3D7"),value=paste0("Pf3D7_",sprintf("%02d", 1:14),"_v3"),pruning.mode="coarse"))

H3K9me3_CutRun_norm<-H3K9me3_CutRun
H3K9me3_ChIP_norm<-H3K9me3_ChIP
HP1_ChIP_norm<-HP1_ChIP
H3K4me3_ChIP_norm<-H3K4me3_ChIP
H3K4me3_CutRun_norm<-H3K4me3_CutRun

H3K9me3_CutRun_norm$score<-H3K9me3_CutRun_norm$score/1.0
H3K9me3_ChIP_norm$score<-H3K9me3_ChIP_norm$score/4.6
HP1_ChIP_norm$score<-HP1_ChIP_norm$score/15
H3K4me3_CutRun_norm$score<-H3K4me3_CutRun_norm$score/0.18
H3K4me3_ChIP_norm$score<-H3K4me3_ChIP_norm$score/2.7

GFF<-sort(import.gff3("../genome/PlasmoDB-51_Pfalciparum3D7.gff"))
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

plotme<-data.frame(sort(subsetByOverlaps(gff[gff$feature%in% c("CDS","ncRNA"),],H3K9me3_CutRun_norm[H3K9me3_CutRun_norm@seqnames == "Chr08",])))
plotme$STRAND<-"plus"
plotme$STRAND[plotme$strand=="-"]<-"minus"

CDS.track<-GeneRegionTrack(
  plotme,
  chr="Chr08",
  #stacking = "squish",
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

genomeAxis <- GenomeAxisTrack(
                      name="chr8",littleTicks=F,col="black",showTitle=T,rotation.title=0,
                      fill.range="transparent",fontcolor="black",labelPos="below",cex.title=0.70,col.border.title="white",font.size=14,cex.title=1
                      ) 

H3K9me3_CutRun <- DataTrack(
                                      H3K9me3_CutRun_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                      type="histogram",col.histogram="green4",window=-1,windowSize = 1000,
                                      name="H3K9me3\nCUT&RUN", ylim=c(-.1,1.1), 
                                      ) 
H3K9me3_ChIP<- DataTrack(
                                            H3K9me3_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                            type="histogram",col.histogram="green2",window=-1, windowSize = 1000,
                                            name="H3K9me3\nChIP",fontsize.legend=22,ylim=c(-.1,1.1), 
                                            )
HP1_ChIP<- DataTrack(
                                        HP1_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                        type="histogram",col.histogram="springgreen2",window=-1, windowSize = 1000,
                                        name="HP1\nChIP",fontsize.legend=22,ylim=c(-.1,1.1),
                                        ) 
H3K4me3_CutRun<- DataTrack(
                                H3K4me3_CutRun_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                type="histogram",col.histogram="purple4",window=-1, windowSize = 1000,
                                name="H3K4me3\nCUT&RUN",fontsize.legend=22,ylim=c(-.1,1.1),
                                )

H3K4me3_ChIP <- DataTrack(
                                    H3K4me3_ChIP_norm,ucscChromosomeNames=FALSE,chromosome = "Chr08",
                                    type="histogram",col.histogram="purple",window=-1, windowSize = 1000,
                                    name="H3K4me3\nChIP",fontsize.legend=22,ylim=c(-.1,1.1), 
                                    )

plotTracks(c(H3K9me3_CutRun,H3K9me3_ChIP,HP1_ChIP,H3K4me3_CutRun,H3K4me3_ChIP,CDS.track,genomeAxis),
           from=398114,to=484205,
           background.panel = "white", background.title = "white",col.title = "black",col.axis = "black",fontsize.legend=18, sizes=c(1.5,1.5,1.5,1.5,1.5,2,0.5)
)

