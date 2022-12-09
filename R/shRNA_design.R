#Ines Patop, 2022
#The things you can use for design of shRNAs

#function to get "not in" from %in% in dplyr
'%!in%' <- function(x,y)!('%in%'(x,y))

#' Generate shRNAs against circRNA back-splice junctions
#'
#'This script contains ONE funciton that will read a table with circRNA/exon-exon splicing junction coordinates and genenames in the format: Name Chr Start End.
#'Optional:a column named Strand
#'So far the following species are available: fly: dm3, dm6 and human: hs19, mice: mm10 and rat: rn4
#'
#'
#' @param sp character with Name of the genome to use, available are: dm6, dm3, hg19, hg38, mm10, rn4. Default is dm6
#' @param input_coordinate document with the coordinate of the circRNAs to generate shRNAs in TSV format tab seppatated. It requires a column with Chromosome, Start and End.
#' @param writetab a boolean, if writetab=T, then the function will not write any table, instead it will return the table to R
#' @param output if writetab = T, it will create the output in this path, default is "shRNAs.tsv"
#'
#' @param circRNAseq boolean, if circRNAseq=T, the final table has the putative circRNA sequence asuming there is no internal alternative splicing (ie. all the possible internal exons one after the other)
#'
#' @param shift boolean, if shift=T, return also 5´and 3´shift sequences
#'
#' @return if writetab=F, It generates a DataFrame with the original table with the shRNA deisng appended
#'
#'
#' @export
#'
#' @examples
#' #There are examples in the folder names "../Test". To run it:
#' OligoDesigner(input_coordinates = "../test/circs_totest.txt", writetab=F)
#'
#'

OligoDesigner<-function(sp="dm6",input_coordinates="circs_totest.txt", writetab=T,output="shRNAs.tsv",shift=T,circRNAseq=T){

  assertthat::assert_that(
    assertthat::see_if(assertthat::is.readable(input_coordinates))
  )


  if (sp =="dm6") {
    require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm6.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm6)
    gen<-BSgenome.Dmelanogaster.UCSC.dm6

  }  else if (sp =="dm3") {
    require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    tx<-TxDb.Dmelanogaster.UCSC.dm3.ensGene
    require(BSgenome.Dmelanogaster.UCSC.dm3)
    gen<-BSgenome.Dmelanogaster.UCSC.dm3

  } else if (sp =="hg19") {
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg19.knownGene
    require(BSgenome.Hsapiens.UCSC.hg19)
    gen<-BSgenome.Hsapiens.UCSC.hg19

  } else if (sp =="hg38") {
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    tx<-TxDb.Hsapiens.UCSC.hg38.knownGene
    require(BSgenome.Hsapiens.UCSC.hg38)
    gen<-BSgenome.Hsapiens.UCSC.hg38

  } else if (sp =="mm10") {
    require(TxDb.Mmusculus.UCSC.mm10.knownGene)
    tx<-TxDb.Mmusculus.UCSC.mm10.knownGene
    require(BSgenome.Mmusculus.UCSC.mm10)
    gen<-BSgenome.Mmusculus.UCSC.mm10

  }  else if (sp =="rn4") {
    require(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    tx<-TxDb.Rnorvegicus.UCSC.rn4.ensGene
    require(BSgenome.Rnorvegicus.UCSC.rn4)
    gen<-BSgenome.Rnorvegicus.UCSC.rn4

  }

  #read the input data
  coordinates<-read.delim(input_coordinates)

  #create exon data frame
  exonic <- GenomicFeatures::exonicParts(tx, linked.to.single.gene.only=TRUE)
  EX <- data.frame(GeneID=exonic$gene_id,
                   Chr=as.character(GenomicRanges::seqnames(exonic)),
                   Start=GenomicRanges::start(exonic),
                   End=GenomicRanges::end(exonic),
                   Strand=GenomicRanges::strand(exonic),
                   stringsAsFactors=FALSE,
                   seq= getSeq(gen, exonic))

  #create the required merging columns
  coordinates$exon_start<-paste0(coordinates$Chr,":",coordinates$Start)
  EX$exon_start<-paste0(EX$Chr,":",EX$Start)
  coordinates$exon_end<-paste0(coordinates$Chr,":",coordinates$End)
  EX$exon_end<-paste0(EX$Chr,":",EX$End)

  coordinates$length<-abs(as.numeric(coordinates$Start)-as.numeric(coordinates$End))

  #extract sequences
  #if strand is not in the inforamtion then we extract it
  if("Strand" %!in% names(coordinates) | "strand" %!in% names(coordinates)){

    coordinates$seq="NA"

    for (i in 1:nrow(coordinates)){

      print(paste0("Generating sequence for circ ",i))

      if(length(unique(EX$exon_end %in% coordinates$exon_end[i])) == 2 & length(unique(EX$exon_start %in% coordinates$exon_start[i])) == 2 ){

        coordinates$Strand[i]<-EX$Strand[ which(EX$exon_start %in% coordinates$exon_start[i])]
      }

      if(length(unique(EX$exon_end %in% coordinates$exon_end[i])) == 2 & length(unique(EX$exon_start %in% coordinates$exon_start[i])) == 2 ){

        coordinates$seq[i] <- paste(EX$seq[ which(EX$exon_start %in% coordinates$exon_start[i]) : which(EX$exon_end %in% coordinates$exon_end[i])],collapse  = "")
      }
    }
  } else {

    for (i in 1:nrow(coordinates)){

      print(paste0("Generating sequence for circ ",i))

      if(length(unique(EX$exon_end %in% coordinates$exon_end[i])) == 2 & length(unique(EX$exon_start %in% coordinates$exon_start[i])) == 2 ){

        coordinates$seq[i]<-paste(EX$seq[ which(EX$exon_start %in% coordinates$exon_start[i]) : which(EX$exon_end %in% coordinates$exon_end[i])],collapse  = "")
      }
    }

  }

  #create primers

  #coordinates$seq_revcomp<-as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$seq)))

  if( any(is.na(coordinates$seq)) ) {warning('Some of the coorrdinates given are not annotated as exons. Those will be removed')}

  coordinates <- coordinates[coordinates$seq %!in% "NA",]
  #create junction
  coordinates$junction <- paste0(substr(coordinates$seq,nchar(coordinates$seq)-9,nchar(coordinates$seq)),substr(coordinates$seq,1,11))

  #create rev comp junction
  coordinates$junction_revcomp <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction)))

  #top and bottom strand of oligo
  coordinates$TopStrand <- paste0("ctagcagt",coordinates$junction,"tagttatattcaagcata",coordinates$junction_revcomp,"gcg")
  coordinates$BotStrand <- paste0("aattcgc",coordinates$junction,"tatgcttgaatataacta",coordinates$junction_revcomp,"actg")

  if(shift){

    #create junction
    coordinates$junction5p <- paste0(substr(coordinates$seq,nchar(coordinates$seq)-6,nchar(coordinates$seq)),substr(coordinates$seq,1,14))

    #create rev comp junction
    coordinates$junction_revcomp5p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction5p)))

    #top and bottom strand of oligo
    coordinates$TopStrand5p <- paste0("ctagcagt",coordinates$junction5p,"tagttatattcaagcata",coordinates$junction_revcomp5p,"gcg")
    coordinates$BotStrand5p <- paste0("aattcgc",coordinates$junction5p,"tatgcttgaatataacta",coordinates$junction_revcomp5p,"actg")


    #create junction
    coordinates$junction3p <- paste0(substr(coordinates$seq,nchar(coordinates$seq)-12,nchar(coordinates$seq)),substr(coordinates$seq,1,8))

    #create rev comp junction
    coordinates$junction_revcomp3p <- as.vector(Biostrings::reverseComplement(Biostrings::DNAStringSet(coordinates$junction3p)))

    #top and bottom strand of oligo
    coordinates$TopStrand3p <- paste0("ctagcagt",coordinates$junction3p,"tagttatattcaagcata",coordinates$junction_revcomp3p,"gcg")
    coordinates$BotStrand3p <- paste0("aattcgc",coordinates$junction3p,"tatgcttgaatataacta",coordinates$junction_revcomp3p,"actg")

  }

  if(circRNAseq==F){
    coordinates <- coordinates[,-c("seq")]
  }

  if(writetab==T){
    write.table(coordinates,quote = F,sep = "\t",file = output,row.names = F)
  } else {
    return(coordinates) }

}




