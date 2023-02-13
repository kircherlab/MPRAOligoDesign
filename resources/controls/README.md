# General controls

Most of them collected by the [UNC MPRA control sequences](https://drive.google.com/drive/u/0/folders/1lRK1ZnC4aldDGwddgjsrhdK6zWrsmHKJ)

All are hg38. Region lengths are 270 (originally derived from 250 bp). There is always a `generate.sh` script how the files where generated.

## noncoding_mendelian_promoter_enhancers


fasta file containing 250 bp sequences centered around the non-coding Mendelian SNVs (enhancer and promoter mutations) supplied by Mike Pazin [here](https://docs.google.com/spreadsheets/d/1c8afZD8OpMWqfY6P86K_qOQdTIl47T8a/edit?usp=sharing&ouid=113585314970533144724&rtpof=true&sd=true)

184 variants (42 enhancer, 142 promoter mutations) resulting in 368 sequences. Original variant positions referenced to hg19, which were lifted over to hg38 using UCSC liftOver tool. Variants where designed on the + strand (not dependent on the gene orientation).  Naming convention used consists of an 11 field (underscore delimited) code:

Chromosome
Start position (1-indexed) of sequence - hg38
End position (1-indexed) of sequence - hg38
Genomic position (1-indexed) of first allele in variant - hg38
Position in sequence of first allele of variant (1-indexed)
Ref
Alt
Indicator if sequence contains ref or alt
Gene (from Pazin spreadsheet)
PMID (from Pazin spreadsheet)
Enhancer or Promoter mutation


## GABA

Here only the onces with ccordiantes are used. Other ones discarded. The original description:

fasta file containing 250 bp sequences derived from those supplied by Chengyu Deng (chengyu.deng@ucsf.edu). Sequences represent MPRA validated active enhancers in four cell types (GABA, Glut, HepG2, and fetal cortical cells). Naming convention for GABAergic and Glutamatergic neuron sequences is cell type, hg38 coordinate, stand, and mean RNA/DNA ratio. For the HepG2 sequences only cell type and mean ratio is supplied. There are two sources of fetal cortical cell data: one sequence specific and one variant-centered. Naming convention for the former is label (da: differentially accessible chromatin regions, positive, negative), genomic coordinates (hg38), mean RNA/DNA ratio, and cell types. ‘cell types’ indicates that the sequence overlaps cell-type-specific scATAC-seq peaks including AstroOligo(astrocyte/oligodendrocyte precursor), dlEN(deep layer excitatory neuron), earlyEN(early excitatory neuron), EndoMural(endothelial and mural cells), IN_CGE(CGE derived interneuron) , IN_MGE(MGE derived interneurons), IPC(intermediate progenitor cell), Microglia, RG(radial glia), ulEN(upper layer excitatory neuron). Naming convention for the variant-centered sequences is rsID, ref/alt allele, genomic coordinate (hg38), and mean RNA/DNA ratio.

Original GABA, glut, and cortical cell sequences (provided in “Oiginal_files” folder) were supplied at 270 bp; those presented here are trimmed to 250 bp removing 10 bp from both ends. Original HepG2 sequences were supplied at 200 bp; those presented here are padded with 25 bp on both ends using random sequences that excluded restriction enzyme recognition sites for MluI, KpnI, XbaI, SpeI, and BsiWI.


## Hon

fasta file containing 250 bp sequences derived from those supplied by Daniel Armendariz in the [Hon lab](https://docs.google.com/spreadsheets/d/18mk-MYRHgZ1dDJB3x8Uby1Pr5G8J2Gmk/edit?usp=sharing&ouid=102313408722438008052&rtpof=true&sd=true). Ten sequences were derived as follows: all coordinates lifted over to hg38; NGFRAP1 region used as supplied trimming 7 bp to get to 250, TBXEnh4 sequence was centered on the supplied region. For the remaining 8 sequences we looked at region overlap with cCREs provided by SCREEN and generated sequences centered on midpoint of cCRE.

TBXEnh1:		EH38E1646128
TBXEnh2:		EH38E1646265
TBXEnh3:		EH38E3042928
TBXEnh5:		Skipped
TBXEnh6:		EH38E1646350 and EH38E3042962 (two sequences)
TBX5_PROM1:	EH38E1646325
TBX5_PROM2:	EH38E1646325 and EH38E1646325 (two sequences)

Naming convention: genomic coordinates (hg38), element (TBXEnh1, e.g.), associated cCRE or “Centered_region” if none.


## AtrialFib

STARR-seq atrial fibrillation allelic effect variants from [“Identification of Function Variant Enhancers Associated with Atrial Fibrillation” Circulation Research (2020) 127:229-243.](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.119.316006)

Variants (24) extracted from Online Table V and lifted over to hg38 (variant data supplied in “Original_files” directory). 250 bp sequences built around variants, with the sequence centered on an ENCODE cCRE (SCREEN v3) when the variant was within the peak or within 100 bp of the center of the region. If centering on cCRE midpoint placed variant within 25 bp of end of sequence, the sequence was shifted towards the midpoint such that the variant remained within 25 bp of end of sequence. Naming convention used: genomic coordinates of sequence (hg38), genomic position of variant, relative position of variant in sequence, ref allele, alt allele, ref/alt indicator which is in sequence, rsID, nearest gene, “STARR-seq-AF”.



## Liang

These are variants in ATAC peaks, that are caQTL, and exhibit significant allele specific chromatin accessibility (ASCA) in heterozygotes, and are eQTL. All of these results are in human neurons (differentiated progeny). Also the QTL are in same direction of effect and the QTL and allelic fold change are in the same direction. Source of variants is [Liang et at 2021](https://pubmed.ncbi.nlm.nih.gov/34017130/). Variant info also provided in “Original_files” directory.

250 bp sequences built around variants, with the sequence centered on an ENCODE cCRE (SCREEN v3) when the variant was within the peak or within 100 bp of the center of the region. If centering on cCRE midpoint placed variant within 25 bp of end of sequence, the sequence was shifted towards the midpoint such that the variant remained within 25 bp of end of sequence. Naming convention used: genomic coordinates of sequence (hg38), genomic position of variant, relative position of variant in sequence, ref allele, alt allele, ref/alt indicator which is in sequence, rsID, “Liang”.

## Engreitz

**NOTE:** All overlap with a TSS. So all get filtered out!

“This is a list of 20 promoters of housekeeping genes that we have used in previous STARR-seq Experiments (see Bergman, Jones et al. Nature 2022) either as promoters (all 20) and/or as enhancers (for 2 of these sequences).  Based on our study, we expect that all of these will act as strong enhancers in MPRA or STARR-seq experiments, and in many/all cell types due to their ubiquitous activity as promoters of housekeeping genes (see paper for further details).  We plan to include these in our future STARR-seq experiments as positive controls for these reasons." [see]](https://docs.google.com/spreadsheets/d/1Yircz6mlGmKiMpUrxZk98-q6TbR_ckQAwyx3nTD__uw/edit#gid=0)

250 bp sequences that represent trimmed versions of the 264 bp sequences that were posted (7 bp removed from both ends). Naming convention used: genomic coordinates of sequence (hg38), gene name, strand, region (hg19; BED format) , promoter strength (poisson), “trim_Engreitz”.

## Selvarajan

These variants (n = 212) were extracted from Selvarajan, et. al. Integrative analysis of liver-specific non-coding regulatory SNPs associated with the risk of coronary artery disease. AJHG (2021) 108. doi: 10.1016/j.ajhg.2021.02.006

These variants were extracted from Table S22 and represent coronary artery disease GWAS variants with allele-specific STARR-seq activity in HepG2 cells. Alleles were correlated to the rsID using Table S6, liftover to hg38 was performed, and allele identities were rectified for three variants undergoing strand flip between builds. Note, in the study variants were tested in common haplotype sequences thus some table entries are for multiple variants. These are presented here individually and sequences generated per variant, not per haplotype combination.

250 bp sequences built around variants, with the sequence centered on either liver ATAC-seq peak (derived from 138 sample study from Mohlke lab, unpublished) or an ENCODE cCRE (SCREEN v3) when the variant was within the peak or within 100 bp of the center of the region. If centering on ATAC summit or cCRE midpoint placed variant within 25 bp of end of sequence, the sequence was shifted such that the variant remained within 25 bp of end of sequence. Naming convention used: genomic coordinates of sequence (hg38), genomic position of variant, relative position of variant in sequence, ref allele, alt allele, ref/alt indicator which is in sequence, rsID, “STARR-seq-HepG2”.


## Mohlke

250 bp sequences built around variants. For hepatocyte variants, the sequence was centered on either liver ATAC-seq peak (derived from 138 sample study from Mohlke lab, unpublished) or an ENCODE cCRE (SCREEN v3) when the variant was within the peak or within 100 bp of the center of the region. For non-hepatocyte variants, the sequence was centered on ENCODE cCRE when the variant was within the peak or within 100 bp of the center of the region. If centering on ATAC summit or cCRE midpoint placed the variant within 25 bp of end of sequence, the sequence was shifted such that the variant remained within 25 bp of end of sequence. Naming convention used: genomic coordinates of sequence (hg38), genomic position of variant, relative position of variant in sequence, ref allele, alt allele, ref/alt indicator which is in sequence, SPDI, “MohlkeHepControls” or “MohlkeNonHepControl”.


## Kircher

Variant-based sequences derived from saturation mutagenesis results from Kircher, et. al. (Kircher M, Xiong C, Martin B, Schubach M, Inoue F, Bell RJA, Costello JF, Shendure J, Ahituv N. Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution. Nature Communication. 2019 Aug 8;10(1):3583. doi: 10.1038/s41467-019-11526-w).

[Source of data](https://github.com/kircherlab/MPRA_SaturationMutagenesis/blob/004b5465f1cf511a6fbd8db858e5504816d7b53d/app.R#L41)

For the LDLR and SORT1 loci, we selected the top 100 variants based on absolute effect size with at least 10 barcodes and p-value < 1e-5, excluding deletions and only considering the top ref/alt if more than one reference allele was listed. Naming convention used: genomic coordinates of sequence (hg38), genomic position of variant, relative position of variant in sequence, ref allele, alt allele, ref/alt indicator which is in sequence, SPDI, “KircherControls”.




## Vista

Elements extracted from a list of cores of elements tested for enhancer activity in mouse embryos provided by Michael Kosicki / Len Pennacchio. The selected elements were selected based on length (<= 300 bp, trimmed or padded to 250), listed as active in liver or forebrain, and annotated indicating they exhibited reliable/high expression. We include a total of 284 sequences, the majority showing activity in the forebrain. Naming convention used: genomic coordinates of sequence (hg38), tested tissues, additional annotation, “vistaElementControl”.

