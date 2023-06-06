# WGSPooledSeq

A compiled R apackage to faciliate identification of candidate genes from pooled sequencing WGS data in C. elegans.

## Description

Ubituitiously overexpressing a GATA-type transcription factor ELT-7 in C. elegans caused pharynx, and somatic gonad transdifferentiates (Td) into 
intestine-like tissue. In addition to Td, animals also arrested on early stages developmental arrest phenotype plausible due to the 
pharyngeal-to-intestine Td. To identify genes involved in the Td and developmental arrest caused by ELT-7 overexpression, 12.5 million 
EMS-mutagenized ELT-7 transgene animals were subjected to the viability selection after ELT-7 overexpression and 660 mutant animals escaped 
ELT-7-mediated developmental arrest were isolated and established into independent lines. The genomic DNA of 660 mutant lines was further pooled 
into six groups and sequenced.

The whole-genomic-sequencing data were then trimmed and mapped to the N2 reference (WBcel235). Finally, variant calling was done with the CRISP 
program [1], and the effect of every SNV was evaluated with Ensembl Variant Effect Predictor (VEP) [2] and SIFT [3].

The WGSPooledSeq R package analyzed the VEP and SIFT data to provide a summarized table of mutation density and the deleterious ratio of coding 
genes. In addition, WGSPooledSeq calculates the mutation distribution with Kolmogorov–Smirnov test and returns a p-value indicating a based (p<.05) 
or uniform distribution.

## Getting Started

### Dependencies

* dplyr
* tidyr
* ggplot2
* ggrepel
* patchwork
* The .vcf file created by CRISP and convert to .txt file.
* The .txt file created by VEP.
* The .csv file created by SIFT4G.

### Installing

* install_github("THYeh44/WGSPooledSeq")

### Executing program

* library(WGSPooledSeq)

* Remove high allele frequency SNVs
```
rmAF(CRISP.vcf.txt, Output.txt, AF threshold (0.0-1.0))
```
* Generate a summary table
```
pooled_analysis(VEP.txt, SIFTannotations.csv, Remove duplicated region(T/F), Output.txt)
```
* Generate a mutation density figure
```
mut_fig(pooled_analysis_output.txt, Chromosome, Left_max, Left_min, Left_ threshold, Right_max, Right_min, Right_ threshold)
```

## Authors

Tsunghan Yeh (tsunghan_yeh@.ucsb.edu)


## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the [GPL-3] License.

## References
* [1] A statistical method for the detection of variants from next-generation resequencing of DNA pools, V. Bansal. Bioinformatics 26(12), 2010 (https://github.com/vibansal/crisp)
* [2] Ensembl Variant Effect Predictor (VEP) (https://useast.ensembl.org/info/docs/tools/vep/index.html)
* [3] SIFT web server: predicting effects of amino acid substitutions on proteins, NL Sim et al., Nucleic Acids Res, 2012 (https://sift.bii.a-star.edu.sg/index.html)
