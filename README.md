# CRISPR4P


[![DOI](https://zenodo.org/badge/45244871.svg)](https://zenodo.org/badge/latestdoi/45244871)

##Introduction

We have created a tool (CRISPR4P) that aids *S. pombe* researchers in the generation of deletion mutants using CRISPR/Cas9 technology. 

The main steps of the CRISPR/Cas9-based method to generate gene deletions are briefly highlighted below. A more detailed methodology is available in the manuscript: https://wellcomeopenresearch.org/articles/1-19/v1 

1. Identify suitable sgRNAs to target region of interest using CRISPR4P tool. 
2. Design of primers required for whole process using CRISPR4P: 1) sgRNA cloning; 2) synthesis of DNA template for homologous recombination (HR template) for gene deletion; and 3) checking primers to confirm gene deletion.
3. Clone sgRNAs into nourseothricin-selectable plasmid pMZ379 that contains Cas9 enzyme gene, the natMX6 selection marker and the rrk1 promoter/leader.
4. Generate HR template by PCR using primers with sequences flanking the region of interest and overlapping at their 3’ ends.
5. Delete region of interest by co-transforming sgRNA/Cas9-plasmid and HR template into S. pombe cells that have been synchronized and cryopreserved to increase transformation efficiency.
6. Select the smallest colonies from selective plate  and check these colonies for deletion junction by colony PCR .

Note that this approach can be adapted for applications other than gene deletions, such as insertion of point mutations or tags. The CRISPR4P tool allows the user to identify possible sgRNA sequences in any region of interest for other applications of the CRISPR/Cas9-based genome editing.

##CRISPR4P primer design tool
Available primer design programs for gene targeting in *S. pombe* allow the manipulation of coding sequences using the standard PCR-based method , or rely on current gene annotations to generate a database that contains primers for deletion of non-coding RNAs, 3’-UTRs or tRNAs. We have designed an online tool, written in Python 2.7 ( www.python.org/), to help with the design of all the different primers required for CRISPR/Cas9-based deletion of virtually any region in the S. pombe genome. This tool, named CRISPR4P (CRISPR ‘for’ Pombe or CRISPR Pombe PCR Primer Program), is freely available from our website ( http://bahlerweb.cs.ucl.ac.uk/cgi-bin/crispr4p/webapp.py). CRISPR4P designs PCR primers for sgRNA cloning and primers to generate the HR template, and also checks primers to verify gene deletions. 

####How to cite this article:
Rodríguez-López M, Cotobal C, Fernández-Sánchez O et al. A CRISPR/Cas9-based method and primer design tool for seamless genome editing in fission yeast [version 1; referees: 1 approved, 1 approved with reservations]. Wellcome Open Res 2016, 1:19
(doi: 10.12688/wellcomeopenres.10038.1) 


