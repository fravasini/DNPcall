# DNPcall
R and samtools based script to call Double Nucleotide Polymorphisms (DNPs) based on reads name.

DNPcall leverages the read name information in samtools pileup file format to finely reconstruct the microhaplotype represented by the DNP.

For issues, please email me: francesco.ravasini@uniroma1.it

## Software required
```samtools``` version 1.18 or higher

```R``` version 4.2.0 or higher

```R``` packages: ```tidyverse``` and ```stringr```

## Usage
```Rscript DNPcall.R --bamlist=list_of_bams --DNPs=list_of_DNPs --reference=reference_genome```

**Files:**
```list_of_bams``` is the list of bamfiles to be analyzed, one per row with the full path to them. An example can be found in ```DNPcall/examples/inputs/bamlist.txt```

```list_of_DNPs``` is the list of DNPs to be called. An example can be found in ```DNPcall/examples/inputs/gnomAD_chr22_DNPs_list.txt```. The columns are in the order: chromosome, position #1 of DNP, position #2 of  DNP, reference allele, alternative allele.

```reference_genome``` is the reference genome in fasta format.

## Outputs
```AllGenotypes.txt``` a file containing all the genotypes for each sample and each DNP analyzed. 0/0 is reference homozygous, 0/1 isheterozygous and 1/1 is alternative homozygous. An example van be found in ```DNPcall/examples/outputs/AllGenotypes.txt```

For each sample, a  file named ```full_output_{name_of_the_sample}``` is also produced. An example van be found in ```DNPcall/examples/outputs/full_output_S_Abkhasian-1.txt```. The columns represent:

```chr```: chromosome

```pos1```: position #1 of DNP

```pos2```: position #2 of DNP

```REF```: reference allele

```ALT```: alternative allele

```N_REF```: number of reads matching reference allele

```N_ALT```: number of reads matching alternative allele

```N_Other```: number of reads matching other alelle, different from both reference and alternative

```allele_ratio```: ratio between reference and alternative allele, in the form N_REF/(N_REF+N_ALT)

```combined_bases```: list of all the combination of bases found in the two position of the DNP (including reference and alternative alleles)

```combined_counts```: number of all the reads matching all the combination of bases, in the same order ase the previous column

```total_count```: total number of reads

## Attention

The allele ratio to call the genotypes are set to: ≥ 0.9 for homozygous reference; ≤ 0.1 for alternative homozygous; 0.4 ≤ x ≤ 0.6 for heterozygous. If these parameters are too strict or too loose for your dataset, change these values in line 242-244 of the ```DNPcall.R``` script.

 

