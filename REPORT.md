# BINF7001 ASSESSMENT 3 (2024) | De novo genome assembly and phylogenomics

# 1. Genome data and assembly

## 1.1 Sequencing data

My dataset is ERR9872452. It was sequenced using the ILLUMINA NextSeq 500 sequencing platform. It produces 2x 150bp reads, with a total read count of 2197746. The G+C content is approximately 38%.


## 1.2 Genome assembly

I used the velvet program to assembly the genome. Velvet requires a kmer parameter to be selected. I used the jellyfish program to count various kmers and graph them. Values tested were k=11, 17, 21, 27, 31, 51. 

By visual inspection, either k=21, k=27 or k=31 would be good candidates. The sharp peak and long tail of k=51 could capture too many long, infrequent occurring kmers, leading to not enough overlaps. Too short would lead to ambiguity, with too many overlaps that may be incorrect or not meaningful.

I then used velveth to generate the metadata for the alignment, and velvetg to do the assembly.

I ran this pipeline for the proposed kmer values:

| kmer | N50 scaffold length | Maximum scaffold length | Total scaffold length | Total number of scaffolds | Config mean length | % assembled |
|------|------|------|------|------|------|------|
| k21 | 82/9.227 KB | 60.132 KB | 2,514,955 | 1,082 | 99.98% | 2060637 / 2197746 (93.76%) |
| k31 | 36/22.561 KB | 105.043 KB | 2,532,370 | 521 | 99.99% | 2188030 / 2197746 (99.56%) |
| k41 | 391/2.038 KB | 10.011 KB | 2,663,113 | 3,166 | 100.00% | 2160007 / 2197746 (98.28%) |
| k51 | 44/18.254 KB | 54.926 KB | 2,506,012 | 464 | 99.99% | 2165246 / 2197746 (98.52%) |
| k61 | 55/14.643 KB | 45.492 KB | 2,455,163 | 401 | 99.85% | 2134797 / 2197746 (97.14%) |
| k71 | 38/20.195 KB | 65.611 KB | 2,439,224 | 254 | 99.83% | 2106757 / 2197746 (95.86%) |


I settled on k=31. Here's why:

- Has a high level of completeness (99.56%).
- While the number of scaffolds that contribute to 50% of the length is low (36), the 50th largest is 22.563KB. This means we have a smaller number of longer scaffolds - this indicates larger, more complete scaffolds.
- The L50 value is very large (105KB). Although only one scaffold, this metric supports the 36/22.561KB figure, and is further evidence for k=31, the final assembly is not fragmented.

# 2. Ab initio gene prediction

## (a) Brief description of your approach

I used the mat file from the practical. It is for Staphylococcus Aureus; my sample is Staphylococcus Pseudintermedius, which is known to be very genetically similar and have many common proteins. I would like to build my own matrix file and see how that compares, if time permits. https://github.com/kuleshov/nanoscope My command:

```sh
genemark -opn -m /opt/BINF7001/2024/Prac8_2024/Staph_aureus_JKD6008.mat ./k31/contigs.fa
```

According to GeneMark, the GC content is 37.5% - this aligns with what we expected.


## (b) Total number of predicted genes

The total number of predicted genes is 4560, found using the following command:

```sh
cat k31/contigs.fa.orf | grep "^>" | wc -l
```

## (c) Average gene length

The average length of the predicted genes is 302 amino acids. This was derived using python from the `protein.fa.orf`, which the proteins as extracted from the GeneMark `configs.fa.orf`.

## (d) Length and function of the longest gene

Longest gene is 1571 amino acids. I ran a BLAST on the protein on NCBI and the best matches suggest this protein is either 

- LPXTG-anchored putative endo-alpha-N-acetylgalactosaminidase SpsG
- YSIRK-type signal peptide-containing protein

Both have been suggested to help secure surface proteins to the cell wall in gram-positive bacteria (Bae, 2003).


# 3. Phylogenomic analysis (maximum two pages; 8)

## (a) Brief description of your approach (name(s) of program(s), key parameters used) 

After using GeneMark to predict genes, I create three blast databases using the `makeblastdb` program. One for the contigs (`contigs.fa`, one for the genes (nucleotides, `nt.fa`) and one for the proteins (`protein.fa`). 

```sh
makeblastdb -dbtype nucl -in inputs/contigs.fa
makeblastdb -dbtype nucl -in inputs/nc.fa
makeblastdb -dbtype prot -in inputs/protein.fa
```

I then used `blastn` query my `nt` and `contigs` databses using the 16s rRNA as the query. The query against `contigs.fa` found 2 results; the query against `nt.fa` provided none. This means two contigs sequences have similarity to the 16s rRNA query.

```sh
echo "⚡️ Query inputs/contigs.fa"
blastn -query gene_query.fa -db inputs/contigs.fa -outfmt 6 -evalue 1e-10

echo "🚀 Query inputs/nc.fa"
blastn -query gene_query.fa -db inputs/nc.fa -outfmt 6 -evalue 1e-10
```

Here's the two hits, shortened for brevity:

```
⚡️ Query inputs/contigs.fa
16S_rRNA	NODE_12_length_1715_cov_724.954529	94.565	1748	73	17	121	1849	1745	1	0.0	2682
16S_rRNA	NODE_977_length_54_cov_585.851868	91.935	62	0	1	1934	1995	1	57	7.97e-16	82.4
🚀 Query inputs/nc.fa
None!
```

The two hits are 94.5% and 91.9% similar representing a similar match. The length is 1748 for the first candidate, and 62 for the second. The first hit is substantially longer, and represents a more complete match.

I extracted it with:

```sh
samtools faidx -i inputs/contigs.fa NODE_12_length_1715_cov_724.954529:1-1745
```

## (b) Total number of homologous protein groups, and number of single-copy groups 

I created a new directory, `my_proteins` and softlinked the 7 protein databases. I included my own (`protein.fa`). I ran `orthofinder -og -f my_proteins`.

The console output after running says "OrthoFinder assigned 17605 genes (92.2% of total) to 2795 orthogroups". There are 1008 single-copy genes.

## (c) 16S rRNA gene tree (8 taxa, rooted using outgroup)

We can infer a tree using MUSCLE and convert to `nex` for usage with MyBayes:

```sh
muscle -in all_16s.fa > 16s_rRNA_gene_tree.aln
readseq -a -f17 16s_rRNA_gene_tree.aln > 16s_rRNA_gene_tree.nex
```

## (d) protein tree of the chosen housekeeping gene (8 taxa, rooted using outgroup)

I chose dnab. It is

```fasta
>tr|A0A7T7NY26|A0A7T7NY26_STAPS Replicative DNA helicase OS=Staphylococcus pseudintermedius OX=283734 GN=dnaB PE=3 SV=1
MDEMYEHNRMPHSHEAEQSVLGAIFLDPELMSSTQEILLPESFYRGAHQHIFRAMMDLNE
DGKDIDIVTVLDRLTQEGVVNEAGGPQYLAEITSNVPTTRNIQYYTDVVFKNAVKRKLIH
TADSIANDGYNDELDLDTVLNDAERRILELSSTRESDGFKDIRDVLGQVYDNAEQLDQNS
GQTPGIPTGYRDLDQMTAGFNRNDLIILAARPSVGKTAFALNIAQKVATHEDQYTVGIFS
LEMGADQLATRMICSSGNVDSNRLRTGTMTEEDWNRFTVAVGKLSRTKIFIDDTPGVRIT
DIRSKCRRLKQEHGLDMIVIDYLQLIQGSGSRASDNRQQEVSEISRMLKAIARELECPVI
ALSQLSRGVEQRQDKRPMMSDIRESGSIEQDADIVAFLYRDDYYNRGDGDDDDDDGGFEP
QTNDENGEIEIIIAKRRNGPTGTVKLHFMKQYNKFTDIDYAHADMG
```

Next we use Orthofinder - the 7 proteins sets provided, and the predicted ones.

```
my_proteins
|-- Macrococcus_caseolyticus.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Macrococcus_caseolyticus.faa
|-- OrthoFinder
|-- Staphylococcus_aureus_subsp_aureus_NCTC8325.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_aureus_subsp_aureus_NCTC8325.faa
|-- Staphylococcus_chromogenes.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_chromogenes.faa
|-- Staphylococcus_epidermidis.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_epidermidis.faa
|-- Staphylococcus_haemolyticus.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_haemolyticus.faa
|-- Staphylococcus_saprophyticus.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_saprophyticus.faa
|-- Staphylococcus_sciuri.faa -> /opt/BINF7001/2024/Prac9_2024/proteins/Staphylococcus_sciuri.faa
`-- protein.fa
```

Output:

```
Writing orthogroups to file
---------------------------
OrthoFinder assigned 17605 genes (92.2% of total) to 2795 orthogroups. Fifty percent of all genes were in orthogroups with 8 or more genes (G50 was 8) and were contained in the largest 1084 orthogroups (O50 was 1084). There were 1194 orthogroups with all species present and 1008 of these consisted entirely of single-copy genes.

2024-10-02 20:58:54 : Done orthogroups

Results:
    /home/s4206685/phylo/my_proteins/OrthoFinder/Results_Oct02/
```

To find the correct orthogroup, I did `blastp` on a few different proteins sets from the 8 I have. I got a few results. All of those had labels that appeared in the same orthogroup, `./Orthogroup_Sequences/OG0000161.fa`.

It has 9 sequences. I was expected 8, since this is supposedly a single copy group. 

Alignment:  

```sh
muscle -in OG0000161.fa > dnaB.aln
readseq -a -f17 dnaB.aln > dnaB_tree.nex
```

## (e) one key difference/similarity between the two trees

## (f) one plausible explanation as to why such a difference occurs


# References

Bae, T., & Schneewind, O. (2003). The YSIRK-G/S motif of staphylococcal protein A and its role in efficiency of signal peptide processing. Journal of bacteriology, 185(9), 2910–2919. https://doi.org/10.1128/JB.185.9.2910-2919.2003
