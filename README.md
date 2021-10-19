resolve-nested-genotypes
=====

`vg deconstruct` ([from vg](https://github.com/vgteam/vg)) uses VCF info fields to describe how sites are nested in the original pangenome graph:

```
##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)
##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph
##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl
```

Given genotypes for a subset of sites from such a VCF, this tool uses these INFO fields to propagate the genotypes down to child sites of the genotyped sites.  For example, consider these two sites from the deconstructed VCF, where the second is nested inside the first. 

```
GRCh38.chr20	139949	>73154428>73154433	AAAGTC	A,AAAGAC	60	.	AC=3,1;AF=0.0337079,0.011236;AN=89;AT=>73154428>73154429>73154430>73154432>73154433,>73154428>73154433,>73154428>73154429<73154431>73154432>73154433;LV=0;
GRCh38.chr20	139953	>73154429>73154432	T	A	60	.	AC=1;AF=0.0116279;AN=86;AT=>73154429>73154430>73154432,>73154429<73154431>73154432;LV=1;PS=>73154428>73154433
```

The allele traversals are:

Parent
```
allele 0: >73154428>73154429>73154430>73154432>73154433
allele 1: >73154428>73154433
allele 2: >73154428>73154429<73154431>73154432>73154433
```

Child
```
allele 0: >73154429>73154430>73154432 (substring of parent allele 0)
allele 1: >73154429<73154431>73154432 (substring of parent allele 2)
```

This gives us the following mapping of all possible parent genotypes to child genotypes:
```
0/0 -> 0/0
0/1 -> 0/.
0/2 -> 0/1
1/1 -> ./.
1/2 -> ./1
2/2 -> 1/1
```

The intended use case is as follows
* extract all top-level sites from a deconstructed VCF (ex zcat `deconstruct.vcf.gz | grep -v PS | bgzip`)
* genotype the top-level sites with [Pangenie](https://bitbucket.org/jana_ebler/pangenie/src/master/)
* infer genotypes of all nested sites from the top-level sites

#### ID INFO Field

The output VCF will be annotated with this field (note that top-level bubbles are whatever's genotyped in the input and need not be LV=0)

##### (A) Top level bubbles (for HPRC: LV=0)
Only LV=0 variants
Will be used as input for PanGenie (PanGenie keeps IDs in INFO field, see below)

chr1    811390  .       AATACAATGTTCTCCAGAAATAGTGTCGCAGAAATAAAGACAGAATATGGAAGGAGACAGGATTGCTGACTAATCTCATATGTACAGGGAGAACGACACTCCATGCAGCTTAGCTCCAACATGTGAGCGTGGGACAGAAAAGCAAGGTGAACCTAAAGACATCCCATGGACACACTGAGCTGGAACCAACTCTGCCCATAGGTGGTGGGGCCAGGTTCAGATGCTTGCAGTAACCGCATCCTCCTGTGAATCTGGGTCTACATGGCTGTGGCTCCCTGGGGGCTACCCACACTGTGGATTCCACAGGTGGAGGTGACCTTCTTGGGCACTGCTGGGACACTTCCTGGGCCTGCCTGGGCAGTGGTAGGGCTCAGAGAACTTCAGCGTTAGGGCCTTGGGGAGTCTATGTCAGACTTGAGGACAATACAATGTTCTCCAGAAATAGTGTCGCAGAAATAAAGACAGAATATGGAAGGAGACAGGATTGCTGACTAATCTCATATGTACAGGGAGAACGACACTCCATGCAGCTTAGCTCCAACATGTGAGTGTGGGACAGAAAAGCAAGGTGAACCTAAAGACATCCCATGGACACACTGAGCTGGAACCAACTCTGCCCATAGGTGGTGGGGCCAGGTTCAGATGCTTGCAGTAACCGCATCCTCCTGTGAATCTGGGTCTACATGGCTGTGGCTCCCTGGGGGCTACCCACACTGTGGATTCCACAGGTGGAGGTGACCTTCTTGGGCACTGCTGGGACACTTCTTGGGCCTGCCTGGGCAGTGGTAGGGCTCAGAGAACTTCAGCGTTAGGGCCTTGGGGAGTCTATGTCAGACTTGAGGAC,A,AATACAATGTTCTCCAGAAATAGTGTCGCAGAAATAAAGACAGAATATGGAAGGAGACAGGATTGCTGACTAATCTCATATGTACAGGGAGAACGACACTCCATGCAGCTTAGCTCCAACATGTGAGCGTGGGACAGAAAAGCAAGGTGAACCTAAAGACATCCCATGGACACACTGAGCTGGAACCAACTCTGCCCATAGGTGGTGGGGCCAGGTTCAGATGCTTGCAGTAACCGCATCCTCCTGTGAATCTGGGTCTACATGGCTGTGTCTCCCTGGGGGCTACCCACACTGTGGATTCCACAGGTGGAGGTGACCTTCTTGGGCACTGCTGGGACACTTCCTGGGCCTGCCTGGGCAGTGGTAGGGCTCAGAGAACTTCAGCGTTAGGGCCTTGGGGAGTCTATGTCAGACTTGAGGAC,AATACAATGTTCTCCAGAAATAGTGTCGCAGAAATAAAGACAGAATATGGAAGGAGACAGGATTGCTGACTAATCTCATATGTACAGGGAGAACGACACTCCATGCAGCTTAGCTCCAACATGTGAGCGTGGGACAGAAAAGCAAGGTGAACCTAAAGACATCCCATGGACACACTGAGCTGGAACCAACTCTGCCCATAGATGGTGGGGCCAGGTTCAGATGCTTGCAGTAACCGCATCCTCCTGTGAATCTGGGTCTACATGGCTGTGGCTCCCTGGGGGCTACCCACACTGTGGATTCCACAGGTGGAGGTGACCTTCTTGGGCACTGCTGGGACACTTCCTGGGCCTGCCTGGGCAGTGGTAGGGCTCAGAGAACTTCAGCGTTAGGGCCTTGGGGAGTCTATGTCAGACTTGAGGAC  .       PASS    ID=chr1-811517-SNV-C-T:chr1-811733-SNV-C-T,chr1-811391-DEL-0-421,chr1-811660-SNV-G-T,chr1-811591-SNV-G-A


The ID field contains one ID per alternative allele. If an allele contains nested alleles, its ID is composed of a sequence of IDs of the lowest level alleles, separated by “:” (see first ALT allele in example above)

The IDs itself are unique, for HGSVC we defined them based on the following pattern: 
<chrom>-<allele start>-<type>-<REF>-<ALT> (for SNV only) 
<chrom>-<allele start>-<type>-<count>-<allele length>  (for other variant types: INS, DEL, COMPLEX)

Note: for variants != SNVs, the <count> field is needed to make sure IDs are unique in case there are multiple alleles with the same length at the same location.

##### (B) Equivalent decomposed representation:

Biallelic (separate record for each alternative allele)
If an alternative allele is composed of multiple IDs (see above), there is a record for each corresponding alternative allele (see below)

chr1    811390  .     AATACAATGTTCTCCAGAAATAGTGTCGCAGAAATAAAGACAGAATATGGAAGGAGACAGGATTGCTGACTAATCTCATATGTACAGGGAGAACGACACTCCATGCAGCTTAGCTCCAACATGTGAGCGTGGGACAGAAAAGCAAGGTGAACCTAAAGACATCCCATGGACACACTGAGCTGGAACCAACTCTGCCCATAGGTGGTGGGGCCAGGTTCAGATGCTTGCAGTAACCGCATCCTCCTGTGAATCTGGGTCTACATGGCTGTGGCTCCCTGGGGGCTACCCACACTGTGGATTCCACAGGTGGAGGTGACCTTCTTGGGCACTGCTGGGACACTTCCTGGGCCTGCCTGGGCAGTGGTAGGGCTCAGAGAACTTCAGCGTTAGGGCCTTGGGGAGTCTATGTCAGACTTGAGGAC  A       .       .       ID=chr1-811391-DEL-0-421       
chr1    811517  .     C       T       .       .       ID=chr1-811517-SNV-C-T
chr1    811591  .     G       A       .       .       ID=chr1-811591-SNV-G-A
chr1    811660  .     G       T       .       .       ID=chr1-811660-SNV-G-T
chr1    811733  .     C       T       .       .       ID=chr1-811733-SNV-C-T

This biallelic representation contains the same information, just in a different representation. It can be used for evaluation and comparison to external truth sets.


### Installing

Download 
```
git clone https://github.com/glennhickey/resolve-nested-genotypes.git
cd resolve-nested-genotypes
```

Build with Rust
```
cargo build --release 
```

### Running

`target/release/resolve-nested-genotypes <deconstructed VCF> <genotyped VCF> > out.vcf`

Note that sites in the `<genotyped VCF>` must have the exact same alleles in the exact same order as sites in the deconstructed VCF.  Only the coordinates are used to map between the two files, so inconsistencies may not be caught.  `<genotyped VCF>` does not need the same IDs or INFO fields. 

### About

Author: Glenn Hickey

Licence: [MIT](LICENCE)

