resolve-nested-genotypes
=====

`vg deconstruct` ([from vg](https://github.com/vgteam/vg)) uses VCF info fields to describe how sites are nested in the original pangenome graph:

```
##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)
##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph
##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl
```

Given genotype for a subset of sites from such a VCF, this tool uses these INFO fields to propagate the genotypes down to child sites of the genotyped sites.  For example, consider these two sites from the deconstructed VCF, where the second is nested inside the first. 

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

### Installing

Download 
```
git clone https://github.com/glennhickey/resolve-nested-genotypes.git
cd resolve-nested-genotypes
```

Build with Rust
```
cargo build
```

### Running

`target/debug/resolve-nested-genotypes <deconstructed VCF> <genotyped VCF>`

Note that sites in the `<genotyped VCF>` must have the exact same alleles in the exact same order as sites in the deconstructed VCF.  Only the coordinates are used to map between the two files, so inconsistencies may not be caught.  `<genotyped VCF>` does not need the same IDs or INFO fields. 

### About

Author: Glenn Hickey
Licence: [MIT](LICENCE)

