# afplot

This is a tool to plot allele frequencies in VCF files. 

There are three plot modes:

* histogram: this will make a histogram with kernel density plot for every chromosome.
* scatter: This will create a scatter plot of allele frequencies per chromosome, with the position on the chromosome on the x-axis
* distance: This will create a scatter plot of the *distance* to the expected theoretical allele frequencies.
 
Multiple VCF files can be supplied simultaneously.
When only a single VCF file is supplied, plots will be colored on call type.
When multiple VCF files are supplied, plots will be colored on label per VCF file. 

Only one sample per VCF file can be plotted. 

It currently assumes the presence of an `AD` column in the `FORMAT` field. 
This column should contain the depth per allele, with the reference allele being first.
 
All VCFs should be indexed with tabix, and should contain contigs in the header.


## Requirements

* Python 3.4+
* numpy
* matplotlib
* pandas
* seaborn
* progressbar2
* pysam
* pyvcf

## Usage


```
usage: afplot [-h] -v VCF -l LABEL [-s SAMPLE] -o OUTPUT [--dpi DPI] [-k]
                 (--scatter | --histogram | --distance) [-e EXCLUDE_PATTERN]

    Create scatter plots or histogram of allele frequencies in vcf files.
    If only one VCF is supplied, plots will be colored on call type (het/hom_ref/hom_alt).
    If multiple VCF files are supplied, plots will be colored per file/label.
    Only *one* sample per VCF file can be plotted.

    Your VCF file *MUST* contain an AD column in the FORMAT field.
    Your VCF file *MUST* have contig names and lengths placed in the header.
    Your VCF file *MUST* be indexed with tabix.

    VCF files preferably have the same contigs,
    i.e. they are produced with the same reference.
    If this is not the case, this script will select the vcf file with the largest number of contigs.

    You may exclude contigs by supplying a regex pattern to the -e parameter.
    This parameter may be repeated.
    

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     Input vcf file(s)
  -l LABEL, --label LABEL
                        Labels to vcf file(s)
  -s SAMPLE, --sample SAMPLE
                        Sample identifiers (1 per vcf). Uses first sample in
                        vcf by default
  -o OUTPUT, --output OUTPUT
                        Path to output png
  --dpi DPI             DPI for output png (default: 300)
  -k, --kde-only        Only show kernel density plot on histogram
  --scatter             Make scatter plot of AFs per chromosome
  --histogram           Make histogram of AFs per chromosome
  --distance            Create scatter plot of distances to expected AFs
  -e EXCLUDE_PATTERN, --exclude-pattern EXCLUDE_PATTERN
                        Regex pattern to exclude from contig list
```


## Examples

### Single VCF

* `afplot -v my.vcf.gz -l my_label -s my_sample --histogram -o mysample.histogram.png`

### Multiple VCFs

* `afplot -v my1.vcf.gz -l my_label1 -s my_sample1 -v my2.vcf.gz -l my_label2 -s my_sample2 --histogram -o both_samples.histogram.png` 

Grouping samples can be achieved by supplying identical labels to samples. E.g.

* `afplot -v 1.vcf.gz -v 2.vcf.gz -v 3.vcf.gz -v 4.vcf.gz -l group1 -l group1 -l group2 -l group2 [...] `

### Excluding contigs

In certain cases, you might not want to plot all contigs.
For instance, when your vcf header contains many small unplaced contigs. 
This can be achieved by supplying a regex pattern to the `-e` flag.
For instance, all contigs containing "gl" can be filtered out by doing:

* `afplot [...] -e '.*gl.*' `

## License

MIT
