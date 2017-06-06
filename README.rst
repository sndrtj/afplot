afplot
======

This is a tool to plot allele frequencies in VCF files.

The two main subcommands that are available are: \* ``regions``: Plot
single regions or regions from a bed file, optionally with a margin. \*
``whole-genome``: Create a single image for every chromosome on the
genome.

Both subcommands have three modes:

-  ``histogram``: This will create a histogram with kernel density plot
   of allele frequencies.
-  ``scatter``: Create a scatter plot of allele frequencies, along the
   region or chromosome.
-  ``distance``: Create a scatter plot of distances to theoretical
   allele frequencies, along the region or chromosome. This only makes
   sense for autosomes of diploid organisms.

By default, colors correspond to call type (hom\_alt/ref/hom\_ref).

Multiple VCF files can be supplied simultaneously for the
``whole-genome`` subcommand, in which case they can be grouped by label.
When multiple VCF files are supplied, plots will be colored on label per
VCF file.

Only one sample per VCF file can be plotted.

We currently assume the presence of an ``AD`` column in the ``FORMAT``
field. This column should contain the depth per allele, with the
reference allele being first.

All VCFs should be indexed with tabix, and should contain contigs in the
header.

Installation
------------

afplot is available through pypi with: ``pip install afplot``

Requirements
------------

-  Python 3.4+
-  click
-  numpy
-  matplotlib
-  pandas
-  seaborn
-  progressbar2
-  pysam
-  pyvcf

Usage
-----

.. code:: text

    Usage: afplot [OPTIONS] COMMAND [ARGS]...

      Plot allele frequencies in VCF files.

      Two basic modes exist:
        - regions: Plot histogram, scatter or distance plots per
          user-specified region.
        - whole-genome: Plot histogram, scatter or distance plots over the
          entire genome.

    Options:
      --help  Show this message and exit.

    Commands:
      regions       Region plots
      whole-genome  Whole-genome plots

Examples
--------

Single VCF on a single region
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``afplot regions histogram -v my.vcf.gz -o output_dir -R chr1:100-200``

Single VCF on a bed file
~~~~~~~~~~~~~~~~~~~~~~~~

-  ``afplot regions histogram -v my.vcf.gz -o output_dir -L regions.bed``

Single VCF whole genome
~~~~~~~~~~~~~~~~~~~~~~~

-  ``afplot whole-genome histogram -v my.vcf.gz -l my_label -s my_sample -o mysample.histogram.png``

Multiple VCFs whole genome
~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``afplot whole-genome histogram -v my1.vcf.gz -l my_label1 -s my_sample1 -v my2.vcf.gz -l my_label2 -s my_sample2 -o both_samples.histogram.png``

Grouping samples can be achieved by supplying identical labels to
samples. E.g.

-  ``afplot whole-genome histogram -v 1.vcf.gz -v 2.vcf.gz -v 3.vcf.gz -v 4.vcf.gz -l group1 -l group1 -l group2 -l group2 [...]``

Excluding contigs on whole genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In certain cases, you may not want to plot all contigs. For instance,
when your vcf header contains many small unplaced contigs. This can be
achieved by supplying a regex pattern to the ``-e`` flag. For instance,
all contigs containing "gl" can be filtered out by doing:

-  ``afplot whole-genome [...] -e '.*gl.*'``

Changelog
---------

0.2
~~~

The entire command line interface was changed to use ``click``, instead
of regular argparse. This allows a more complex CLI. In stead of having
flags for plot mode, ``afplot`` now uses subcommands.

While the CLI has changed, and the internals of ``afplot`` have been
refactored, the old-style (version 0.1) API remains in place for now.
This may be deprecated in the future.

Support for plotting regions was added. Region plotting outputs on a
directory, rather than on a single file.

License
-------

MIT
