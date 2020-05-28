sangeranalyseR's tutorial
*************************

Why sangeranalyseR
==================
sangeranalseR is an R package that provides fast, flexible, and reproducible workflows for assembling your sanger seuqencing data into contigs.

It adds to a list of already widely-used tools, like `Geneious <https://www.geneious.com>`_, `CodonCode Aligner <https://www.codoncode.com/aligner/>`_ and `Phred-Phrap-Consed <http://www.phrap.org/phredphrapconsed.html>`_;. What makes it different from these tools is that it's free, it's open source, and it's in R.

|

Main features
=============
* **Pure R environment**: As far as we know, this is the first package that allows end-to-end analysis of Sanger sequencing data in a pure R environment.
* **Automated data analysis**: Given appropriately-named input files, a lot of the data analysis can be automated. Once you've set up an appropriate workflow for your data, you can run it again in seconds.
* **Interactive Shiny apps**: Local Shiny apps mean you visualize the data at many levels, view chromatograms, and adjust things like trimming parameters.
* **Exporting and importing FASTA files**: sangeranalyseR is primarily designed with loading raw :code:`ab1` files in mind, but it can also load sequencesin **FASTA** format. Aligned results and trimmed reads can be written into **FASTA** file format.
* **Thorough report**: A single command creates a comprehensive interactive HTML report that provides a huge amount of detail on the analysis.

|

What sangeranalyseR **doesn't** do
==================================

One really important feature that sangeranalyseR doesn't have is the ability to edit bases by hand. R is just not the right language for this. If you need to edit your reads by hand, we suggest doing that in another tool like `Geneious <https://www.geneious.com>`_, then exporting your reads as FASTA files and following the instructions for using sangeranalyseR with FASTA input.

|

User Manual
===========
If you are already familiar with sangeranalyseR and want to have a quick look at function signatures, please refer to `sangeranalyseR user manual <https://bioconductor.org/packages/devel/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_

|

User support
============
Please go through the :ref:`Documentation` below first. If you have questions about using the package, a bug report, or a feature request, please use the GitHub issue tracker here:

https://github.com/roblanf/sangeranalyseR/issues

|

Key contributors
================

The first (and not very good) version of the package was written by Rob Lanfear (at ANU in Australia), in collaboration with Kirston Barton and Sarah Palmer (then both at the University of Sydney). The second and far far better version of the package was written by Kuan-Hao (Howard) Chao at ANU. (This section was written by Rob Lanfear, lest you think Howard wrote it!)

|

Documentation
=============

.. toctree::
   :maxdepth: 3

   content/installation
   content/quickstart
   content/beginner
   content/advance_user_SangerRead_ab1
   content/advance_user_SangerContig_ab1
   content/advance_user_SangerAlignment_ab1

   content/advance_user_SangerRead_fasta
   content/advance_user_SangerContig_fasta
   content/advance_user_SangerAlignment_fasta
   content/how_to_page
   content/conclusion
   content/license
   content/contact
   content/help
