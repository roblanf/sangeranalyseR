sangeranalyseR's tutorial
*************************

Why sangeranalyseR
==================
Sanger sequencing was first proposed in 1977 and is still widely used in sequencing high-quality DNA. There are some widely used tools e.g. `Geneious <https://www.geneious.com>`_, `CodonCode Aligner <https://www.codoncode.com/aligner/>`_ and `Phred-Phrap-Consed <http://www.phrap.org/phredphrapconsed.html>`_; however, these are all commercial software which require expensive license fee. Therefore, we develop sangeranalyseR allowing users to do Sanger sequencing data analysis in pure R environment. sangeranalyseR is an open source software which provides another option for biologists and clinical researchers to do Sanger sequencing data analysis in an easy way.

|

Main features
=============
* **Well-structured S4 classes**: There are three S4 classes in sangeranalyseR which are *SangerRead*, *SangerContig* and *SangerAlignment* representing different levels of Sanger sequencing data analysis and storing user inputs and results in each level. They are the main structure of sangeranalyseR for the following analysis.
* **Automated data analysis**: The S4 instance creation process is automated with systematic named input files by default parameters. The downstream analysis function is maximally simplified.
* **Interactive Shiny apps**: Local Shiny apps for *SangerContig* and *SangerAlignment* are provided to visualize the S4 instance. Users are allowed to change trimming and chromatogram parameters inside Shiny apps.
* **Exporting reads to FASTA**: Aligned results and trimmed reads can be written into **FASTA** file format.
* **Thorough report**: A comprehensive report can be created with a single command.
* Pure R environment: Users can do Sanger sequencing data analysis in pure R environment.

|

User support
============
Please go through the :ref:`Documentation` below first. If you have further questions, feedback, new ideas, and bug reports, please sign up the following Google group and post a topic to the

https://groups.google.com/d/forum/iqtree

|

Documentation
=============

.. toctree::
   :maxdepth: 3

   content/getting_start
   content/quick_commands
   content/beginner
   content/advance_user_SangerRead_ab1
   content/advance_user_SangerContig_ab1
   content/advance_user_SangerAlignment_ab1

   content/advance_user_SangerRead_fasta
   content/advance_user_SangerContig_fasta
   content/advance_user_SangerAlignment_fasta
   content/conclusion
   content/license
   content/contact
   content/help
