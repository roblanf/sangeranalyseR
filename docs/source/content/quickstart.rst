Quick Start Guide
=================

This page provides simple quick-start information for using sangeranalyseR with :code:`.ab1` files. Please read the :ref:`Beginners Guide` page for more details on each step.

If you haven't already, please follow the steps in the :ref:`Installation` page to install and load sangeranalyseR.

Super-Quick Start (3 lines of code)
+++++++++++++++++++++++++++++++++++

The most minimal example gets the job done in three lines of code. More details below.

.. code-block:: R

   my_aligned_contigs <- SangerAlignment(parentDirectory     = "./my_data/",
                                         suffixForwardRegExp = "_[0-9]\*_F.ab1$",
                                         suffixReverseRegExp = "_[0-9]\*_F.ab1$")

   writeFasta(my_aligned_contigs)

   generateReport(my_aligned_contigs)


Step 1: Prepare your input files
++++++++++++++++++++++++++++++++

Put all your :code:`AB1` files in a directory :code:`./my_data/`. The directory can be called anything.

Name your files according to the convention :code:`contig_index_direction.ab1`. E.g. :code:`Drosophila_COI_1_F.ab1` and :code:`Drosophila_COI_2_R.ab1` describes a forward and reverse read to assemble into one contig. You can have as many files and contigs as you like in one directory.

Step 2: Load and analyse your data
++++++++++++++++++++++++++++++++++

.. code-block:: R

   my_aligned_contigs <- SangerAlignment(parentDirectory     = "./my_data/",
                                         suffixForwardRegExp = "_[0-9]\*_F.ab1$",
                                         suffixReverseRegExp = "_[0-9]\*_F.ab1$")


This command loads, trims, builds contigs, and aligns contigs. All of these are done with sensible default values, which can be changed. I


Step 3 (optional): Explore your data
++++++++++++++++++++++++++++++++++++

.. code-block:: R

   launchApp(my_aligned_contigs)

This launches an interactive Shiny app where you can view your analysis, change the default settings, etc.


Step 4: Output your aligned contigs
+++++++++++++++++++++++++++++++++++

.. code-block:: R

   writeFasta(my_aligned_contigs)

This will save your aligned contigs as a FASTA file.

Step 5 (optional): Generate an interactive report
+++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: R

   generateReport(my_aligned_contigs)

This will save a detailed interactive HTML report that you can explore.

|

For more detailed analysis steps, please choose one the following topics :

* :ref:`Beginners Guide`

* :ref:`Advanced User Guide - *SangerRead* (**AB1**)`

* :ref:`Advanced User Guide - *SangerContig* (**AB1**)`

* :ref:`Advanced User Guide - *SangerAlignment* (**AB1**)`

* :ref:`Advanced User Guide - *SangerRead* (**FASTA**)`

* :ref:`Advanced User Guide - *SangerContig* (**FASTA**)`

* :ref:`Advanced User Guide - *SangerAlignment* (**FASTA**)`
