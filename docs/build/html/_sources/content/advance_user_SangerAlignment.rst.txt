Advanced User Guide - *SangerAlignment*
=======================================



*SangerAlignment* is the highest class level in sangeranalyseR showed in :ref:`Figure_1<SangerAlignment_hierachy>`. It contains *SangerContig* list and the contigs alignment result. Users can access to *SangerContig* and *SangerRead* instance inside *SangerAlignment* instance. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerAlignment* level.

.. _SangerAlignment_hierachy:
.. figure::  ../image/SangerAlignment_hierachy.png
   :align:   center
   :scale:   20 %

   Figure 1. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

|

Input files preparation
-----------------------
Users can choose to input **ab1** or **FASTA** as their input file format.

ab1 files
+++++++++
The main input file format to create *SangerAlignment* instance is **ab1**. Before starting the analysis, users needs to prepare a directory which contains all the **ab1** files. Here are some filename regulations:

.. note::

    *  All the input files must have **ab1** as its file extension
    *  The reads that belong to the same contig must have the same contig name in its filename.
    *  Forward or reverse direction also needs to be specified in the filename.


There are three parameters, :code:`parentDirectory`, :code:`suffixForwardRegExp`, and :code:`suffixReverseRegExp`, that users need to provide so that program can automatically group all **ab1** files.

.. note::

  * :code:`parentDirectory` is the root directory that contains all the **ab1** files. It can be absolute or relative path. We suggest users to put only target **ab1** files inside this directory without other unrelated files.
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all filename of forward reads. :code:`grepl` function in R is used to select forward reads from all **ab1** files.

  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all filename of reverse reads. :code:`grepl` function in R is used to select reverse reads from all **ab1** files.


It is important to carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`.
The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. :code:`_[0-9]*_F.ab1$`, :code:`_[0-9]*_R.ab1$` for matching forward and reverse reads are highly suggested and are used as default. Users are strongly recommended to follow the default file naming convention to reduce any change of error.


For basic preparing input files example, please go to :ref:`Beginner Guide`. Here, we have another more complicated example.


.. _SangerAlignment_file_structure:
.. figure::  ../image/SangerAlignment_file_structure.png
   :align:   center
   :scale:   50 %

   Figure 1. Input ab1 files inside the parent directory, :code:`./tmp/`.


:ref:`Figure_1<SangerAlignment_file_structure>` shows the directory that users need to prepare and the file naming convention. Here, we explain the three parameters in this example.



FASTA files
+++++++++++

|

*SangerAlignment* creation
--------------------------

Inputs Definition
+++++++++++++++++

Slots Definition
++++++++++++++++

|

Update quality trimming parameters
----------------------------------

|

*SangerAlignment* Shiny app
---------------------------

|


Writing FASTA files
-------------------

|

Generating *SangerAlignment* report
-----------------------------------
