Beginner Guide
==============

Please first install :code:`sangeranalyseR` and load it in your R console. If you haven't finished the above steps, please refer to :ref:`Getting Start`.

Input files preparation
-----------------------
sangeranalyseR takes ABIF files as input and constructs contigs for each species. All the contigs will be aligned and a



The results contain contigs alignment, evolution tree,

All reads in each contig will be aligned first to get a contig,

All the files must have ab1 in its suffix. The reads that belong to the same group should have the same contig name. Forward and reverse read also need to specify in the filename.

There are three parameters that users need to provide.
First is parent directory :code:`parentDirectory`
Second is forward suffix :code:`suffixForwardRegExp`.
Third is reverse suffix :code:`suffixReverseRegExp`.

All files should be placed under the :code:`parentDirectory`. They do not need to be in the same layer of the directory. sangeranalyseR will recursively search all the directories inside :code:`parentDirectory` and find all files that end with ab1 and match the forward or reverse suffix. Files with the same contig names will be categorized in the same group and aligned into a contig. Therefore, it is very important to make sure that filenames are correctly and systematically named.

.. _Figure_1:
.. figure::  ../image/figure_1.png
   :align:   center


Figure_1_ is the example of the hierarchy of input ab1 files.

:code:`parentDirectory`: "./tmp"
:code:`suffixForwardRegExp`: "_F_[0-9]*.ab1".
:code:`suffixReverseRegExp`: "_R_[0-9]*.ab1".

|

*SangerAlignment* creation
--------------------------
.. code-block:: R

   rawDir <- system.file("extdata", package="sangeranalyseR")
   parentDir <- file.path(rawDir, "Allolobophora_chlorotica", "RBNII")

.. code-block:: R

   sangerAlignment <- SangerAlignment(parentDirectory = parentDir,
                                      suffixForwardRegExp = "_[F]_[0-9]*.ab1",
                                      suffixReverseRegExp = "_[R]_[0-9]*.ab1")

|

Launch Shiny App
----------------
.. code-block:: R

   launchAppSA(sangerAlignment)

.. _Figure_2:
.. figure::  ../image/figure_2.png
   :align:   center

Figure_2_ is the

|

Writing FASTA file
------------------
.. code-block:: R

   writeFastaSA(sangerAlignment)

|

Generating report
-----------------
.. code-block:: R

   generateReportSA(sangerAlignment)

.. _Figure_3:
.. figure::  ../image/figure_3.png
   :align:   center

Figure_3_ is the

.. _Figure_4:
.. figure::  ../image/figure_4.png
   :align:   center

Figure_4_ is the

|

Where to go from here ?
-----------------------
Congratulation, you have finished the :ref:`Beginner Guide`. For learning more

:ref:`Beginner Guide`
