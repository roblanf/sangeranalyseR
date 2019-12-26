Advanced User Guide - *SangerContig*
====================================

*SangerContig* is the second level in sangeranalyseR showed in :ref:`Figure_1<SangerContig_hierarchy>` which corresponds to a contig in Sanger sequencing. Among slots inside it, there are two lists, forward and reverse read list, storing *SangerRead* in the corresponding direction. In this section, we are going to go through details about sangeranalyseR data analysis in *SangerContig* level.

.. _SangerContig_hierarchy:
.. figure::  ../image/SangerContig_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerContig* level.

|

*SangerContig* input files preparation
--------------------------------------

Users can choose to input **ab1** or **FASTA** as their input file format.

ab1 files
+++++++++

The main input file format to create *SangerRead* instance is **ab1**. Before starting the analysis, users need to prepare all **ab1** files inside one directory. This directory is the parent directory and all **ab1** files must be in the first layer of it; in other words, there should not be any directory containing any **ab1** files inside the parent directory. Because sangeranalyseR will group **ab1** files based on their direction automatically, users have to follow the filename regulations below:

.. note::

    *  All the input files must have **.ab1** as its file extension
    *  All the input files must have the same contig name in its filename.
    *  Forward or reverse direction also has to be specified in the filename.


There are four parameters, :code:`parentDirectory`, :code:`contigName`, :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`,that users need to provide so that program can automatically match correct **ab1** files and divide them into forward and reverse direction.

.. note::

  * :code:`parentDirectory`: The root directory that contains all the **ab1** files. It can be absolute or relative path. We suggest users to put only target **ab1** files inside this directory without other unrelated files.
  * :code:`contigName`: The value of this parameter is a regular expression that matches filenames that are going to be included in the *SangerContig* level analysis. :code:`grepl` function in R is used.
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all filenames in forward direction. :code:`grepl` function in R is used to select forward reads from all **ab1** files.
  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all filenames in reverse direction. :code:`grepl` function in R is used to select reverse reads from all **ab1** files.

Here, we have an example:




.. _SangerContig_file_structure:
.. figure::  ../image/SangerContig_file_structure.png
   :align:   center
   :scale:   90 %

   Figure 2. *SangerContig* filename regulation.

:ref:`Figure_2<SangerContig_file_structure>` shows the file naming convention and hierarchy. In this example, :code:`ACHLO` is the parent directory that contains all **ab1** files. They must be in the first layer of the directory.

sangeranalyseR will first match the :code:`contigName` to exclude unrelated files and then separate the forward and reverse reads by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. Therefore, it is important to make sure all target **ab1** files share the same :code:`contigName` and carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1$"`, :code:`"_[0-9]*_R.ab1$"` for matching forward and reverse reads are highly suggested and are used as default. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerContig:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **ab1** file naming convention - *SangerContig*.

:ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig>` shows the suggested **ab1** file naming convention. Users are strongly recommended to follow this file naming convention and use the default :code:`suffixForwardRegExp` : ":code:`_[0-9]*_F.ab1$`" and :code:`suffixReverseRegExp` : ":code:`_[0-9]*_R.ab1$`" to reduce any chance of error.


FASTA files
+++++++++++

|


*SangerContig* instance creation
--------------------------------

After preparing the input directory, we can create the *SangerContig* S4 instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list all possible parameters.

.. code-block:: R

    sangerContig <- SangerContig(inputSource            = "ABIF",
                                 parentDirectory        = "./tmp/",
                                 contigName             = "ACHLO006-09[LCO1490_t1.HCO2198_t1]",
                                 suffixForwardRegExp    = "[0-9]*_F.ab1",
                                 suffixReverseRegExp    = "[0-9]*_R.ab1",
                                 TrimmingMethod         = "M1",
                                 M1TrimmingCutoff       = 0.0001,
                                 M2CutoffQualityScore   = NULL,
                                 M2SlidingWindowSize    = NULL,
                                 baseNumPerRow          = 100,
                                 heightPerRow           = 200,
                                 signalRatioCutoff      = 0.33,
                                 showTrimmed            = TRUE,
                                 refAminoAcidSeq        = "",
                                 minReadsNum            = 2,
                                 minReadLength          = 20,
                                 minFractionCall        = 0.5,
                                 maxFractionLost        = 0.5,
                                 geneticCode            = GENETIC_CODE,
                                 acceptStopCodons       = TRUE,
                                 readingFrame           = 1,
                                 processorsNum          = NULL)


In this example, :code:`contigName` is set to :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]"`, so only :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_1_F.ab1"` and :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_2_R.ab1"` will be selected to align to a contig.


*SangerContig* inputs & slots definition
++++++++++++++++++++++++++++++++++++++++

The inputs of :code:`SangerContig` constructor function and :code:`new` method are same. For more details, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_.

|


*SangerContig* quality trimming parameters updating
---------------------------------------------------

In the previous :ref:`*SangerContig* instance creation` part, the constructor function will apply the quality trimming parameters to all reads. After creating the SangerContig S4 instance, users can change the trimming parameters by running updateQualityParam function which will update all reads with the new trimming parameters and redo reads alignment. If users want to do quality trimming read by read instead all at once, please read :ref:`*SangerContig* Shiny app`.

.. code-block:: R

   newSangerContig <- updateQualityParam(sangerContig,
                                         TrimmingMethod       = "M2",
                                         M1TrimmingCutoff     = NULL,
                                         M2CutoffQualityScore = 29,
                                         M2SlidingWindowSize  = 15)

|

*SangerContig* Shiny app
------------------------

We create an interactive local Shiny app for users to go into each *SangerRead* in *SangerContig* instance. Users only need to run one function with previously created instance as input and the *SangerContig* Shiny app will pop up. Here, we will go through pages in the two levels, *SangerRead* and *SangerContig* pages.

*SangerContig* page
+++++++++++++++++++
*SangerContig* page is the initial page of *SangerContig* Shiny app.

.. _SangerContig_shiny_SangerContig_page:
.. figure::  ../image/SangerContig_shiny_SangerContig_page.png
   :align:   center
   :scale:   30 %

   Figure 4. *SangerContig* Shiny app initial page - *SangerContig* page.

.. _SangerContig_shiny_alignment_differenceDF:
.. figure::  ../image/SangerContig_shiny_alignment_differenceDF.png
   :align:   center
   :scale:   30 %

   Figure 5. *SangerContig* page - *SangerRead* alignment.


.. _SangerContig_shiny_dendrogram:
.. figure::  ../image/SangerContig_shiny_dendrogram.png
   :align:   center
   :scale:   30 %

   Figure 6. *SangerContig* page - dendrogram.


.. _SangerContig_shiny_samples_distance:
.. figure::  ../image/SangerContig_shiny_samples_distance.png
   :align:   center
   :scale:   30 %

   Figure 7. *SangerContig* page - samples distance.

.. _SangerContig_shiny_indelsDF_stopcodonsDF:
.. figure::  ../image/SangerContig_shiny_indelsDF_stopcodonsDF.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - indels and stop codons data frame.


*SangerRead* page
+++++++++++++++++
.. _SangerContig_shiny_SangerRead_page:
.. figure::  ../image/SangerContig_shiny_SangerRead_page.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - indels and stop codons data frame.

.. _SangerContig_shiny_SangerRead_page:
.. figure::  ../image/SangerContig_shiny_SangerRead_page.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - indels and stop codons data frame.

|
