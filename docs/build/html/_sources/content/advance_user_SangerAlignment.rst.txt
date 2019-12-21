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
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all filename in forward direction. :code:`grepl` function in R is used to select forward reads from all **ab1** files.
  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all filename in reverse direction. :code:`grepl` function in R is used to select reverse reads from all **ab1** files.





For basic input files preparation example, please go to :ref:`Beginner Guide`. Here, we have another more complicated example.


.. _SangerAlignment_file_structure_complex:
.. figure::  ../image/SangerAlignment_file_structure_complex.png
   :align:   center
   :scale:   120 %

   Figure 1. Input ab1 files inside the parent directory, :code:`./tmp/`.


:ref:`Figure_1<SangerAlignment_file_structure_complex>` shows the file naming convention and directory hierarchy. In this example, the parent directory is :code:`extdata` and the directories in first layer are :code:`Allolobophora_chlorotica` and :code:`Drosophila_melanogaster`. All target **ab1** files need to be inside parent directory but it is not necessary to put them in the same level of directory. sangeranalyseR will recursively search all files with **.ab1** file extension and automatically group reads with the same contig name. The direction of reads in each contig will be grouped by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp` with filenames. Therefore, it is important to carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1$"`, :code:`"_[0-9]*_R.ab1$"` for matching forward and reverse reads are highly suggested and are used as default. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _SangerAlignment_filename_convention:
.. figure::  ../image/SangerAlignment_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 2. Suggested **ab1** file naming convention.

:ref:`Figure_2<SangerAlignment_filename_convention>` shows the suggested **ab1** file naming convention. Users are strongly recommended to follow this file naming convention and use the default :code:`suffixForwardRegExp` : :code:`"_[0-9]*_F.ab1$"` and :code:`suffixReverseRegExp` : :code:`"_[0-9]*_R.ab1$"` to reduce any chance of error.



FASTA files
+++++++++++

|

*SangerAlignment* creation
--------------------------
After preparing the input directory, we can create the *SangerAlignment* S4 instance by running :code:`SangerAlignment` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it make instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list all possible parameters that users can assign. For a simpler command, please go to :ref:`Quick Commands`.

.. code-block:: R

   sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                      parentDirectory      = "./tmp/",
                                      suffixForwardRegExp  = "_[0-9]*_F.ab1$",
                                      suffixReverseRegExp  = "_[0-9]*_R.ab1$",
                                      refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                      TrimmingMethod        = "M1",
                                      M1TrimmingCutoff      = 0.0001,
                                      M2CutoffQualityScore  = NULL,
                                      M2SlidingWindowSize   = NULL,
                                      baseNumPerRow         = 100,
                                      heightPerRow          = 200,
                                      signalRatioCutoff     = 0.33,
                                      showTrimmed           = TRUE,
                                      minReadsNum           = 2,
                                      minReadLength         = 20,
                                      minFractionCall       = 0.5,
                                      maxFractionLost       = 0.5,
                                      geneticCode           = GENETIC_CODE,
                                      acceptStopCodons      = TRUE,
                                      readingFrame          = 1,
                                      minFractionCallSA     = 0.5,
                                      maxFractionLostSA     = 0.5,
                                      processorsNum         = NULL)

Inputs & Slots Definition
+++++++++++++++++++++++++
The inputs of :code:`SangerAlignment` constructor function and :code:`new` method are same. For more details, please refer to `sangeranalyseR reference manul <http://packages.python.org/an_example_pypi_project/>`_.

|

Update quality trimming parameters
----------------------------------
:code:`SangerAlignment` constructor function will apply the given quality trimming parameters to all reads. After creating the *SangerAlignment* S4 instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will change the trimming parameters for each read and redo alignment in both *SangerContig* and *SangerAlignment* levels. If users want to do quality trimming read by read, please
read :ref:`*SangerAlignment* Shiny app`.

.. code-block:: R

   newSangerAlignment <- updateQualityParam(sangerAlignment,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 29,
                                            M2SlidingWindowSize  = 15)

|

*SangerAlignment* Shiny app
---------------------------
We create an interactive local Shiny app for users to go inside each *SangerRead* and *SangerContig* in *SangerAlignment*.

Users only need to run one function with previously created *SangerAlignment* instance as input and the Shiny app will pop up.

.. code-block:: R

   launchAppSA(sangerAlignment)


*SangerAlignment* Page
++++++++++++++++++++++

.. _SangerAlignment_ShinyApp_1:
.. figure::  ../image/SangerAlignment_ShinyApp_1.png
   :align:   center
   :scale:   40 %

   Figure 2. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_2:
.. figure::  ../image/SangerAlignment_ShinyApp_2.png
   :align:   center
   :scale:   40 %

   Figure 3. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_3:
.. figure::  ../image/SangerAlignment_ShinyApp_3.png
   :align:   center
   :scale:   40 %

   Figure 4. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_4:
.. figure::  ../image/SangerAlignment_ShinyApp_4.png
   :align:   center
   :scale:   40 %

   Figure 5. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_5:
.. figure::  ../image/SangerAlignment_ShinyApp_5.png
   :align:   center
   :scale:   40 %

   Figure 6. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_6:
.. figure::  ../image/SangerAlignment_ShinyApp_6.png
   :align:   center
   :scale:   40 %

   Figure 7. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_7:
.. figure::  ../image/SangerAlignment_ShinyApp_7.png
   :align:   center
   :scale:   40 %

   Figure 8. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_8:
.. figure::  ../image/SangerAlignment_ShinyApp_8.png
   :align:   center
   :scale:   40 %

   Figure 9. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_9:
.. figure::  ../image/SangerAlignment_ShinyApp_9.png
   :align:   center
   :scale:   40 %

   Figure 10. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_10:
.. figure::  ../image/SangerAlignment_ShinyApp_10.png
   :align:   center
   :scale:   40 %

   Figure 11. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_11:
.. figure::  ../image/SangerAlignment_ShinyApp_11.png
   :align:   center
   :scale:   40 %

   Figure 12. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_12:
.. figure::  ../image/SangerAlignment_ShinyApp_12.png
   :align:   center
   :scale:   40 %

   Figure 13. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_13:
.. figure::  ../image/SangerAlignment_ShinyApp_13.png
   :align:   center
   :scale:   40 %

   Figure 14. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_14:
.. figure::  ../image/SangerAlignment_ShinyApp_14.png
   :align:   center
   :scale:   40 %

   Figure 15. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_15:
.. figure::  ../image/SangerAlignment_ShinyApp_15.png
   :align:   center
   :scale:   40 %

   Figure 16. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_16:
.. figure::  ../image/SangerAlignment_ShinyApp_16.png
   :align:   center
   :scale:   40 %

   Figure 17. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_17:
.. figure::  ../image/SangerAlignment_ShinyApp_17.png
   :align:   center
   :scale:   40 %

   Figure 18. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

.. _SangerAlignment_ShinyApp_18:
.. figure::  ../image/SangerAlignment_ShinyApp_18.png
   :align:   center
   :scale:   40 %

   Figure 19. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.


*SangerContig* Page
+++++++++++++++++++

*SangerRead* Page
+++++++++++++++++

|


Inside Shiny app, users can change the quality trimming and chromatogram parameters for each read and redo alignment.

Writing FASTA files
-------------------

|

Generating *SangerAlignment* report
-----------------------------------
