Advanced User Guide - *SangerAlignment* (**AB1**)
=================================================


*SangerAlignment* is the highest class level in sangeranalyseR showed in :ref:`Figure_1<SangerAlignment_hierachy>`. It contains *SangerContig* list and the contigs alignment result. Users can access to *SangerContig* and *SangerRead* instance inside *SangerAlignment* instance. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerAlignment* level.

.. _SangerAlignment_hierachy:
.. figure::  ../image/SangerAlignment_hierachy.png
   :align:   center
   :scale:   20 %

   Figure 1. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

|

Preparing *SangerAlignment* **AB1** input
-----------------------------------------
Users can choose to input **AB1** or **FASTA** as their input file format.

The main input file format to create *SangerAlignment* instance is **AB1**. Before starting the analysis, users need to prepare a directory which contains all the **AB1** files. Here are some filename regulations:

.. note::

    *  All the input files must have **.ab1** as its file extension
    *  The reads that belong to the same contig must have the same contig name in its filename.
    *  Forward or reverse direction also has to be specified in the filename.


There are three parameters, :code:`parentDirectory`, :code:`suffixForwardRegExp`, and :code:`suffixReverseRegExp`, that users need to provide so that program can automatically group all **AB1** files.

.. note::

  * :code:`parentDirectory`: The root directory that contains all the **AB1** files. It can be absolute or relative path. We suggest users to put only target **AB1** files inside this directory without other unrelated files.
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all filename in forward direction. :code:`grepl` function in R is used to select forward reads from all **AB1** files.
  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all filename in reverse direction. :code:`grepl` function in R is used to select reverse reads from all **AB1** files.





For basic input files preparation example, please go to :ref:`Beginner Guide`. Here, we have another more complicated example.


.. _SangerAlignment_file_structure_complex:
.. figure::  ../image/SangerAlignment_file_structure_complex.png
   :align:   center
   :scale:   120 %

   Figure 2. Input ab1 files inside the parent directory, :code:`./tmp/`.


:ref:`Figure_2<SangerAlignment_file_structure_complex>` shows the file naming regulation and directory hierarchy. In this example, the parent directory is :code:`extdata` and the directories in first layer are :code:`Allolobophora_chlorotica` and :code:`Drosophila_melanogaster`. All target **AB1** files need to be inside parent directory but it is not necessary to put them in the same level of directory. sangeranalyseR will recursively search all files with **.ab1** file extension and automatically group reads with the same contig name. The direction of reads in each contig will be grouped by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp` with filenames. Therefore, it is important to carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1$"`, :code:`"_[0-9]*_R.ab1$"` for matching forward and reverse reads are highly suggested and are used as default. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerAlignment:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file naming regulation - *SangerAlignment*.

:ref:`Figure_3<sangeranalyseR_filename_convention_SangerAlignment>` shows the suggested **AB1** file naming regulation. Users are strongly recommended to follow this file naming regulation and use the default :code:`suffixForwardRegExp` : ":code:`_[0-9]*_F.ab1$`" and :code:`suffixReverseRegExp` : ":code:`_[0-9]*_R.ab1$`" to reduce any chance of error.

|

Creating *SangerAlignment* instance from **AB1**
------------------------------------------------
After preparing the input directory, we can create the *SangerAlignment* S4 instance by running :code:`SangerAlignment` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list important parameters. For a simpler command, please go to :ref:`Quick Commands`.

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

The inputs of :code:`SangerAlignment` constructor function and :code:`new` method are same. For more details about *SangerAlignment* inputs & slots definition, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_.

|

Updating *SangerAlignment* quality trimming parameters
------------------------------------------------------
In the previous :ref:`Creating *SangerAlignment* instance from **AB1**` part, the constructor function will apply the quality trimming parameters to all reads. After creating the *SangerAlignment* S4 instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will update all reads with the new trimming parameters and redo reads alignment in *SangerContig* and contigs alignment in *SangerAlignment*. If users want to do quality trimming read by read instead all at once, please read :ref:`Launching *SangerAlignment* Shiny app`.

.. code-block:: R

   newSangerAlignment <- updateQualityParam(sangerAlignment,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 29,
                                            M2SlidingWindowSize  = 15)

|

Launching *SangerAlignment* Shiny app
-------------------------------------
We create an interactive local Shiny app for users to go into each *SangerRead* and *SangerContig* in *SangerAlignment* instance. Users only need to run one function with previously created instance as input and the *SangerAlignment* Shiny app will pop up. Here, we will go through pages in the three levels.

.. code-block:: R

   launchAppSA(sangerAlignment)


*SangerAlignment* page (SA app)
+++++++++++++++++++++++++++++++
:ref:`Figure 4<SangerAlignment_ShinyApp_1>` is the initial page and the toppest layer of *SangerAlignment* App. It provides basic parameters in *SangerAlignment* instance, contigs alignment result and phylogenetic tree etc. After updating all the trimmed reads, users need to click “Re-calculate Contigs Alignment” button to do contigs alignment again. From the left-hand side panel, we can clearly see the hierarchy of the *SangerAlignment* S4 instance and easily access to all reads and contigs in it.

.. _SangerAlignment_ShinyApp_1:
.. figure::  ../image/SangerAlignment_ShinyApp_1.png
   :align:   center
   :scale:   25 %

   Figure 4. *SangerAlignment* Shiny app initial page - *SangerAlignment* Page.

Scroll down a bit, users can see the contigs alignment result generated by `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ R package embedded in *SangerAlignment* page. :ref:`Figure 5<SangerAlignment_ShinyApp_2>` shows the contigs alignment result.

.. _SangerAlignment_ShinyApp_2:
.. figure::  ../image/SangerAlignment_ShinyApp_2.png
   :align:   center
   :scale:   30 %

   Figure 5. *SangerAlignment* Page - contigs alignment result.

In *SangerAlignment* page, the phylogenetic tree result is provided as well (:ref:`Figure 6<SangerAlignment_ShinyApp_3>`). The tree is generated by `ape <https://cran.r-project.org/web/packages/ape/index.html>`_ R package which uses NJ algorithm.

.. _SangerAlignment_ShinyApp_3:
.. figure::  ../image/SangerAlignment_ShinyApp_3.png
   :align:   center
   :scale:   60 %

   Figure 6. *SangerAlignment* Page - phylogenetic tree result.


*SangerContig* page (SA app)
++++++++++++++++++++++++++++
Now, let's go to the page in the next level, *SangerContig* page. Users can click into all contigs and check their results. :ref:`Figure 7<SangerAlignment_ShinyApp_5>` shows the overview page of Contig 1. Notice that there is a red “Re-calculate Contig” button. After changing the quality trimming parameters, users need to click the button before checking the results below in order to get the updated information.

.. _SangerAlignment_ShinyApp_5:
.. figure::  ../image/SangerAlignment_ShinyApp_5.png
   :align:   center
   :scale:   25 %

   Figure 7. *SangerAlignment* Shiny app - *SangerContig* page.

The information provided in this page includes : “input parameters”, “genetic code table”, “reference amino acid sequence”, “reads alignment”, “difference data frame”, “dendrogram”, “sample distance heatmap”, “indels data frame”, “stop codons data frame”. :ref:`Figure 8<SangerAlignment_ShinyApp_6>` and :ref:`Figure 9<SangerAlignment_ShinyApp_7>` show part of the results in the *SangerContig* page. The results are dynamic based on the trimming parameters from user inputs.

.. _SangerAlignment_ShinyApp_6:
.. figure::  ../image/SangerAlignment_ShinyApp_6.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - contig-related parameters, genetic code and reference amino acid sequence.

.. _SangerAlignment_ShinyApp_7:
.. figure::  ../image/SangerAlignment_ShinyApp_7.png
   :align:   center
   :scale:   30 %

   Figure 9. *SangerContig* page - reads alignment and difference data frame.


*SangerRead* page (SA app)
++++++++++++++++++++++++++
Now, let's go to the page in the lowest level, *SangerRead* page. *SangerRead* page contains all details of a read including its trimming and chromatogram inputs and results. All reads are in "forward" or "reverse" direction. Under *SangerContig* page, there are two expendable tabs, “Forward Reads” and “Reverse Reads” storing the corresponding reads on the left-hand side navigation panel in :ref:`Figure 10<SangerAlignment_ShinyApp_8>`. In this example, there are one read in each tab and :ref:`Figure 10<SangerAlignment_ShinyApp_8>` shows the “1 - 1 Forward Read” page. It provides basic information, quality trimming inputs, chromatogram plotting inputs etc. Primary/secondary sequences in this figure are dynamic based on the :code:`signalRatioCutoff` value for base calling and the length of them are always same. Another thing to mention is that primary/secondary sequences and the sequences in the chromatogram in :ref:`Figure 15<SangerAlignment_ShinyApp_14>` below will always be same after trimming and their color codings for A/T/C/G are same as well.

.. _SangerAlignment_ShinyApp_8:
.. figure::  ../image/SangerAlignment_ShinyApp_8.png
   :align:   center
   :scale:   25 %

   Figure 10. *SangerAlignment* Shiny app - *SangerRead* page.

In quality trimming steps, we removes fragment at both ends of sequencing reads with low quality score. It is important because trimmed reads would improves alignment results. :ref:`Figure 11<SangerAlignment_ShinyApp_9>` shows the UI for Trimming Method 1 (M1): ‘Modified Mott Trimming’. This method is implemented in `Phred <http://www.phrap.org/phredphrapconsed.html>`_. Users can change the cutoff score and click “Apply Trimming Parameters" button to update the UI. The value of input must be between 0 and 1. If the input is invalid, the cutoff score will be set to default 0.0001.

.. _SangerAlignment_ShinyApp_9:
.. figure::  ../image/SangerAlignment_ShinyApp_9.png
   :align:   center
   :scale:   45 %

   Figure 11. *SangerRead* page - Trimming Method 1 (M1): ‘Modified Mott Trimming’ UI.

:ref:`Figure 12<SangerAlignment_ShinyApp_10>` shows another quality trimming methods for users to choose from, Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’. This method is implemented in `Trimmomatics <http://www.usadellab.org/cms/?page=trimmomatic>`_. Users can change the cutoff quality score as well as sliding window size and click “Apply Trimming Parameters" button to update the UI. The value of cutoff quality score must be between 0 and 60 (default 20); the value of sliding window size must be between 0 and 40 (default 10). If the inputs are invalid, their values will be set to default.

.. _SangerAlignment_ShinyApp_10:
.. figure::  ../image/SangerAlignment_ShinyApp_10.png
   :align:   center
   :scale:   45 %

   Figure 12. *SangerRead* page - Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’ UI.

:ref:`Figure 13<SangerAlignment_ShinyApp_11>` shows the quality report before and after trimming. After clicking the “Apply Trimming Parameters” button, the values of these information boxes will be updated to the latest values.

.. _SangerAlignment_ShinyApp_11:
.. figure::  ../image/SangerAlignment_ShinyApp_11.png
   :align:   center
   :scale:   45 %

   Figure 13. *SangerRead* page - read quality report before / after trimming.

In :ref:`Figure 14<SangerAlignment_ShinyApp_13>`, the x-axis is the index of the base pairs; the y-axis is the Phred quality score. The green horizontal bar at the top is the raw read region and the orange horizontal bar represents the trimmed read region. Both :ref:`Figure 14<SangerAlignment_ShinyApp_13>` timming plot and :ref:`Figure 15<SangerAlignment_ShinyApp_14>` chromatogram will be updated once users change the quality trimming parameters and click the “Apply Trimming Parameters" button in :ref:`Figure 15<SangerAlignment_ShinyApp_14>`.

.. _SangerAlignment_ShinyApp_13:
.. figure::  ../image/SangerAlignment_ShinyApp_13.png
   :align:   center
   :scale:   50 %

   Figure 14. *SangerRead* page - quality trimming plot.

If we only see primary and secondary sequences in the table, we will loose some variations. Chromatogram is very helpful to check the peak resolution. :ref:`Figure 15<SangerAlignment_ShinyApp_14>` shows the panel of plotting chromatogram. Users can change four parameters: :code:`Base Number Per Row`, :code:`Height Per Row`, :code:`Signal Ratio Cutoff`, and :code:`Show Trimmed Region`. Among them, :code:`Signal Ratio Cutoff` is a key parameter. If its value is default value 0.33, it indicates that the lower peak should be at least 1/3rd as high as the higher peak for it count as a secondary peak.

.. _SangerAlignment_ShinyApp_14:
.. figure::  ../image/SangerAlignment_ShinyApp_14.png
   :align:   center
   :scale:   45 %

   Figure 15. *SangerRead* page - chromatogram panel.

Here is an example of applying new chromatogram parameters. We click “Show Trimmed Region” to set its value from FALSE to TRUE. :ref:`Figure 16<SangerAlignment_ShinyApp_15>` shows the loading notification popup during base calling and chromatogram plotting.

.. _SangerAlignment_ShinyApp_15:
.. figure::  ../image/SangerAlignment_ShinyApp_15.png
   :align:   center
   :scale:   45 %

   Figure 16. *SangerRead* page - loading notification popup during replotting chromatogram.

After replotting the chromatogram, trimmed region is showed in red striped region. :ref:`Figure 17<SangerAlignment_ShinyApp_16>` shows part of the the chromatogram (1 bp ~ 240 bp). Moreover, chromatogram will be replotted when trimmed positions or chromatogram parameters are updated.

.. _SangerAlignment_ShinyApp_16:
.. figure::  ../image/SangerAlignment_ShinyApp_16.png
   :align:   center
   :scale:   50 %

   Figure 17. *SangerRead* page - chromatogram with trimmed region showed.

To let users browse the trimmed primary/secondary sequences without finding “Trimming Start Point” and “Trimming End Point” by themselves, we provide the final trimmed primary/secondary sequences that will be used for reads alignment in table format with quality scores in :ref:`Figure 18<SangerAlignment_ShinyApp_17>`. Frameshift amino acid sequences are also provided.

.. _SangerAlignment_ShinyApp_17:
.. figure::  ../image/SangerAlignment_ShinyApp_17.png
   :align:   center
   :scale:   45 %

   Figure 18. *SangerRead* page - trimmed primary/secondary sequences and Phred quality score in table format.

We have updated the trimming and chromatogram parameters for each read. Now, we need to click “Re-calculate contig” button to do alignment again. Last but not least, we can save all data into a new ‘SangerContig’ S4 instance by clicking “Save S4 Object button”. New S4 instance would be saved in Rda format. Users can run :code:`readRDS` function to load it into current R environment. :ref:`Figure 19<SangerAlignment_ShinyApp_18>` shows some hints in the save notification popup.

.. _SangerAlignment_ShinyApp_18:
.. figure::  ../image/SangerAlignment_ShinyApp_18.png
   :align:   center
   :scale:   40 %

   Figure 19. *SangerRead* page - saving notification popup.


|


Writing *SangerAlignment* FASTA files :sub:`(AB1)`
--------------------------------------------------
Users can write the *SangerAlignment* instance to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`contigs_unalignment`: Writing contigs into a single **FASTA** file.
* :code:`contigs_alignment`: Writing contigs alignment and contigs consensus read to a single **FASTA** file.
* :code:`all_reads`: Writing all reads to a single **FASTA** file.
* :code:`all`: Writing contigs, contigs alignment, and all reads into three different files.

Below is the one-line function that users need to run. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFastaSA` function.

.. code-block:: R

   writeFastaSA(sangerAlignment,
                outputDir         = tempdir(),
                compress          = FALSE,
                compression_level = NA,
                selection         = "all")

|

Generating *SangerAlignment* report :sub:`(AB1)`
------------------------------------------------
Last but not least, users can save *SangerAlignment* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.

There are two parameters, :code:`includeSangerContig` and :code:`includeSangerRead`, for users to decide to which level the *SangerAlignment* report will go. Moreover, after the reports are generated, users can easily navigate through reports in different levels within the **HTML** file.

* :code:`includeSangerContig`: Whether users want to generate the report of each *SangerContig* in *SangerAlignment*.
* :code:`includeSangerRead`: If :code:`includeSangerContig` is :code:`TRUE`, then users can set this value to decide whether they want to include *SangerRead* reports in each *SangerContig*.

One thing to pay attention to is that if users have many reads, it would take quite a long time to write out all reports. If users only want to generate the contigs alignment, remember to set :code:`includeSangerContig` and :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReportSA(sangerAlignment,
                    outputDir           = tempdir(),
                    includeSangerContig = TRUE,
                    includeSangerRead   = TRUE)
