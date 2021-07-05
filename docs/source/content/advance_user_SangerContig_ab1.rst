Advanced User Guide - *SangerContig* (**AB1**)
==============================================

*SangerContig* is the intermediate level in sangeranalyseR (:ref:`Figure_1<SangerContig_hierarchy>`), and each *SangerContig* instance corresponds to a contig in a Sanger sequencing experiment. Among its slots, there are two lists, forward and reverse read list, storing *SangerRead* in the corresponding direction. In this section, we are going to go through details about sangeranalyseR data analysis in *SangerContig* level with **AB1** file input.

.. _SangerContig_hierarchy:
.. figure::  ../image/SangerContig_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerContig* level.

|

Preparing *SangerContig* **AB1** inputs
++++++++++++++++++++++++++++++++++++++++

The main input file format to create *SangerRead* instance is **AB1**. Before starting the analysis, users need to prepare all **AB1** files inside one directory, and all of them must be in the first layer of that directory. In other words, there should  be no directories containing any **AB1** files inside it. There are two ways for users to group their **AB1** files which are **"regular expression matching"** and **"CSV file matching"**, and following are instructions of how to prepare and name your **AB1** input files.



"regular expression matching" inputs preparation
-------------------------------------------------
For regular expression matching method, sangeranalyseR will group **AB1** files based on their contig names and read directions in their filenames automatically; therefore, users have to follow the file naming regulations below:

.. note::

    *  All the input files must have **.ab1** as its file extension
    *  All the input files must have the same contig name in its filename.
    *  Forward or reverse direction also has to be specified in the filename.




There are four parameters, :code:`ABIF_Directory`, :code:`contigName`, :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`,that users need to provide so that program can automatically match correct **AB1** files and divide them into forward and reverse direction.

.. note::

  * :code:`ABIF_Directory`: The root directory that contains all the **AB1** files. It can be absolute or relative path. We suggest users to put only target **AB1** files inside this directory and don't include any other unrelated files.
  * :code:`contigName`: It is a regular expression that matches filenames that are going to be included in the *SangerContig* level analysis. :code:`grepl` function in R is used.
  * :code:`REGEX_SuffixForward`: It is a regular expression that matches all filenames in forward direction. :code:`grepl` function in R is used to select forward reads from all **AB1** files.
  * :code:`REGEX_SuffixReverse`: It is a regular expression that matches all filenames in reverse direction. :code:`grepl` function in R is used to select reverse reads from all **AB1** files.

Here is an example:


.. _SangerContig_file_structure:
.. figure::  ../image/SangerContig_file_structure.png
   :align:   center
   :scale:   60 %

   Figure 2. *SangerContig* filename regulation.

:ref:`Figure_2<SangerContig_file_structure>` shows the file naming regulation and hierarchy. In this example, :code:`ACHLO` is the parent directory that contains all **AB1** files. They must be in the first layer of the directory.

sangeranalyseR will first match the :code:`contigName` to exclude unrelated files and then separate the forward and reverse reads by matching :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. Therefore, it is important to make sure that all target **AB1** files share the same :code:`contigName` and carefully select :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1"`, :code:`"_[0-9]*_R.ab1"` for matching forward and reverse reads are highly suggested and are used as default. Moreover, it is a good habit to index your reads in the same contig group since there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerContig:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file naming regulation - *SangerContig*.

:ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig>` shows the suggested **AB1** file naming regulation. Users are strongly recommended to follow this file naming regulation and use the default :code:`REGEX_SuffixForward` : ":code:`_[0-9]*_F.ab1`" and :code:`REGEX_SuffixReverse` : ":code:`_[0-9]*_R.ab1`" to reduce any chance of error.



"CSV file matching" inputs preparation
---------------------------------------
For CSV file matching method, sangeranalyseR will group **AB1** files based on the information in the CSV file automatically; therefore, users have to prepare a CSV with the following regulation:

.. note::

    *  It needs to contain three columns: "reads", "direction", and "contig".
    *  "reads" column stores the filename of the **AB1** files.
    *  "direction" column stores the direction of the reads. It must be "F" (forward) or "R" (reverse).
    *  "contig" column stores the contig name that each read blongs.


There are three parameters, :code:`ABIF_Directory`, :code:`contigName`, and :code:`CSV_NamesConversion`,that users need to provide so that program can automatically match correct **AB1** files and divide them into forward and reverse direction.

.. note::

  * :code:`ABIF_Directory`: The root directory that contains all the **AB1** files. It can be absolute or relative path. We suggest users to put only target **AB1** files inside this directory and don't include any other unrelated files.
  * :code:`contigName`: It is a regular expression that matches filenames that are going to be included in the *SangerContig* level analysis. :code:`grepl` function in R is used.
  * :code:`CSV_NamesConversion`: It is a regular expression that matches all filenames in forward direction. :code:`grepl` function in R is used to select forward reads from all **AB1** files.


|


Creating *SangerContig* instance from **AB1**
++++++++++++++++++++++++++++++++++++++++++++++

After preparing the input directory, we can create the *SangerContig* S4 instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list important parameters.

.. code-block:: R

   # using `constructor` function to create SangerRead instance
    sangerContig <- SangerContig(inputSource            = "ABIF",
                                 ABIF_Directory         = "./tmp/",
                                 contigName             = "Achl_ACHLO006-09",
                                 REGEX_SuffixForward    = "_[0-9]*_F.ab1",
                                 REGEX_SuffixReverse    = "_[0-9]*_R.ab1",
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

   # using `new` method to create SangerRead instance
   sangerContig <- new("SangerContig",
                       inputSource            = "ABIF",
                       ABIF_Directory         = "./tmp/",
                       contigName             = "Achl_ACHLO006-09",
                       REGEX_SuffixForward    = "_[0-9]*_F.ab1",
                       REGEX_SuffixReverse    = "_[0-9]*_R.ab1",
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


In this example, :code:`contigName` is set to :code:`"Achl_ACHLO006-09"`, so only :code:`"Achl_ACHLO006-09_1_F.ab1"` and :code:`"Achl_ACHLO006-09_2_R.ab1"` will be selected to align to a contig.

The inputs of :code:`SangerContig` constructor function and :code:`new` method are same. For more details about *SangerContig* inputs and slots definition, please refer to `sangeranalyseR reference manual (need update) <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. The created *SangerContig* instance, :code:`sangerContig`, is used as the input for the following functions.

After the *SangerContig* instance is created, you can run :code:`sangerContig` to get basic information of the instance or run :code:`sangerContig@objectResults@readResultTable` to check the creation result of every Sanger read.

Here is the output of :code:`sangerContig`::

   SangerContig S4 instance
           Input Source :  ABIF 
         Process Method :  REGEX 
         ABIF Directory :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica/ACHLO 
   REGEX Suffix Forward :  _[0-9]*_F.ab1 
   REGEX Suffix Reverse :  _[0-9]*_R.ab1 
            Contig Name :  Achl_ACHLO006-09 
          'minReadsNum' :  2 
        'minReadLength' :  20 
      'minFractionCall' :  0.5 
      'maxFractionLost' :  0.5 
     'acceptStopCodons' :  TRUE 
         'readingFrame' :  1 
        Contig Sequence :  TTATATTTTATTCTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGGATTACGTCTTGAACGAATTCCCCTGTTTGTCTGAGCTGTGCTAATTACAGTTGTTCTTCTACTTCTATCTTTACCAGTGCTAGCAGGTGCCATTACCATACTTCTTACCGACCGAAACCTCAATACTTCATTCTTTGATCCTGCCGGTGGTGGAGACCCCATCCTC 
  Forward reads in the contig >>  1 
  Reverse reads in the contig >>  1 
   SUCCESS [2021-01-07 11:40:39] 'Achl_ACHLO006-09' is successfully created!


Here is the output of :code:`sangerContig@objectResults@readResultTable`::


                       readName  creationResult  errorType  errorMessage  inputSource     direction

   1   Achl_ACHLO006-09_1_F.ab1            TRUE       None          None         ABIF  Forward Read
   2   Achl_ACHLO006-09_2_R.ab1            TRUE       None          None         ABIF  Reverse Read


|


Updating *SangerContig* quality trimming parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++

In the previous :ref:`Creating *SangerContig* instance from **AB1**` part, the constructor function will apply the quality trimming parameters to all reads. After creating the SangerContig S4 instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will update all reads with the new trimming parameters and redo reads alignment. If users want to do quality trimming read by read instead of all at once, please read :ref:`Launching *SangerContig* Shiny app` page.

.. code-block:: R

   newSangerContig <- updateQualityParam(sangerContig,
                                         TrimmingMethod       = "M2",
                                         M1TrimmingCutoff     = NULL,
                                         M2CutoffQualityScore = 20,
                                         M2SlidingWindowSize  = 15)

|

Launching *SangerContig* Shiny app
+++++++++++++++++++++++++++++++++++

We create an interactive local Shiny app for users to go into each *SangerRead* in *SangerContig* instance. Users only need to run one function, :code:`launchApp`, with previously created instance as input and the *SangerContig* Shiny app will pop up. Here, we will go through *SangerRead* and *SangerContig* pages.

.. code-block:: R

  launchApp(newSangerContig)


*SangerContig* page (SC app)
-----------------------------
*SangerContig* page is the initial page of *SangerContig* Shiny app. :ref:`Figure 4<SangerContig_shiny_SangerContig_page>` shows the overview page of the contig. Notice that there is a red "Re-calculate Contig" button. Users need to click the button after changing the quality trimming parameters in order to get the updated information. In SangerContig page, there are two expendable tabs, “Forward Reads” and “Reverse Reads” storing the corresponding reads on the left-hand side navigation panel in :ref:`Figure 4<SangerContig_shiny_SangerContig_page>`. See :ref:`*SangerRead* page (SC app)` for more details of the subpage.

.. _SangerContig_shiny_SangerContig_page:
.. figure::  ../image/SangerContig_shiny_SangerContig_page.png
   :align:   center
   :scale:   20 %

   Figure 4. *SangerContig* Shiny app initial page - *SangerContig* page.

The information provided in this page are input parameters and contig results including “genetic code table”, “reference amino acid sequence”, “reads alignment”, “difference data frame”, “dendrogram”, “sample distance heatmap”, “indels data frame”, and “stop codons data frame”.

:ref:`Figure 5<SangerContig_shiny_alignment_differenceDF>` shows reads alignment result and difference data frame. The alignment is generated by :code:`AlignSeqs` or :code:`AlignTranslation` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package.


.. _SangerContig_shiny_alignment_differenceDF:
.. figure::  ../image/SangerContig_shiny_alignment_differenceDF.png
   :align:   center
   :scale:   30 %

   Figure 5. *SangerContig* page - reads alignment and difference data frame.

:ref:`Figure 6<SangerContig_shiny_dendrogram>` shows dendrogram result in both plot and in data frame. The results are generated by :code:`IdClusters` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package.

.. _SangerContig_shiny_dendrogram:
.. figure::  ../image/SangerContig_shiny_dendrogram.png
   :align:   center
   :scale:   30 %

   Figure 6. *SangerContig* page - dendrogram.

:ref:`Figure 7<SangerContig_shiny_samples_distance>` shows distance between **AB1** files. The results are generated by :code:`DistanceMatrix` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package. The heatmap is generated by :code:`plot_ly` function in `plotly <https://plot.ly/r/>`_ package.

.. _SangerContig_shiny_samples_distance:
.. figure::  ../image/SangerContig_shiny_samples_distance.png
   :align:   center
   :scale:   30 %

   Figure 7. *SangerContig* page - samples distance.

:ref:`Figure 8<SangerContig_shiny_indelsDF_stopcodonsDF>` shows insertions, deletions and stop codons data frame.

.. _SangerContig_shiny_indelsDF_stopcodonsDF:
.. figure::  ../image/SangerContig_shiny_indelsDF_stopcodonsDF.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - indels and stop codons data frame.


*SangerRead* page (SC app)
--------------------------

Now, let's go to the next level which is also the lowest level, *SangerRead* page. *SangerRead* page contains all details of a read including its trimming and chromatogram inputs and results. All reads are in "forward" or "reverse" direction. In this example, there is one read in each direction and :ref:`Figure 9<SangerContig_shiny_SangerRead_page>` shows "1 Forward Read" page. This page provides basic information, quality trimming inputs, chromatogram plotting inputs etc. Primary/secondary sequences and quality Phred scores table in this figure are dynamic based on the :code:`signalRatioCutoff` value for base calling and the length of them are always same. Another thing to mention is that primary/secondary sequences and the sequences in the chromatogram in :ref:`Figure 14<SangerContig_shiny_chromatogram_panel>` below will always be same after trimming and their color codings for A/T/C/G are same as well.

.. _SangerContig_shiny_SangerRead_page:
.. figure::  ../image/SangerContig_shiny_SangerRead_page.png
   :align:   center
   :scale:   20 %

   Figure 9. *SangerContig* Shiny app - *SangerRead* page

In quality trimming steps, we removes fragment at both ends of sequencing reads with low quality score. It is important because trimmed reads will improves alignment results. :ref:`Figure 10<SangerContig_shiny_trimming_1>` shows the UI for Trimming Method 1 (M1): ‘Modified Mott Trimming’. This method is implemented in `Phred <http://www.phrap.org/phredphrapconsed.html>`_. Users can change the cutoff score and click “Apply Trimming Parameters" button to update the UI. The value of input must be between 0 and 1. If the input is invalid, the cutoff score will be set to default 0.0001.

.. _SangerContig_shiny_trimming_1:
.. figure::  ../image/SangerContig_shiny_trimming_1.png
   :align:   center
   :scale:   30 %

   Figure 10. *SangerRead* page - Trimming Method 1 (M1): ‘Modified Mott Trimming’ UI.

:ref:`Figure 11<SangerContig_shiny_trimming_2>` shows another quality trimming method for users to choose from, Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’. This method is implemented in `Trimmomatics <http://www.usadellab.org/cms/?page=trimmomatic>`_. Users can change the cutoff quality score as well as sliding window size and click “Apply Trimming Parameters" button to update the UI. The value of cutoff quality score must be between 0 and 60 (default 20); the value of sliding window size must be between 0 and 40 (default 10). If the inputs are invalid, their values will be set to default.

.. _SangerContig_shiny_trimming_2:
.. figure::  ../image/SangerContig_shiny_trimming_2.png
   :align:   center
   :scale:   30 %

   Figure 11. *SangerRead* page - Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’ UI.

:ref:`Figure 12<SangerContig_shiny_trimmed_before_after>` shows the quality report before and after trimming. After clicking the “Apply Trimming Parameters” button in :ref:`Figure 10<SangerContig_shiny_trimming_1>` or :ref:`Figure 11<SangerContig_shiny_trimming_2>`, the values of these information boxes will be updated to the latest values.

.. _SangerContig_shiny_trimmed_before_after:
.. figure::  ../image/SangerContig_shiny_trimmed_before_after.png
   :align:   center
   :scale:   30 %

   Figure 12. *SangerRead* page - read quality report before / after trimming.

In :ref:`Figure 13<SangerContig_shiny_bp_quality_plot>`, the x-axis is the index of the base pairs; the y-axis is the Phred quality score. The green horizontal bar at the top of the plot is the raw read region and the orange horizontal bar represents the remaining read region. Both :ref:`Figure 13<SangerContig_shiny_bp_quality_plot>` trimming plot and :ref:`Figure 14<SangerContig_shiny_chromatogram_panel>` chromatogram will be updated once users change the quality trimming parameters and click the “Apply Trimming Parameters" button in :ref:`Figure 14<SangerContig_shiny_chromatogram_panel>`.

.. _SangerContig_shiny_bp_quality_plot:
.. figure::  ../image/SangerContig_shiny_bp_quality_plot.png
   :align:   center
   :scale:   30 %

   Figure 13. *SangerContig* page - quality trimming plot.

If we only see primary and secondary sequences in the table, we will loose some variations. Chromatogram is very helpful to check the peak resolution. :ref:`Figure 14<SangerContig_shiny_chromatogram_panel>` shows the panel of plotting chromatogram. Users can change four parameters: :code:`Base Number Per Row`, :code:`Height Per Row`, :code:`Signal Ratio Cutoff`, and :code:`Show Trimmed Region`. Among them, :code:`Signal Ratio Cutoff` is a key parameter. If its value is default value 0.33, it indicates that the lower peak should be at least 1/3rd as high as the higher peak for it count as a secondary peak.

.. _SangerContig_shiny_chromatogram_panel:
.. figure::  ../image/SangerContig_shiny_chromatogram_panel.png
   :align:   center
   :scale:   30 %

   Figure 14. *SangerContig* page - chromatogram panel.

Here is an example of applying new chromatogram parameters. We click “Show Trimmed Region” to set its value from :code:`FALSE` to :code:`TRUE` and click the "Apply Chromatogram Parameters" button. :ref:`Figure 15<SangerContig_plotting_popup>` shows the loading notification popup during base calling and chromatogram plotting.


.. _SangerContig_plotting_popup:
.. figure::  ../image/SangerContig_plotting_popup.png
   :align:   center
   :scale:   30 %

   Figure 15. *SangerContig* page - loading notification popup during replotting chromatogram.

After replotting the chromatogram, we can see that trimmed region is showed in red striped region. :ref:`Figure 16<SangerContig_shiny_chromatogram>` shows part of the the chromatogram (1 bp ~ 240 bp). Moreover, chromatogram will be replotted when trimmed positions or chromatogram parameters are updated.


.. _SangerContig_shiny_chromatogram:
.. figure::  ../image/SangerContig_shiny_chromatogram.png
   :align:   center
   :scale:   30 %

   Figure 16. *SangerContig* page - chromatogram with trimmed region showed.

To let users browse the trimmed primary/secondary sequences without finding “Trimming Start Point” and “Trimming End Point” by themselves, we provide the final trimmed primary/secondary sequences that will be used for reads alignment with quality scores in table format in :ref:`Figure 17<SangerContig_shiny_trimmed_sequences>`. Frameshift amino acid sequences are also provided.


.. _SangerContig_shiny_trimmed_sequences:
.. figure::  ../image/SangerContig_shiny_trimmed_sequences.png
   :align:   center
   :scale:   28 %

   Figure 17. *SangerContig* page - trimmed primary/secondary sequences and Phred quality score in table format.

We have updated the trimming and chromatogram parameters for each read. Now, we need to click “Re-calculate contig” button to do alignment again. Last but not least, we can save all data into a new ‘SangerContig’ S4 instance by clicking “Save S4 Instance button”. New S4 instance will be saved in **Rda** format. Users can run :code:`readRDS` function to load it into current R environment. :ref:`Figure 18<SangerContig_shiny_save_popup>` shows some hints in the save notification popup.

.. _SangerContig_shiny_save_popup:
.. figure::  ../image/SangerContig_shiny_save_popup.png
   :align:   center
   :scale:   25 %

   Figure 18. *SangerContig* page - saving notification popup.

|

Writing *SangerContig* FASTA files :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++++++
Users can write the *SangerContig* instance to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`reads_unalignment`: Writing reads into a single **FASTA** file (only trimmed without alignment).
* :code:`reads_alignment`: Writing reads alignment and contig read to a single **FASTA** file.
* :code:`contig`: Writing the contig to a single **FASTA** file.
* :code:`all`: Writing reads, reads alignment, and the contig into three different files.

Below is the oneliner for writing out **FASTA** files. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFasta` function.

.. code-block:: R

   writeFasta(newSangerContig,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA,
              selection         = "all")

Users can download the output FASTA file of this example through the following three links:

(1) :download:`Achl_ACHLO006-09_reads_unalignment.fa <../files/SangerContig_ab1/Achl_ACHLO006-09_reads_unalignment.fa>`
(2) :download:`Achl_ACHLO006-09_reads_alignment.fa <../files/SangerContig_ab1/Achl_ACHLO006-09_reads_alignment.fa>`
(3) :download:`Achl_ACHLO006-09_contig.fa <../files/SangerContig_ab1/Achl_ACHLO006-09_contig.fa>`


|

Generating *SangerContig* report :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save *SangerContig* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.


Users can set :code:`includeSangerRead` parameter to decide to which level the *SangerContig* report will go. Moreover, after the reports are generated,
users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it will take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReport(newSangerContig,
                  outputDir           = tempdir(),
                  includeSangerRead   = TRUE)

Users can access to '*Basic Information*', '*SangerContig Input Parameters*', '*Contig Sequence*' and '*Contig Results*' sections inside the generated `SangerContig html report of this example <https://howardchao.github.io/sangeranalyseR_report/SangerContig/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]/SangerContig_Report.html>`_. Furthermore, users can also navigate through html reports of all forward and reverse *SangerRead* in this *SangerContig* report.

-----

|
|


A Reproducible Example (*SangerContig*, **ab1**)
++++++++++++++++++++++++++++++++++++++++++++++++++


1. Preparing *SangerContig* **AB1** inputs
---------------------------------------------
The data of this example is in the sangeranalyseR package; thus, you can simply get its path from the library.

.. code-block:: R

    rawDataDir <- system.file("extdata", package = "sangeranalyseR")
    parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")

|

2. Creating *SangerContig* instance from **AB1**
----------------------------------------------------
Run the following on-liner to create the *SangerContig* instance.


.. code-block:: R

   # using `constructor` function to create SangerRead instance
   my_sangerContig <- SangerContig(inputSource           = "ABIF",
                                   processMethod         = "REGEX",
                                   ABIF_Directory        = parentDir,
                                   contigName            = "Achl_RBNII384-13",
                                   REGEX_SuffixForward   = "_[0-9]*_F.ab1",
                                   REGEX_SuffixReverse   = "_[0-9]*_R.ab1",
                                   refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")
   
   # using `new` method to create SangerRead instance
   my_sangerContig <- new("SangerContig",
                          inputSource           = "ABIF",
                          processMethod         = "REGEX",
                          ABIF_Directory        = parentDir,
                          contigName            = "Achl_RBNII384-13",
                          REGEX_SuffixForward   = "_[0-9]*_F.ab1",
                          REGEX_SuffixReverse   = "_[0-9]*_R.ab1",
                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")

.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

         INFO [2021-29-06 17:44:10] ========================================================
         INFO [2021-29-06 17:44:10] ================ Creating 'SangerContig' ===============
         INFO [2021-29-06 17:44:10] ========================================================
         INFO [2021-29-06 17:44:10]   >> Contig Name: 'Achl_RBNII384-13'
         INFO [2021-29-06 17:44:10]   >> You are using Regular Expression Method to group AB1 files!
         INFO [2021-29-06 17:44:10] >> Your contig name is Achl_RBNII384-13
         SUCCESS [2021-29-06 17:44:11] --------------------------------------------------------
         SUCCESS [2021-29-06 17:44:11] -------- 'SangerRead' S4 instance is created !! --------
         SUCCESS [2021-29-06 17:44:11] --------------------------------------------------------
         SUCCESS [2021-29-06 17:44:11]    >> 'Achl_RBNII384-13_1_F.ab1' is created (Forward Read; ABIF).
         SUCCESS [2021-29-06 17:44:12] --------------------------------------------------------
         SUCCESS [2021-29-06 17:44:12] -------- 'SangerRead' S4 instance is created !! --------
         SUCCESS [2021-29-06 17:44:12] --------------------------------------------------------
         SUCCESS [2021-29-06 17:44:12]    >> 'Achl_RBNII384-13_2_R.ab1' is created (Reverse Read; ABIF).
         INFO [2021-29-06 17:44:12]    >> The number of reads detected: 2
         INFO [2021-29-06 17:44:12] Correcting frameshifts in reads using amino acidreference sequence
         Assessing frameshifts in nucleotide sequences:
         |=============================================================================| 100%

         Time difference of 0.2 secs
         SUCCESS [2021-29-06 17:44:16] ==========================================================
         SUCCESS [2021-29-06 17:44:16] ======== 'SangerContig' S4 instance is created !! ========
         SUCCESS [2021-29-06 17:44:16] ==========================================================
         INFO [2021-29-06 17:44:16]    >> 2 read(s) created from ABIF file.
         INFO [2021-29-06 17:44:16]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
         INFO [2021-29-06 17:44:16]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
         INFO [2021-29-06 17:44:16]    >> Trimmed by 'M1 - Mott’s trimming algorithm'.
         DEBUG [2021-29-06 17:44:16]    >> For more information, please run 'object'
         DEBUG [2021-29-06 17:44:16]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads
         
|


3. Updating *SangerContig* quality trimming parameters
-------------------------------------------------------

.. code-block:: R

      newSangerContig <- updateQualityParam(my_sangerContig,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 20,
                                            M2SlidingWindowSize  = 15)

|


4. Launching *SangerContig* Shiny app
-----------------------------------------

Launch an interactive Shiny app to view your analysis, change the default settings, and so on.

.. code-block:: R

   launchApp(my_sangerContig)

|


5. Writing *SangerContig* FASTA files :sub:`(AB1)`
--------------------------------------------------

The following function can write the *SangerContig* instance into FASTA files. You just need to tell it where with the :code:`outputDir` argument.
.. code-block:: R

   writeFasta(my_sangerContig)


.. container:: toggle

     .. container:: header

        Following is the R shell output that you will get.

     .. code-block::

         INFO [2021-29-06 17:47:19] Your input is 'SangerContig' S4 instance
         INFO [2021-29-06 17:47:19] >>> outputDir : /private/var/folders/33/7v38jdjd2874jcxb6l71m00h0000gn/T/RtmpRAPaMV
         INFO [2021-29-06 17:47:19] Start to write 'Achl_RBNII384-13' to FASTA format ...
         INFO [2021-29-06 17:47:19] >> Writing alignment to FASTA ...
         INFO [2021-29-06 17:47:19] >> Writing all single reads to FASTA ...
         INFO [2021-29-06 17:47:19] >> Writing consensus read to FASTA ...
         INFO [2021-29-06 17:47:19] Finish writing 'Achl_RBNII384-13' to FASTA format

|

And you will get three FASTA files:

(1) :download:`Achl_RBNII384-13_reads_unalignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_unalignment.fa>`
(2) :download:`Achl_RBNII384-13_reads_alignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_alignment.fa>`
(3) :download:`Achl_RBNII384-13_contig.fa <../files/SangerContig_ab1/Achl_RBNII384-13_contig.fa>`

|

6. Generating *SangerContig* report :sub:`(AB1)`
-------------------------------------------------

Last but not least, generate an Rmarkdown report to store all the sequence information.

.. code-block:: R

   generateReport(my_sangerContig)

-----

|
|
