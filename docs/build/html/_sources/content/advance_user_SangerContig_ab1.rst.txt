Advanced User Guide - *SangerContig* (**AB1**)
==============================================

*SangerContig* is in the intermediate level of sangeranalyseR (:ref:`Figure_1<SangerContig_hierarchy>`), and each *SangerContig* instance corresponds to a contig in a Sanger sequencing experiment. Among its slots, there are two lists, forward and reverse read list, storing *SangerRead* in the corresponding direction. 

In this section, we are going to go through details about a reproducible *SangerContig* analysis example with the **AB1** file input in sangeranalyseR. By running the following example codes, you will get an end-to-end *SangerContig* analysis result. 


.. _SangerContig_hierarchy:
.. figure::  ../image/SangerContig_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerContig* level.

|

Preparing *SangerContig* **AB1** inputs
++++++++++++++++++++++++++++++++++++++++

The main input file format to create *SangerContig* instance is **AB1**. Before starting the analysis, users need to prepare one directory containing all **AB1** files, and all of them must be in the first layer of that directory. In other words, there should be no subdirectories. In this example, the data is in the sangeranalyseR package; thus, you can simply get its path by running the following codes:


.. code-block:: R

    rawDataDir <- system.file("extdata", package = "sangeranalyseR")
    parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")

The value of :code:`parentDir` is where all **AB1** files are placed. If your operating system is macOS, then its value should look like this:

.. code-block:: 

    /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica/RBNII

And we showed the files under :code:`parentDir` in :ref:`Figure_2<SangerContig_file_structure_ab1>`:

.. _SangerContig_file_structure_ab1:
.. figure::  ../image/SangerContig_file_structure.png
   :align:   center
   :scale:   60 %

   Figure 2. *SangerContig* filename regulation.

:ref:`Figure_2<SangerContig_file_structure_ab1>` shows the file-naming regulation and hierarchy. In this example, :code:`RBNII` is the parent directory, and all **AB1** files must be under its first layer. There are two ways for users to group their **AB1** files which are **"regular expression matching"** and **"CSV file matching"**, and following are instructions of how to prepare and name your **AB1** input files.


(1) "regular expression matching" *SangerContig* inputs (**AB1**)
-------------------------------------------------------------------
For regular expression matching method, sangeranalyseR will group **AB1** files based on their contig names and read directions in their filenames automatically; therefore, users have to follow the file-naming regulations below:

.. note::

    *  All input files must have **.ab1** as its file extension
    *  All input files must have the same contig name in their filename.
    *  Forward or reverse direction has to be specified in the filename.




There are four parameters, :code:`ABIF_Directory`, :code:`contigName`, :code:`REGEX_SuffixForward`, and :code:`REGEX_SuffixReverse`, that define the grouping rule to let sangeranalyseR automatically match correct **AB1** files and divide them into forward and reverse directions.

.. note::

  * :code:`ABIF_Directory`: this is the directory that contains all **AB1** files, and it can be either an absolute or relative path. We suggest users to put only target **AB1** files inside this directory and do not include any other unrelated files.

  * :code:`contigName`: this is a regular expression that matches filenames that are going to be included in the *SangerContig* analysis. :code:`grepl` function in R is used.

  * :code:`REGEX_SuffixForward`: this is a regular expression that matches all filenames in forward direction. :code:`grepl` function in R is used.

  * :code:`REGEX_SuffixReverse`: this is a regular expression that matches all filenames in reverse direction. :code:`grepl` function in R is used.

If you don't know what regular expression is, don't panic - it's just a way of recognising text. Please refer to :ref:`What is a regular expression?` for more details. Here is an example of how it works in sangeranalseR:

So how sangeranalyseR works is that it first matches the :code:`contigName` to exclude unrelated files and then separate the forward and reverse reads by matching :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. Therefore, it is important to make sure that all target **AB1** files share the same :code:`contigName` and carefully select your :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. The bad file-naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate wrong results. Therefore, it is important to have a consistent naming strategy. So, how should we systematically name **AB1** files? We suggest users to follow the file-naming regulation in :ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig>`. 

.. _sangeranalyseR_filename_convention_SangerContig:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file-naming regulation - *SangerContig*.


As you can see, the first part of the regulation is a consensus read name (or contig name), which narrows down the scope of **AB1** files to those we are going to examine. The second part of the regulation is an index. Since there might be more than one read that is in the forward or reverse direction, we recommend you to number your reads in the same contig group. The third part is a direction which is either 'F' (forward) or 'R' (reverse). Last but not least, files have to end with **.ab1** file extension.

To make it more specific, let's go back to the true example. In :ref:`Figure_2<SangerContig_file_structure_ab1>`, there are a lot of **AB1** files from different contigs in :code:`RBNII` (:code:`ABIF_Directory`). 
First, we set :code:`contigName` to :code:`"Achl_RBNII384-13"` to reduce candidates from eight to two **AB1** files, :code:`Achl_RBNII384-13_1_F.ab1` and :code:`Achl_RBNII384-13_2_R.ab1`. Then, we set :code:`REGEX_SuffixForward` to :code:`"_[0-9]*_F.ab1$"` and :code:`REGEX_SuffixReverse` to :code:`"_[0-9]*_R.ab1$"` to let sangeranalyseR match and group forward and reverse reads automatically. By the regular expression rule, :code:`Achl_RBNII384-13_1_F.ab1` and :code:`Achl_RBNII384-13_2_R.ab1` will be categorized into "forward read list" and "reverse read list" respectively. The reason why we strongly recommend you to follow this file-naming regulation is that by doing so, you can directly adopt the example regular expression matching values, :code:`"_[0-9]*_F.ab1$"` and :code:`"_[0-9]*_R.ab1$"`, to group reads and reduce chances of error. 

After understanding how parameters work, please refer to :ref:`Creating *SangerContig* instance from **AB1**` below to see how to create 'Achl_RBNII384-13' *SangerContig* instance.

(2) "CSV file matching" *SangerContig* inputs (**AB1**)
--------------------------------------------------------
For those who are not familiar with regular expression, we provide a second grouping approach, CSV file matching method. sangeranalyseR will group **AB1** files based on the information in a CSV file automatically; therefore, users have to follow the regulations below:

.. note::

    Here is an :download:`example CSV file <../files/SangerContig_ab1/names_conversion_2.csv>` (:ref:`Figure_4<sangeranalyseR_csv_file_SangerContig_ab1>`)

      .. _sangeranalyseR_csv_file_SangerContig_ab1:
      .. figure::  ../image/sangeranalyseR_csv_file_sangercontig_ab1.png
         :align:   center
         :scale:   70 %

         Figure 4. Example CSV file for *SangerContig* instance creation.  

    *  There must be three columns, "**reads**", "**direction**", and "**contig**", in the CSV file.
    *  The "**reads**" column stores the filename of **AB1** files that are going to be included in the analysis.
    *  The "**direction**" column stores the direction of the reads. It must be "F" (forward) or "R" (reverse).
    *  The "**contig**" column stores the contig name that each read blongs. Reads in the same contig have to have the same contig name, and they will be grouped into the same *SangerContig* instance.


There are three parameters, :code:`ABIF_Directory`, :code:`contigName`, and :code:`CSV_NamesConversion`,that define the grouping rule to help sangeranalseR to automatically match correct **AB1** files and divide them into forward and reverse directions.

.. note::

  * :code:`ABIF_Directory`: this is the directory that contains all **AB1** files, and it can be either an absolute or relative path. We suggest users to put only target AB1 files inside this directory and do not include any other unrelated files.
  * :code:`contigName`: this is a regular expression that matches filenames that are going to be included in the *SangerContig* analysis. :code:`grepl` function in R is used.
  * :code:`CSV_NamesConversion`: this is the path to the **CSV** file. It can be either an absolute or relative path.



The main difference between "CSV file matching" and "regular expression matching" is where the grouping rule is written. For "regular expression matching", rules are writtein in filenames, and thus there are more naming requirements for users to follow. In contrast, rules of "CSV file matching" are written in an additional **CSV** file so it is more flexible on **AB1** file-naming.


So how sangeranalyseR works is that it first reads in the **CSV** file (with *"reads"*, *"direction"*, and *"contig"* columns), filter out rows whose *"contig"* is not the value of :code:`contigName` parameter, find the names of **AB1** files listed in *"reads"*, and assign directions to them based on *"direction"*.

To make it more specific, let's go back to the true example. First, we prepare a :download:`CSV file <../files/SangerContig_ab1/names_conversion_2.csv>` (:code:`CSV_NamesConversion`) and a file directory like :ref:`Figure_2<SangerContig_file_structure_ab1>` (:code:`ABIF_Directory`) with some **AB1** files from different contigs. In the **CSV** file, both rows have the contig name :code:`"Achl_RBNII384-13"`, which is what we need to assign to the :code:`contigName` parameter. sangeranalyseR then checks and matches *"reads"* of these two rows, :code:`"Achl_RBNII384-13_1_F.ab1"` and :code:`"Achl_RBNII384-13_2_R.ab1"`, in :code:`RBNII` directory and reduce candidates from eight to two **AB1** files. Last, these two reads are assigned into "forward read list" and "reverse read list" respectively by the *"direction"* column.

After understanding how parameters work, please refer to :ref:`Creating *SangerContig* instance from **AB1**` below to see how to create 'Achl_RBNII384-13' *SangerContig* instance.

|


Creating *SangerContig* instance from **AB1**
++++++++++++++++++++++++++++++++++++++++++++++

After preparing the input directory, we can create a *SangerContig* instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Their input parameters are same, and all of them have their default values. For more details about *SangerContig* inputs and slots definition, please refer to `sangeranalyseR reference manual <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. We will explain two *SangerContig* instance creation methods, "regular expression matching" and "CSV file matching".




(1) "regular expression matching" *SangerContig* creation (**AB1**)
---------------------------------------------------------------------

The consturctor function and :code:`new` method below contain four parameters, :code:`ABIF_Directory`, :code:`contigName`, :code:`REGEX_SuffixForward`, and :code:`REGEX_SuffixReverse`, that we mentioned in the previous section. It also includes important parameters like quality trimming, chromatogram visualization, consensus alignment, and so on. Run the following code and create :code:`my_sangerContig` instance. 


.. code-block:: R
                     
   # using `constructor` function to create SangerContig instance
   my_sangerContig <- SangerContig(inputSource           = "ABIF",
                                   processMethod         = "REGEX",
                                   ABIF_Directory        = parentDir,
                                   contigName            = "Achl_RBNII384-13",
                                   REGEX_SuffixForward   = "_[0-9]*_F.ab1$",
                                   REGEX_SuffixReverse   = "_[0-9]*_R.ab1$",
                                   TrimmingMethod        = "M1",
                                   M1TrimmingCutoff      = 0.0001,
                                   M2CutoffQualityScore  = NULL,
                                   M2SlidingWindowSize   = NULL,
                                   baseNumPerRow         = 100,
                                   heightPerRow          = 200,
                                   signalRatioCutoff     = 0.33,
                                   showTrimmed           = TRUE,
                                   refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                   minReadsNum           = 2,
                                   minReadLength         = 20,
                                   minFractionCall       = 0.5,
                                   maxFractionLost       = 0.5,
                                   geneticCode           = GENETIC_CODE,
                                   acceptStopCodons      = TRUE,
                                   readingFrame          = 1,
                                   processorsNum         = 1)

   # using `new` method to create SangerContig instance
   my_sangerContig <- new("SangerContig",
                          inputSource           = "ABIF",
                          processMethod         = "REGEX",
                          ABIF_Directory        = parentDir,
                          contigName            = "Achl_RBNII384-13",
                          REGEX_SuffixForward   = "_[0-9]*_F.ab1$",
                          REGEX_SuffixReverse   = "_[0-9]*_R.ab1$",
                          TrimmingMethod        = "M1",
                          M1TrimmingCutoff      = 0.0001,
                          M2CutoffQualityScore  = NULL,
                          M2SlidingWindowSize   = NULL,
                          baseNumPerRow         = 100,
                          heightPerRow          = 200,
                          signalRatioCutoff     = 0.33,
                          showTrimmed           = TRUE,
                          refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                          minReadsNum           = 2,
                          minReadLength         = 20,
                          minFractionCall       = 0.5,
                          maxFractionLost       = 0.5,
                          geneticCode           = GENETIC_CODE,
                          acceptStopCodons      = TRUE,
                          readingFrame          = 1,
                          processorsNum         = 1)



In this example, :code:`contigName` is set to :code:`Achl_RBNII384-13`, so only :code:`Achl_RBNII384-13_1_F.ab1` and :code:`Achl_RBNII384-13_2_R.ab1` are selected. Moreover, by regular expression pattern matching, :code:`Achl_RBNII384-13_1_F.ab1` is categorized into the forward list, and :code:`Achl_RBNII384-13_2_R.ab1` is categorized into the reverse read. Both reads are aligned into a contig, :code:`my_sangerContig`, and it will be used as the input for the following functions.

Inside the R shell, you can run :code:`my_sangerContig` to get basic information of the instance or run :code:`my_sangerContig@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerContig` is successfully created.

Here is the output of :code:`my_sangerContig`::

   SangerContig S4 instance
            Input Source :  ABIF 
            Process Method :  REGEX 
            ABIF Directory :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica/RBNII 
      REGEX Suffix Forward :  _[0-9]*_F.ab1$ 
      REGEX Suffix Reverse :  _[0-9]*_R.ab1$ 
               Contig Name :  Achl_RBNII384-13 
            'minReadsNum' :  2 
         'minReadLength' :  20 
         'minFractionCall' :  0.5 
         'maxFractionLost' :  0.5 
      'acceptStopCodons' :  TRUE 
            'readingFrame' :  1 
         Contig Sequence :  AGCAGGATAGTAGGGGCTGGTATAAGACTCCTAATTCGAATTGAGCTAAGACAGCCGGGAGCATTTCTAGGAAGGGATCAACTCTATAACACTATTGTAACTGCTCACGCATTTGTAATAATTTTCTTTCTAGTAATACCTGTATTTATTGGGGGGTTCGGTAATTGACTTCTACCTTTAATACTTGGAGCCCCTGACATGGCATTCCCACGTCTTAACAACATAAGATTTTGACTCCTTCCCCCATCACTAATCCTTCTAGTATCCTCTGCTGCAGTAGAAAAGGGTGCGGGAACTGGATGAACTGTTTATCCACCCCTAGCAAGAAACATTGCTCATGCCGGCCCATCTGTAGACTTAGCTATTTTTTCTCTTCATTTAGCAGGTGCTTCATCAATCTTGGGTGCCATTAATTTTATTACTACTGTTATTAACATACGATGAAGAGGCTTACGACTTGAACGAATCCCATTATTCGTTTGAGCCGTACTAATTACAGTGGTCCTTCTACTCTTATCTTTACCAGTATTAGCCGGTGCAATTACTATACTACTTACCGATCGAAATCTAAATACCTCCTTCTTTGACCCTGCTGGAGGCGGAGAT 
   Forward reads in the contig >>  1 
   Reverse reads in the contig >>  1 
   SUCCESS [2021-12-07 17:01:18] 'Achl_RBNII384-13' is successfully created!

Here is the output of :code:`my_sangerContig@objectResults@readResultTable`::


                        readName creationResult errorType errorMessage inputSource    direction
      1 Achl_RBNII384-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
      2 Achl_RBNII384-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read




(2) "CSV file matching" *SangerContig* creation (**AB1**)
----------------------------------------------------------
The consturctor function and :code:`new` method below contain three parameters, :code:`ABIF_Directory`, :code:`contigName`, and :code:`CSV_NamesConversion`, that we mentioned in the previous section. It also includes important parameters like quality trimming, chromatogram visualization, consensus alignment, and so on. Run the following code and create :code:`my_sangerContig` instance. 

.. code-block:: R

      csv_namesConversion <- file.path(rawDataDir, "ab1", "SangerContig", "names_conversion_2.csv")

      # using `constructor` function to create SangerContig instance
      my_sangerContig <- SangerContig(inputSource            = "ABIF",
                                      processMethod          = "CSV",
                                      ABIF_Directory         = parentDir,
                                      contigName             = "Achl_RBNII384-13",
                                      CSV_NamesConversion    = csv_namesConversion,
                                      TrimmingMethod         = "M1",
                                      M1TrimmingCutoff       = 0.0001,
                                      M2CutoffQualityScore   = NULL,
                                      M2SlidingWindowSize    = NULL,
                                      baseNumPerRow          = 100,
                                      heightPerRow           = 200,
                                      signalRatioCutoff      = 0.33,
                                      showTrimmed            = TRUE,
                                      refAminoAcidSeq        = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                      minReadsNum            = 2,
                                      minReadLength          = 20,
                                      minFractionCall        = 0.5,
                                      maxFractionLost        = 0.5,
                                      geneticCode            = GENETIC_CODE,
                                      acceptStopCodons       = TRUE,
                                      readingFrame           = 1,
                                      processorsNum          = 1)


      # using `new` method to create SangerContig instance
      my_sangerContig <- new("SangerContig",
                             inputSource           = "ABIF",
                             processMethod         = "CSV",
                             ABIF_Directory        = parentDir,
                             contigName            = "Achl_RBNII384-13",
                             CSV_NamesConversion   = csv_namesConversion,
                             TrimmingMethod         = "M1",
                             M1TrimmingCutoff       = 0.0001,
                             M2CutoffQualityScore   = NULL,
                             M2SlidingWindowSize    = NULL,
                             baseNumPerRow          = 100,
                             heightPerRow           = 200,
                             signalRatioCutoff      = 0.33,
                             showTrimmed            = TRUE,
                             refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                             minReadsNum            = 2,
                             minReadLength          = 20,
                             minFractionCall        = 0.5,
                             maxFractionLost        = 0.5,
                             geneticCode            = GENETIC_CODE,
                             acceptStopCodons       = TRUE,
                             readingFrame           = 1,
                             processorsNum          = 1)


First, you need to load the **CSV** file into the R environment. If you are still don't know how to prepare it, please check :ref:`(2) "CSV file matching" *SangerContig* inputs (**AB1**)`. Then, it will follow rules in the CSV file and create :code:`my_sangerContig`. After it's created, inside the R shell, you can run :code:`my_sangerContig` to get basic information of the instance or run :code:`my_sangerContig@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerContig` is successfully created.

Here is the output of :code:`my_sangerContig`::

   SangerContig S4 instance
            Input Source :  ABIF 
            Process Method :  CSV 
            ABIF Directory :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica/RBNII 
      CSV Names Conversion :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv 
               Contig Name :  Achl_RBNII384-13 
            'minReadsNum' :  2 
         'minReadLength' :  20 
         'minFractionCall' :  0.5 
         'maxFractionLost' :  0.5 
      'acceptStopCodons' :  TRUE 
            'readingFrame' :  1 
         Contig Sequence :  AGCAGGATAGTAGGGGCTGGTATAAGACTCCTAATTCGAATTGAGCTAAGACAGCCGGGAGCATTTCTAGGAAGGGATCAACTCTATAACACTATTGTAACTGCTCACGCATTTGTAATAATTTTCTTTCTAGTAATACCTGTATTTATTGGGGGGTTCGGTAATTGACTTCTACCTTTAATACTTGGAGCCCCTGACATGGCATTCCCACGTCTTAACAACATAAGATTTTGACTCCTTCCCCCATCACTAATCCTTCTAGTATCCTCTGCTGCAGTAGAAAAGGGTGCGGGAACTGGATGAACTGTTTATCCACCCCTAGCAAGAAACATTGCTCATGCCGGCCCATCTGTAGACTTAGCTATTTTTTCTCTTCATTTAGCAGGTGCTTCATCAATCTTGGGTGCCATTAATTTTATTACTACTGTTATTAACATACGATGAAGAGGCTTACGACTTGAACGAATCCCATTATTCGTTTGAGCCGTACTAATTACAGTGGTCCTTCTACTCTTATCTTTACCAGTATTAGCCGGTGCAATTACTATACTACTTACCGATCGAAATCTAAATACCTCCTTCTTTGACCCTGCTGGAGGCGGAGAT 
   Forward reads in the contig >>  1 
   Reverse reads in the contig >>  1 
   SUCCESS [2021-12-07 17:11:48] 'Achl_RBNII384-13' is successfully created!


Here is the output of :code:`my_sangerContig@objectResults@readResultTable`::


                     readName creationResult errorType errorMessage inputSource    direction
   1 Achl_RBNII384-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   2 Achl_RBNII384-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read



|


Updating *SangerContig* quality trimming parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++

In the previous :ref:`Creating *SangerContig* instance from **AB1**` part, the constructor function will apply the quality trimming parameters to all reads. After creating a *SangerContig* instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will update all reads with the new trimming parameters and redo reads alignment. If users want to do quality trimming read by read instead of all at once, please move on to the next section, :ref:`Launching *SangerContig* Shiny app` page.

.. code-block:: R

   newSangerContig <- updateQualityParam(my_sangerContig,
                                         TrimmingMethod       = "M2",
                                         M1TrimmingCutoff     = NULL,
                                         M2CutoffQualityScore = 20,
                                         M2SlidingWindowSize  = 15)

|

Launching *SangerContig* Shiny app
+++++++++++++++++++++++++++++++++++

We create an interactive local Shiny app for users to go into each *SangerRead* in *SangerContig* instance. Users only need to run one function, :code:`launchApp`, with previously created instance as input and the *SangerContig* Shiny app will pop up. Here, we will go through *SangerRead* and *SangerContig* pages.

.. code-block:: R

  launchApp(my_sangerContig)


*SangerContig* page (SC app)
-----------------------------
*SangerContig* page is the initial page of *SangerContig* Shiny app. :ref:`Figure 5<SangerContig_shiny_SangerContig_page>` shows the overview page of the contig. Notice that there is a red "Re-calculate Contig" button. Users need to click the button after changing the quality trimming parameters in order to get the updated information. In SangerContig page, there are two expendable tabs, “Forward Reads” and “Reverse Reads” storing the corresponding reads on the left-hand side navigation panel in :ref:`Figure 5<SangerContig_shiny_SangerContig_page>`. See :ref:`*SangerRead* page (SC app)` for more details of the subpage.

.. _SangerContig_shiny_SangerContig_page:
.. figure::  ../image/SangerContig_shiny_SangerContig_page.png
   :align:   center
   :scale:   20 %

   Figure 5. *SangerContig* Shiny app initial page - *SangerContig* page.

The information provided in this page are input parameters and contig results including “genetic code table”, “reference amino acid sequence”, “reads alignment”, “difference data frame”, “dendrogram”, “sample distance heatmap”, “indels data frame”, and “stop codons data frame”.

:ref:`Figure 6<SangerContig_shiny_alignment_differenceDF>` shows reads alignment result and difference data frame. The alignment is generated by :code:`AlignSeqs` or :code:`AlignTranslation` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package.


.. _SangerContig_shiny_alignment_differenceDF:
.. figure::  ../image/SangerContig_shiny_alignment_differenceDF.png
   :align:   center
   :scale:   30 %

   Figure 6. *SangerContig* page - reads alignment and difference data frame.

:ref:`Figure 7<SangerContig_shiny_dendrogram>` shows dendrogram result in both plot and in data frame. The results are generated by :code:`IdClusters` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package.

.. _SangerContig_shiny_dendrogram:
.. figure::  ../image/SangerContig_shiny_dendrogram.png
   :align:   center
   :scale:   30 %

   Figure 7. *SangerContig* page - dendrogram.

:ref:`Figure 8<SangerContig_shiny_samples_distance>` shows distance between **AB1** files. The results are generated by :code:`DistanceMatrix` function in `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ package. The heatmap is generated by :code:`plot_ly` function in `plotly <https://plot.ly/r/>`_ package.

.. _SangerContig_shiny_samples_distance:
.. figure::  ../image/SangerContig_shiny_samples_distance.png
   :align:   center
   :scale:   30 %

   Figure 8. *SangerContig* page - samples distance.

:ref:`Figure 9<SangerContig_shiny_indelsDF_stopcodonsDF>` shows insertions, deletions and stop codons data frame.

.. _SangerContig_shiny_indelsDF_stopcodonsDF:
.. figure::  ../image/SangerContig_shiny_indelsDF_stopcodonsDF.png
   :align:   center
   :scale:   30 %

   Figure 9. *SangerContig* page - indels and stop codons data frame.


*SangerRead* page (SC app)
--------------------------

Now, let's go to the next level which is also the lowest level, *SangerRead* page. *SangerRead* page contains all details of a read including its trimming and chromatogram inputs and results. All reads are in "forward" or "reverse" direction. In this example, there is one read in each direction and :ref:`Figure 10<SangerContig_shiny_SangerRead_page>` shows "1 Forward Read" page. This page provides basic information, quality trimming inputs, chromatogram plotting inputs etc. Primary/secondary sequences and quality Phred scores table in this figure are dynamic based on the :code:`signalRatioCutoff` value for base calling and the length of them are always same. Another thing to mention is that primary/secondary sequences and the sequences in the chromatogram in :ref:`Figure 15<SangerContig_shiny_chromatogram_panel>` below will always be same after trimming and their color codings for A/T/C/G are same as well.

.. _SangerContig_shiny_SangerRead_page:
.. figure::  ../image/SangerContig_shiny_SangerRead_page.png
   :align:   center
   :scale:   20 %

   Figure 10. *SangerContig* Shiny app - *SangerRead* page

In quality trimming steps, we removes fragment at both ends of sequencing reads with low quality score. It is important because trimmed reads will improves alignment results. :ref:`Figure 11<SangerContig_shiny_trimming_1>` shows the UI for Trimming Method 1 (M1): ‘Modified Mott Trimming’. This method is implemented in `Phred <http://www.phrap.org/phredphrapconsed.html>`_. Users can change the cutoff score and click “Apply Trimming Parameters" button to update the UI. The value of input must be between 0 and 1. If the input is invalid, the cutoff score will be set to default 0.0001.

.. _SangerContig_shiny_trimming_1:
.. figure::  ../image/SangerContig_shiny_trimming_1.png
   :align:   center
   :scale:   30 %

   Figure 11. *SangerRead* page - Trimming Method 1 (M1): ‘Modified Mott Trimming’ UI.

:ref:`Figure 12<SangerContig_shiny_trimming_2>` shows another quality trimming method for users to choose from, Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’. This method is implemented in `Trimmomatics <http://www.usadellab.org/cms/?page=trimmomatic>`_. Users can change the cutoff quality score as well as sliding window size and click “Apply Trimming Parameters" button to update the UI. The value of cutoff quality score must be between 0 and 60 (default 20); the value of sliding window size must be between 0 and 40 (default 10). If the inputs are invalid, their values will be set to default.

.. _SangerContig_shiny_trimming_2:
.. figure::  ../image/SangerContig_shiny_trimming_2.png
   :align:   center
   :scale:   30 %

   Figure 12. *SangerRead* page - Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’ UI.

:ref:`Figure 13<SangerContig_shiny_trimmed_before_after>` shows the quality report before and after trimming. After clicking the “Apply Trimming Parameters” button in :ref:`Figure 11<SangerContig_shiny_trimming_1>` or :ref:`Figure 12<SangerContig_shiny_trimming_2>`, the values of these information boxes will be updated to the latest values.

.. _SangerContig_shiny_trimmed_before_after:
.. figure::  ../image/SangerContig_shiny_trimmed_before_after.png
   :align:   center
   :scale:   30 %

   Figure 13. *SangerRead* page - read quality report before / after trimming.

In :ref:`Figure 14<SangerContig_shiny_bp_quality_plot>`, the x-axis is the index of the base pairs; the y-axis is the Phred quality score. The green horizontal bar at the top of the plot is the raw read region and the orange horizontal bar represents the remaining read region. Both :ref:`Figure 14<SangerContig_shiny_bp_quality_plot>` trimming plot and :ref:`Figure 15<SangerContig_shiny_chromatogram_panel>` chromatogram will be updated once users change the quality trimming parameters and click the “Apply Trimming Parameters" button in :ref:`Figure 15<SangerContig_shiny_chromatogram_panel>`.

.. _SangerContig_shiny_bp_quality_plot:
.. figure::  ../image/SangerContig_shiny_bp_quality_plot.png
   :align:   center
   :scale:   30 %

   Figure 14. *SangerContig* page - quality trimming plot.

If we only see primary and secondary sequences in the table, we will loose some variations. Chromatogram is very helpful to check the peak resolution. :ref:`Figure 15<SangerContig_shiny_chromatogram_panel>` shows the panel of plotting chromatogram. Users can change four parameters: :code:`Base Number Per Row`, :code:`Height Per Row`, :code:`Signal Ratio Cutoff`, and :code:`Show Trimmed Region`. Among them, :code:`Signal Ratio Cutoff` is a key parameter. If its value is default value 0.33, it indicates that the lower peak should be at least 1/3rd as high as the higher peak for it count as a secondary peak.

.. _SangerContig_shiny_chromatogram_panel:
.. figure::  ../image/SangerContig_shiny_chromatogram_panel.png
   :align:   center
   :scale:   30 %

   Figure 15. *SangerContig* page - chromatogram panel.

Here is an example of applying new chromatogram parameters. We click “Show Trimmed Region” to set its value from :code:`FALSE` to :code:`TRUE` and click the "Apply Chromatogram Parameters" button. :ref:`Figure 16<SangerContig_plotting_popup>` shows the loading notification popup during base calling and chromatogram plotting.


.. _SangerContig_plotting_popup:
.. figure::  ../image/SangerContig_plotting_popup.png
   :align:   center
   :scale:   30 %

   Figure 16. *SangerContig* page - loading notification popup during replotting chromatogram.

After replotting the chromatogram, we can see that trimmed region is showed in red striped region. :ref:`Figure 17<SangerContig_shiny_chromatogram>` shows part of the the chromatogram (1 bp ~ 240 bp). Moreover, chromatogram will be replotted when trimmed positions or chromatogram parameters are updated.


.. _SangerContig_shiny_chromatogram:
.. figure::  ../image/SangerContig_shiny_chromatogram.png
   :align:   center
   :scale:   30 %

   Figure 17. *SangerContig* page - chromatogram with trimmed region showed.

To let users browse the trimmed primary/secondary sequences without finding “Trimming Start Point” and “Trimming End Point” by themselves, we provide the final trimmed primary/secondary sequences that will be used for reads alignment with quality scores in table format in :ref:`Figure 18<SangerContig_shiny_trimmed_sequences>`. Frameshift amino acid sequences are also provided.


.. _SangerContig_shiny_trimmed_sequences:
.. figure::  ../image/SangerContig_shiny_trimmed_sequences.png
   :align:   center
   :scale:   28 %

   Figure 18. *SangerContig* page - trimmed primary/secondary sequences and Phred quality score in table format.

We have updated the trimming and chromatogram parameters for each read. Now, we need to click “Re-calculate contig” button to do alignment again. Last but not least, we can save all data into a new ‘SangerContig’ S4 instance by clicking “Save S4 Instance button”. New S4 instance will be saved in **Rda** format. Users can run :code:`readRDS` function to load it into current R environment. :ref:`Figure 19<SangerContig_shiny_save_popup>` shows some hints in the save notification popup.

.. _SangerContig_shiny_save_popup:
.. figure::  ../image/SangerContig_shiny_save_popup.png
   :align:   center
   :scale:   25 %

   Figure 19. *SangerContig* page - saving notification popup.

|

Writing *SangerContig* FASTA files :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++++++
Users can write the *SangerContig* instance, :code:`my_sangerContig`, to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`reads_unalignment`: Writing reads into a single **FASTA** file (only trimmed without alignment).
* :code:`reads_alignment`: Writing reads alignment and contig read to a single **FASTA** file.
* :code:`contig`: Writing the contig to a single **FASTA** file.
* :code:`all`: Writing reads, reads alignment, and the contig into three different files.

Below is the oneliner for writing out **FASTA** files. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFasta` function.

.. code-block:: R

   writeFasta(my_sangerContig,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA,
              selection         = "all")

Users can download the output FASTA file of this example through the following three links:

(1) :download:`Achl_RBNII384-13_reads_unalignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_unalignment.fa>`
(2) :download:`Achl_RBNII384-13_reads_alignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_alignment.fa>`
(3) :download:`Achl_RBNII384-13_contig.fa <../files/SangerContig_ab1/Achl_RBNII384-13_contig.fa>`


|

Generating *SangerContig* report :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save *SangerContig* instance, :code:`my_sangerContig`, into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.

Users can set :code:`includeSangerRead` parameter to decide to which level the *SangerContig* report will go. Moreover, after the reports are generated,
users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it will take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReport(my_sangerContig,
                  outputDir           = tempdir(),
                  includeSangerRead   = FALSE)

Users can access to '*Basic Information*', '*SangerContig Input Parameters*', '*Contig Sequence*' and '*Contig Results*' sections inside the generated `SangerContig html report of this example <https://howardchao.github.io/sangeranalyseR_report/SangerContig/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]/SangerContig_Report.html>`_. Furthermore, users can also navigate through html reports of all forward and reverse *SangerRead* in this *SangerContig* report.

-----

|
|


Code summary (*SangerContig*, **AB1**)
++++++++++++++++++++++++++++++++++++++++++++++++++


(1) Preparing *SangerContig* **AB1** inputs
---------------------------------------------

.. code-block:: R

    rawDataDir <- system.file("extdata", package = "sangeranalyseR")
    parentDir <- file.path(rawDataDir, "Allolobophora_chlorotica", "RBNII")

|

(2) Creating *SangerContig* instance from **AB1**
----------------------------------------------------

(2.1) "Regular Expression Method" *SangerContig* creation (**AB1**)
***********************************************************************

.. code-block:: R

   # using `constructor` function to create SangerContig instance
   my_sangerContig <- SangerContig(inputSource           = "ABIF",
                                   processMethod         = "REGEX",
                                   ABIF_Directory        = parentDir,
                                   contigName            = "Achl_RBNII384-13",
                                   REGEX_SuffixForward   = "_[0-9]*_F.ab1$",
                                   REGEX_SuffixReverse   = "_[0-9]*_R.ab1$",
                                   refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")

   # using `new` method to create SangerContig instance
   my_sangerContig <- new("SangerContig",
                          inputSource           = "ABIF",
                          processMethod         = "REGEX",
                          ABIF_Directory        = parentDir,
                          contigName            = "Achl_RBNII384-13",
                          REGEX_SuffixForward   = "_[0-9]*_F.ab1$",
                          REGEX_SuffixReverse   = "_[0-9]*_R.ab1$",
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

(2.2) "CSV file matching" *SangerContig* creation (**AB1**)
*************************************************************


.. code-block:: R

   csv_namesConversion <- file.path(rawDataDir, "ab1", "SangerContig", "names_conversion_2.csv")

   # using `constructor` function to create SangerContig instance
   my_sangerContig <- SangerContig(inputSource           = "ABIF",
                                   processMethod         = "CSV",
                                   ABIF_Directory        = parentDir,
                                   contigName            = "Achl_RBNII384-13",
                                   CSV_NamesConversion   = csv_namesConversion,
                                   refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


   # using `new` method to create SangerContig instance
   my_sangerContig <- new("SangerContig",
                          inputSource           = "ABIF",
                          processMethod         = "CSV",
                          ABIF_Directory        = parentDir,
                          contigName            = "Achl_RBNII384-13",
                          CSV_NamesConversion   = csv_namesConversion,
                          refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

      WARN [2021-09-07 01:47:45] 'Achl_RBNII395-13_1_F.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      WARN [2021-09-07 01:47:45] 'Achl_RBNII395-13_2_R.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      WARN [2021-09-07 01:47:45] 'Achl_RBNII396-13_1_F.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      WARN [2021-09-07 01:47:45] 'Achl_RBNII396-13_2_R.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      WARN [2021-09-07 01:47:45] 'Achl_RBNII397-13_1_F.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      WARN [2021-09-07 01:47:45] 'Achl_RBNII397-13_2_R.ab1' is not in the csv file (/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerContig/names_conversion_2.csv)
      INFO [2021-09-07 01:47:45] ========================================================
      INFO [2021-09-07 01:47:45] ================ Creating 'SangerContig' ===============
      INFO [2021-09-07 01:47:45] ========================================================
      INFO [2021-09-07 01:47:45]   >> Contig Name: 'Achl_RBNII384-13'
      INFO [2021-09-07 01:47:45]   >> You are using CSV Name Conversion Method to group AB1 files!
      INFO [2021-09-07 01:47:45] >> Your contig name is Achl_RBNII384-13
      SUCCESS [2021-09-07 01:47:46] --------------------------------------------------------
      SUCCESS [2021-09-07 01:47:46] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-09-07 01:47:46] --------------------------------------------------------
      SUCCESS [2021-09-07 01:47:46]    >> 'Achl_RBNII384-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-09-07 01:47:47] --------------------------------------------------------
      SUCCESS [2021-09-07 01:47:47] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-09-07 01:47:47] --------------------------------------------------------
      SUCCESS [2021-09-07 01:47:47]    >> 'Achl_RBNII384-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-09-07 01:47:47]    >> The number of reads detected: 2
      INFO [2021-09-07 01:47:47] Correcting frameshifts in reads using amino acidreference sequence
      Assessing frameshifts in nucleotide sequences:
        |================================================================================================| 100%

      Time difference of 0.08 secs
      SUCCESS [2021-09-07 01:47:51] ==========================================================
      SUCCESS [2021-09-07 01:47:51] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-09-07 01:47:51] ==========================================================
      INFO [2021-09-07 01:47:51]    >> 2 read(s) created from ABIF file.
      INFO [2021-09-07 01:47:51]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-09-07 01:47:51]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-09-07 01:47:51]    >> Trimmed by 'M1 - Mott’s trimming algorithm'.
      DEBUG [2021-09-07 01:47:51]    >> For more information, please run 'object'
      DEBUG [2021-09-07 01:47:51]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads

|


(3) Updating *SangerContig* quality trimming parameters
-------------------------------------------------------

.. code-block:: R

      newSangerContig <- updateQualityParam(my_sangerContig,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 20,
                                            M2SlidingWindowSize  = 15)

|


(4) Launching *SangerContig* Shiny app
-----------------------------------------

.. code-block:: R

   launchApp(my_sangerContig)

|


(5) Writing *SangerContig* FASTA files :sub:`(AB1)`
-----------------------------------------------------

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

You will get three FASTA files:

(1) :download:`Achl_RBNII384-13_reads_unalignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_unalignment.fa>`
(2) :download:`Achl_RBNII384-13_reads_alignment.fa <../files/SangerContig_ab1/Achl_RBNII384-13_reads_alignment.fa>`
(3) :download:`Achl_RBNII384-13_contig.fa <../files/SangerContig_ab1/Achl_RBNII384-13_contig.fa>`

|

(6) Generating *SangerContig* report :sub:`(AB1)`
-------------------------------------------------

.. code-block:: R

   generateReport(my_sangerContig)

-----

|
|
