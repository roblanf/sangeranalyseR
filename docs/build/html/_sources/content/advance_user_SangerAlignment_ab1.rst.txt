Advanced User Guide - *SangerAlignment* (**AB1**)
=================================================
*SangerAlignment* is in the toppest level of sangeranalyseR (:ref:`Figure_1<SangerAlignment_hierachy>`), and each **SangerAlignment** instance corresponds to an alignment of contigs in a Sanger sequencing experiment. Among its slots, there is a *SangerContig* list which will be aligned into a consensus contig. Users can access to each *SangerContig* and *SangerRead* inside a *SangerAlignment* instance. 


In this section, we are going to go through details about a reproducible *SangerAlignment* analysis example with the **AB1** file input in sangeranalyseR. By running the following example codes, you will get an end-to-end *SangerAlignment* analysis result.

.. _SangerAlignment_hierachy:
.. figure::  ../image/SangerAlignment_hierachy.png
   :align:   center
   :scale:   20 %

   Figure 1. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

|

Preparing *SangerAlignment* **AB1** input
+++++++++++++++++++++++++++++++++++++++++++++

The main input file format to create *SangerAlignment* instance is **AB1**. Before starting the analysis, users need to prepare one directory containing all **AB1** files, and they can be either all placed in the first layer of that directory or be distributed in different subdirectories. In this example, the data are in the sangeranalyseR package; thus, you can simply get its path by running the following codes:


.. code-block:: R

   rawDataDir <- system.file("extdata", package = "sangeranalyseR")
   parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica')

The value of :code:`parentDir` is where all **AB1** files are placed. If your operating system is macOS, then its value should look like this:

.. code-block:: 

   /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica/

And we showed the files under :code:`parentDir` in :ref:`Figure_2<SangerAlignment_file_structure_ab1>`:

.. _SangerAlignment_file_structure_ab1:
.. figure::  ../image/SangerAlignment_file_structure.png
   :align:   center
   :scale:   50 %

   Figure 2. *SangerAlignment* filename regulation.

:ref:`Figure_2<SangerAlignment_file_structure_ab1>` shows the file-naming regulation and hierarchy. In this example, :code:`Allolobophora_chlorotica` is the parent directory, and **AB1** files are separated into :code:`ACHLO` and :code:`RBNII` directories. There are two ways for users to group their **AB1** files which are **"regular expression matching"** and **"CSV file matching"**, and following are instructions of how to prepare and name your **AB1** input files.



(1) "regular expression matching" *SangerAlignment* inputs (**AB1**)
---------------------------------------------------------------------
For regular expression matching method, sangeranalyseR will group **AB1** files based on their contig names and read directions in their filenames automatically; therefore, users have to follow the file-naming regulations below:

.. note::

    *  All input files must have **.ab1** as its file extension.
    *  Input files that are in the same contig group must have the same contig name in their filenames.
    *  Forward or reverse direction has to be specified in the filename.

There are three parameters, :code:`ABIF_Directory`, :code:`REGEX_SuffixForward`, and :code:`REGEX_SuffixReverse`, that define the grouping rule to let sangeranalyseR automatically match correct **AB1** files and divide them into forward and reverse directions.

.. note::

  * :code:`ABIF_Directory`: this is the directory that contains all **AB1** files, and it can be either an absolute or relative path. We suggest users to put only target **AB1** files inside this directory and do not include any other unrelated files.

  * :code:`REGEX_SuffixForward`: this is a regular expression that matches all filenames in forward direction. :code:`grepl` function in R is used.

  * :code:`REGEX_SuffixReverse`: this is a regular expression that matches all filenames in reverse direction. :code:`grepl` function in R is used.

If you don't know what regular expression is, don't panic - it's just a way of recognising text. Please refer to :ref:`What is a regular expression?` for more details. Here is an example of how it works in sangeranalseR:

So how sangeranalyseR works is that it first matches the forward and reverse reads by matching :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. Then, sangeranalyseR uses the :code:`str_split` function to split and vectorize their filenames into "contig name" and "direction-suffix" two parts. For those having the same "contig name" will be grouped into the same contig. 


Therefore, it is important to have a consistent naming strategy. You need to make sure that **AB1** files in the same contig group share the same contig name and carefully select your :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. The bad file-naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate wrong results. So, how should we systematically name **AB1** files? We suggest users to follow the file-naming regulation in :ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig_SangerAlignment>`. 

.. _sangeranalyseR_filename_convention_SangerContig_SangerAlignment:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file-naming regulation - *SangerContig*.


As you can see, the first part of the regulation is a consensus read name (or contig name), which helps sangeranalseR to identify which reads should be grouped into the same contig automatically. The second part of the regulation is an index; since there might be more than one read that is in the forward or reverse direction, we recommend you to number your reads in the same contig group. The third part is a direction which is either 'F' (forward) or 'R' (reverse). Last but not least, files have to end with **.ab1** file extension.

To make it more specific, let's go back to the true example. In :ref:`Figure_2<SangerAlignment_file_structure_ab1>`, there are two subdirectories, :code:`ACHLO` and :code:`RBNII`, containing lots of **AB1** files from different contigs in the root directory, :code:`Allolobophora_chlorotica` (:code:`ABIF_Directory`). 


First, we set :code:`REGEX_SuffixForward` to :code:`"_[0-9]*_F.ab1$"` and :code:`REGEX_SuffixReverse` to :code:`"_[0-9]*_R.ab1$"` to let sangeranalyseR match and group forward and reverse reads automatically. By the regular expression rule, :code:`Achl_ACHLO006-09_1_F.ab1`, :code:`Achl_ACHLO007-09_1_F.ab1`, :code:`Achl_ACHLO040-09_1_F.ab1`, :code:`Achl_ACHLO041-09_1_F.ab1`, :code:`Achl_RBNII384-13_1_F.ab1`, :code:`Achl_RBNII395-13_1_F.ab1`, :code:`Achl_RBNII396-13_1_F.ab1`, and :code:`Achl_RBNII397-13_1_F.ab1` are categorized into forward reads, and :code:`Achl_ACHLO006-09_1_R.ab1`, :code:`Achl_ACHLO007-09_1_R.ab1`, :code:`Achl_ACHLO040-09_1_R.ab1`, :code:`Achl_ACHLO041-09_1_R.ab1`, :code:`Achl_RBNII384-13_1_R.ab1`, :code:`Achl_RBNII395-13_1_R.ab1`, :code:`Achl_RBNII396-13_1_R.ab1`, and :code:`Achl_RBNII397-13_1_R.ab1` are categorized into reverse reads. Then, :code:`str_split` function is used to split each filename above into "contig name" and "direction-suffix". Eight contig names are detected in this example which are :code:`Achl_ACHLO006-09`, :code:`Achl_ACHLO007-09`, :code:`Achl_ACHLO040-09`, :code:`Achl_ACHLO041-09`, :code:`Achl_RBNII384-13`, :code:`Achl_RBNII395-13`, :code:`Achl_RBNII396-13`, and :code:`Achl_RBNII397-13`. Last, a loop iterates through all contigs, and sangeranalseR creates each of them into a *SangerContig* instance. You can check :ref:`Advanced User Guide - *SangerContig* (**AB1**)` to see how sangeranalyseR creates a *SangerContig* instance.

The reason why we strongly recommend you to follow this file-naming regulation is that by doing so, you can directly adopt the example regular expression matching values, :code:`"_[0-9]*_F.ab1$"` and :code:`"_[0-9]*_R.ab1$"`, to group reads and reduce chances of error. Everything mentioned above will be done automatically. 



After understanding how parameters work, please refer to :ref:`Creating *SangerAlignment* instance from **AB1**` below to see how sangeranalseR creates *SangerAlignment* instance.


(2) "CSV file matching" *SangerAlignment* inputs (**AB1**)
-----------------------------------------------------------
For those who are not familiar with regular expression, we provide a second grouping approach, CSV file matching method. sangeranalyseR will group **AB1** files based on the information in a CSV file automatically. The note below shows the regulations:

.. note::

    Here is an :download:`example CSV file <../files/SangerAlignment_ab1/names_conversion.csv>` (:ref:`Figure 4<sangeranalyseR_csv_file_SangerAlignment_ab1>`)

      .. _sangeranalyseR_csv_file_SangerAlignment_ab1:
      .. figure::  ../image/sangeranalyseR_csv_file_sangeralignment_ab1.png
         :align:   center
         :scale:   80 %

         Figure 4. Example CSV file for *SangerAlignment* instance creation.  

    *  There must be three columns, "**reads**", "**direction**", and "**contig**", in the CSV file.
    *  The "**reads**" column stores the filename of **AB1** files that are going to be included in the analysis.
    *  The "**direction**" column stores the direction of the reads. It must be "F" (forward) or "R" (reverse).
    *  The "**contig**" column stores the contig name that each read blongs. Reads in the same contig have to have the same contig name, and they will be grouped into the same contig.


There are two parameters, :code:`ABIF_Directory` and :code:`CSV_NamesConversion`,that define the grouping rule to help sangeranalseR to automatically match correct **AB1** files and divide them into forward and reverse directions.

.. note::

  * :code:`ABIF_Directory`: this is the directory that contains all **AB1** files, and it can be either an absolute or relative path. We suggest users to put only target AB1 files inside this directory and do not include any other unrelated files.
  * :code:`CSV_NamesConversion`: this is the path to the CSV file. It can be either an absolute or relative path.

The main difference between "CSV file matching" and "regular expression matching" is where the grouping rule is written. For "regular expression matching", rules are writtein in filenames, and thus more naming requirements are required. In contrast, rules of "CSV file matching" are written in an additional CSV file so it is more flexible on **AB1** file-naming.

So how sangeranalyseR works is that it first reads in the CSV file (with *"reads"*, *"direction"*, and *"contig"* columns), find the names of **AB1** files listed in *"reads"*, group them based on *"contig"*, and assign directions to them based on *"direction"*.

To make it more specific, let's go back to the true example. First, we prepare a :download:`CSV file <../files/SangerAlignment_ab1/names_conversion.csv>` (:code:`CSV_NamesConversion`) and a file directory like :ref:`Figure_2<SangerAlignment_file_structure_ab1>` (:code:`ABIF_Directory`) with **AB1** files from different contigs. In the CSV file, there are 16 rows and 8 distinct contig names. sangeranalyseR matches *"reads"* of these 16 rows to filenames in :code:`Allolobophora_chlorotica` directory. Then sangeranalyseR groups all matched reads, :code:`Achl_ACHLO006-09_1_F.ab1`, :code:`Achl_ACHLO007-09_1_F.ab1`, :code:`Achl_ACHLO040-09_1_F.ab1`, :code:`Achl_ACHLO041-09_1_F.ab1`, :code:`Achl_RBNII384-13_1_F.ab1`, :code:`Achl_RBNII395-13_1_F.ab1`, :code:`Achl_RBNII396-13_1_F.ab1`, :code:`Achl_RBNII397-13_1_F.ab1`,  :code:`Achl_ACHLO006-09_1_R.ab1`, :code:`Achl_ACHLO007-09_1_R.ab1`, :code:`Achl_ACHLO040-09_1_R.ab1`, :code:`Achl_ACHLO041-09_1_R.ab1`, :code:`Achl_RBNII384-13_1_R.ab1`, :code:`Achl_RBNII395-13_1_R.ab1`, :code:`Achl_RBNII396-13_1_R.ab1`, and :code:`Achl_RBNII397-13_1_R.ab1`, into 8 distinct contig names which are :code:`Achl_ACHLO006-09`, :code:`Achl_ACHLO007-09`, :code:`Achl_ACHLO040-09`, :code:`Achl_ACHLO041-09`, :code:`Achl_RBNII384-13`, :code:`Achl_RBNII395-13`, :code:`Achl_RBNII396-13`, and :code:`Achl_RBNII397-13`, by the *"contig"* column. Last, the directions of reads in each contig are assigned by the *"direction"* column. Take :code:`Achl_ACHLO041-09` contig as an example. Its "forward read list" will include :code:`Achl_ACHLO041-09_1_F.ab1`, and its "reverse read list" will include :code:`Achl_ACHLO041-09_1_R.ab1`.


After understanding how parameters work, please refer to :ref:`Creating *SangerAlignment* instance from **AB1**` below to see how sangeranalseR creates *SangerAlignment* instance.

|

Creating *SangerAlignment* instance from **AB1**
+++++++++++++++++++++++++++++++++++++++++++++++++
After preparing the input directory, we can create a *SangerAlignment* instance by running :code:`SangerAlignment` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Their input parameters are same, and all of them have their default values. For more details about *SangerAlignment* inputs and slots definition, please refer to `sangeranalyseR reference manual <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. We will explain two *SangerAlignment* instance creation methods, "regular expression matching" and "CSV file matching".




(1) "regular expression matching" *SangerAlignment* creation (**AB1**)
-----------------------------------------------------------------------

The consturctor function and :code:`new` method below contain three parameters, :code:`ABIF_Directory`, :code:`REGEX_SuffixForward`, and :code:`REGEX_SuffixReverse`, that we mentioned in the previous section. It also includes important parameters like quality trimming, chromatogram visualization, consensus alignment, contigs alignment, and so on. Run the following code and create :code:`my_sangerAlignment` instance. 


.. code-block:: R
                     
   # using `constructor` function to create SangerAlignment instance
   my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                         processMethod        = "REGEX",
                                         ABIF_Directory       = parentDir,
                                         REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                                         REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                                         TrimmingMethod       = "M1",
                                         M1TrimmingCutoff     = 0.0001,
                                         M2CutoffQualityScore = NULL,
                                         M2SlidingWindowSize  = NULL,
                                         baseNumPerRow        = 100,
                                         heightPerRow         = 200,
                                         signalRatioCutoff    = 0.33,
                                         showTrimmed          = TRUE,
                                         refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                         minReadsNum          = 2,
                                         minReadLength        = 20,
                                         minFractionCall      = 0.5,
                                         maxFractionLost      = 0.5,
                                         geneticCode          = GENETIC_CODE,
                                         acceptStopCodons     = TRUE,
                                         readingFrame         = 1,
                                         processorsNum        = 2)


   # using `new` method to create SangerAlignment instance
   my_sangerAlignment <- new("SangerAlignment",
                             inputSource          = "ABIF",
                             processMethod        = "REGEX",
                             ABIF_Directory       = parentDir,
                             REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                             REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                             TrimmingMethod       = "M1",
                             M1TrimmingCutoff     = 0.0001,
                             M2CutoffQualityScore = NULL,
                             M2SlidingWindowSize  = NULL,
                             baseNumPerRow        = 100,
                             heightPerRow         = 200,
                             signalRatioCutoff    = 0.33,
                             showTrimmed          = TRUE,
                             refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                             minReadsNum          = 2,
                             minReadLength        = 20,
                             minFractionCall      = 0.5,
                             maxFractionLost      = 0.5,
                             geneticCode          = GENETIC_CODE,
                             acceptStopCodons     = TRUE,
                             readingFrame         = 1,
                             processorsNum        = 2)



In this example, 16 reads are detected and 8 distinct *SangerContig* instances are created. These *SangerContig* instances are stored in a "contig list" in :code:`my_sangerAlignment`, which will be used as the input for the following functions.

Inside the R shell, you can run :code:`my_sangerAlignment` to get basic information of the instance or run :code:`my_sangerAlignment@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerAlignment` is successfully created.

Here is the output of :code:`my_sangerAlignment`::

   SangerAlignment S4 instance
            Input Source :  ABIF 
            Process Method :  REGEX 
            ABIF Directory :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica 
      REGEX Suffix Forward :  _[0-9]*_F.ab1$ 
      REGEX Suffix Reverse :  _[0-9]*_R.ab1$ 
         Contigs Consensus :  TTATAYTTTATTYTRGGCGTCTGAAGCAGGATAGTAGGAGCYGGTATAAGACTCCTAATTCGAATTGAGCTAAGACARCCGGGAGCATTCCTAGGAAGRGATCAACTCTATAACACTATTGTAACTGCTCACGCATTTGTAATAATTTTCTTTCTAGTAATACCTGTATTTATTGGGGGGTTCGGTAATTGACTTCTACCTTTAATACTTGGAGCCCCTGACATGGCATTCCCACGACTTAACAACATAAGATTCTGACTCCTTCCCCCATCACTAATCCTTCTAGTGTCCTCTGCTGCAGTAGAAAAAGGTGCBGGAACTGGATGAACTGTTTATCCRCCCCTAGCAAGAAATATTGCTCATGCCGGCCCATCTGTAGACTTAGCTATYTTTTCTCTTCATTTAGCAGGTGCTTCATCAATCTTAGGKGCYATTAATTTTATYACTACTGTTATTAACATACGATGAAGAGGCTTACGACTTGAACGAATCCCATTATTCGTTTGAGCCGTACTAATTACAGTGGTHCTTCTACTCCTATCYTTACCAGTATTAGCCGGTGCRATTACYATACTACTTACCGATCGAAATCTAAATACCTCCTTCTTTGAYCCTGCTGGAGGTGGAGATCCCATCCTCTACCAACACTTATTCTGATTTTTTGGTCACCCTGAG 
   SUCCESS [2021-13-07 23:16:16] 'SangerAlignment' is successfully created!

Here is the output of :code:`my_sangerAlignment@objectResults@readResultTable`::

                     readName creationResult errorType errorMessage inputSource    direction
   1  Achl_ACHLO006-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   2  Achl_ACHLO006-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   3  Achl_ACHLO007-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   4  Achl_ACHLO007-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   5  Achl_ACHLO040-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   6  Achl_ACHLO040-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   7  Achl_ACHLO041-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   8  Achl_ACHLO041-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   9  Achl_RBNII384-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   10 Achl_RBNII384-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   11 Achl_RBNII395-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   12 Achl_RBNII395-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   13 Achl_RBNII396-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   14 Achl_RBNII396-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   15 Achl_RBNII397-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   16 Achl_RBNII397-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read



(2) "CSV file matching" *SangerAlignment* creation (**AB1**)
-------------------------------------------------------------
The consturctor function and :code:`new` method below contain two parameters, :code:`ABIF_Directory`, and :code:`CSV_NamesConversion`, that we mentioned in the previous section. It also includes important parameters like quality trimming, chromatogram visualization, consensus alignment, contigs alignment, and so on. Run the following code and create :code:`my_sangerAlignment` instance. 

.. code-block:: R

   csv_namesConversion <- file.path(rawDataDir, "ab1", "SangerAlignment", "names_conversion_all.csv")

   # using `constructor` function to create SangerAlignment instance
   my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                         processMethod        = "CSV",
                                         ABIF_Directory       = parentDir,
                                         CSV_NamesConversion  = csv_namesConversion,
                                         TrimmingMethod       = "M1",
                                         M1TrimmingCutoff     = 0.0001,
                                         M2CutoffQualityScore = NULL,
                                         M2SlidingWindowSize  = NULL,
                                         baseNumPerRow        = 100,
                                         heightPerRow         = 200,
                                         signalRatioCutoff    = 0.33,
                                         showTrimmed          = TRUE,
                                         refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                         minReadsNum          = 2,
                                         minReadLength        = 20,
                                         minFractionCall      = 0.5,
                                         maxFractionLost      = 0.5,
                                         geneticCode          = GENETIC_CODE,
                                         acceptStopCodons     = TRUE,
                                         readingFrame         = 1,
                                         processorsNum        = 1)


   # using `new` method to create SangerAlignment instance
   my_sangerAlignment <- new("SangerAlignment",
                             processMethod        = "CSV",
                             ABIF_Directory       = parentDir,
                             CSV_NamesConversion  = csv_namesConversion,
                             TrimmingMethod       = "M1",
                             M1TrimmingCutoff     = 0.0001,
                             M2CutoffQualityScore = NULL,
                             M2SlidingWindowSize  = NULL,
                             baseNumPerRow        = 100,
                             heightPerRow         = 200,
                             signalRatioCutoff    = 0.33,
                             showTrimmed          = TRUE,
                             refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                             minReadsNum          = 2,
                             minReadLength        = 20,
                             minFractionCall      = 0.5,
                             maxFractionLost      = 0.5,
                             geneticCode          = GENETIC_CODE,
                             acceptStopCodons     = TRUE,
                             readingFrame         = 1,
                             processorsNum        = 1)


First, you need to load the CSV file into the R environment. If you are still don't know how to prepare it, please check :ref:`(2) "CSV file matching" *SangerAlignment* inputs (**AB1**)`. Then, it will follow rules in the CSV file and create :code:`my_sangerAlignment`. After it's created, inside the R shell, you can run :code:`my_sangerAlignment` to get basic information of the instance or run :code:`my_sangerAlignment@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerAlignment` is successfully created.

Here is the output of :code:`my_sangerAlignment`::

   SangerAlignment S4 instance
            Input Source :  ABIF 
            Process Method :  CSV 
            ABIF Directory :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/Allolobophora_chlorotica 
      CSV Names Conversion :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/ab1/SangerAlignment/names_conversion_all.csv 
         Contigs Consensus :  TTATAYTTTATTYTRGGCGTCTGAAGCAGGATAGTAGGAGCYGGTATAAGACTCCTAATTCGAATTGAGCTAAGACARCCGGGAGCATTCCTAGGAAGRGATCAACTCTATAACACTATTGTAACTGCTCACGCATTTGTAATAATTTTCTTTCTAGTAATACCTGTATTTATTGGGGGGTTCGGTAATTGACTTCTACCTTTAATACTTGGAGCCCCTGACATGGCATTCCCACGACTTAACAACATAAGATTCTGACTCCTTCCCCCATCACTAATCCTTCTAGTGTCCTCTGCTGCAGTAGAAAAAGGTGCBGGAACTGGATGAACTGTTTATCCRCCCCTAGCAAGAAATATTGCTCATGCCGGCCCATCTGTAGACTTAGCTATYTTTTCTCTTCATTTAGCAGGTGCTTCATCAATCTTAGGKGCYATTAATTTTATYACTACTGTTATTAACATACGATGAAGAGGCTTACGACTTGAACGAATCCCATTATTCGTTTGAGCCGTACTAATTACAGTGGTHCTTCTACTCCTATCYTTACCAGTATTAGCCGGTGCRATTACYATACTACTTACCGATCGAAATCTAAATACCTCCTTCTTTGAYCCTGCTGGAGGTGGAGATCCCATCCTCTACCAACACTTATTCTGATTTTTTGGTCACCCTGAG 
   SUCCESS [2021-14-07 01:48:28] 'SangerAlignment' is successfully created!


Here is the output of :code:`my_sangerAlignment@objectResults@readResultTable`::

                     readName creationResult errorType errorMessage inputSource    direction
   1  Achl_ACHLO006-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   2  Achl_ACHLO006-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   3  Achl_ACHLO007-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   4  Achl_ACHLO007-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   5  Achl_ACHLO040-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   6  Achl_ACHLO040-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   7  Achl_ACHLO041-09_1_F.ab1           TRUE      None         None        ABIF Forward Read
   8  Achl_ACHLO041-09_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   9  Achl_RBNII384-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   10 Achl_RBNII384-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   11 Achl_RBNII395-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   12 Achl_RBNII395-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   13 Achl_RBNII396-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   14 Achl_RBNII396-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read
   15 Achl_RBNII397-13_1_F.ab1           TRUE      None         None        ABIF Forward Read
   16 Achl_RBNII397-13_2_R.ab1           TRUE      None         None        ABIF Reverse Read

|

Updating *SangerAlignment* quality trimming parameters
++++++++++++++++++++++++++++++++++++++++++++++++++++++++
In the previous :ref:`Creating *SangerAlignment* instance from **AB1**` part, the constructor function will apply the quality trimming parameters to all reads. After creating a *SangerAlignment* S4 instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will update all reads with the new trimming parameters and redo reads alignment in *SangerContig* and contigs alignment in *SangerAlignment*. If users want to do quality trimming read by read instead all at once, please read :ref:`Launching *SangerAlignment* Shiny app`.

.. code-block:: R

   newSangerAlignment <- updateQualityParam(my_sangerAlignment,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 29,
                                            M2SlidingWindowSize  = 15)

|

Launching *SangerAlignment* Shiny app
++++++++++++++++++++++++++++++++++++++
We create an interactive local Shiny app for users to go into each *SangerRead* and *SangerContig* in *SangerAlignment* instance. Users only need to run one function with previously created instance as input, :code:`my_sangerAlignment`, and the *SangerAlignment* Shiny app will pop up. Here, we will go through pages in the three levels.

.. code-block:: R

   launchApp(my_sangerAlignment)


*SangerAlignment* page (SA app)
-------------------------------
:ref:`Figure 5<SangerAlignment_ShinyApp_1>` is the initial page and the toppest layer of *SangerAlignment* App. It provides basic parameters in *SangerAlignment* instance, contigs alignment result and phylogenetic tree etc. Before checking the results, users need to click “Re-calculate Contigs Alignment” button to do contigs alignment in order to get the updated results. From the left-hand side panel, we can clearly see the hierarchy of the *SangerAlignment* S4 instance and easily access to all reads and contigs in it.

.. _SangerAlignment_ShinyApp_1:
.. figure::  ../image/SangerAlignment_ShinyApp_1.png
   :align:   center
   :scale:   40 %

   Figure 5. *SangerAlignment* Shiny app initial page - *SangerAlignment* Page.

Scroll down a bit, users can see the contigs alignment result generated by `DECIPHER <https://bioconductor.org/packages/release/bioc/html/DECIPHER.html>`_ R package embedded in *SangerAlignment* page. :ref:`Figure 6<SangerAlignment_ShinyApp_2>` shows the contigs alignment result.

.. _SangerAlignment_ShinyApp_2:
.. figure::  ../image/SangerAlignment_ShinyApp_2.png
   :align:   center
   :scale:   50 %

   Figure 6. *SangerAlignment* Page - contigs alignment result.

In *SangerAlignment* page, the phylogenetic tree result is provided as well (:ref:`Figure 7<SangerAlignment_ShinyApp_3>`). The tree is generated by `ape <https://cran.r-project.org/web/packages/ape/index.html>`_ R package which uses neighbor-joining algorithm.

.. _SangerAlignment_ShinyApp_3:
.. figure::  ../image/SangerAlignment_ShinyApp_3.png
   :align:   center
   :scale:   50 %

   Figure 7. *SangerAlignment* Page - phylogenetic tree result.


*SangerContig* page (SA app)
-----------------------------
Now, let's go to the page in the next level, *SangerContig* page. Users can click into all contigs and check their results. :ref:`Figure 8<SangerAlignment_ShinyApp_5>` shows the overview page of Contig 1. Notice that there is a red “Re-calculate Contig” button. After changing the quality trimming parameters, users need to click the button before checking the results below in order to get the updated information.

.. _SangerAlignment_ShinyApp_5:
.. figure::  ../image/SangerAlignment_ShinyApp_5.png
   :align:   center
   :scale:   37 %

   Figure 8. *SangerAlignment* Shiny app - *SangerContig* page.

The information provided in this page includes : “input parameters”, “genetic code table”, “reference amino acid sequence”, “reads alignment”, “difference data frame”, “dendrogram”, “sample distance heatmap”, “indels data frame”, “stop codons data frame”. :ref:`Figure 9<SangerAlignment_ShinyApp_6>` and :ref:`Figure 10<SangerAlignment_ShinyApp_7>` show part of the results in the *SangerContig* page. The results are dynamic based on the trimming parameters from user inputs.

.. _SangerAlignment_ShinyApp_6:
.. figure::  ../image/SangerAlignment_ShinyApp_6.png
   :align:   center
   :scale:   50 %

   Figure 9. *SangerContig* page - contig-related parameters, genetic code and reference amino acid sequence.

.. _SangerAlignment_ShinyApp_7:
.. figure::  ../image/SangerAlignment_ShinyApp_7.png
   :align:   center
   :scale:   50 %

   Figure 10. *SangerContig* page - reads alignment and difference data frame.


*SangerRead* page (SA app)
---------------------------
Now, let's go to the page in the lowest level, *SangerRead* page. *SangerRead* page contains all details of a read including its trimming and chromatogram inputs and results. All reads are in "forward" or "reverse" direction. Under "Contig Overview" tab (*SangerContig* page), there are two expendable tabs, “Forward Reads” and “Reverse Reads” storing corresponding reads on the left-hand side navigation panel in :ref:`Figure 11<SangerAlignment_ShinyApp_8>`. In this example, there are one read in each tab and :ref:`Figure 11<SangerAlignment_ShinyApp_8>` shows the “1 - 1 Forward Read” page. It provides basic information, quality trimming inputs, chromatogram plotting inputs etc. Primary/secondary sequences in this figure are dynamic based on the :code:`signalRatioCutoff` value for base calling and the length of them are always same. Another thing to mention is that primary/secondary sequences and the sequences in the chromatogram in :ref:`Figure 16<SangerAlignment_ShinyApp_14>` below will always be same after trimming and their color codings for A/T/C/G are same as well.

.. _SangerAlignment_ShinyApp_8:
.. figure::  ../image/SangerAlignment_ShinyApp_8.png
   :align:   center
   :scale:   35 %

   Figure 11. *SangerAlignment* Shiny app - *SangerRead* page.

In quality trimming steps, we removes fragment at both ends of sequencing reads with low quality score. It is important because trimmed reads will improves alignment results. :ref:`Figure 12<SangerAlignment_ShinyApp_9>` shows the UI for Trimming Method 1 (M1): ‘Modified Mott Trimming’. This method is implemented in `Phred <http://www.phrap.org/phredphrapconsed.html>`_. Users can change the cutoff score and click “Apply Trimming Parameters" button to update the UI. The value of input must be between 0 and 1. If the input is invalid, the cutoff score will be set to default 0.0001.

.. _SangerAlignment_ShinyApp_9:
.. figure::  ../image/SangerAlignment_ShinyApp_9.png
   :align:   center
   :scale:   30 %

   Figure 12. *SangerRead* page - Trimming Method 1 (M1): ‘Modified Mott Trimming’ UI.

:ref:`Figure 13<SangerAlignment_ShinyApp_10>` shows another quality trimming methods for users to choose from, Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’. This method is implemented in `Trimmomatics <http://www.usadellab.org/cms/?page=trimmomatic>`_. Users can change the cutoff quality score as well as sliding window size and click “Apply Trimming Parameters" button to update the UI. The value of cutoff quality score must be between 0 and 60 (default 20); the value of sliding window size must be between 0 and 40 (default 10). If the inputs are invalid, their values will be set to default.

.. _SangerAlignment_ShinyApp_10:
.. figure::  ../image/SangerAlignment_ShinyApp_10.png
   :align:   center
   :scale:   30 %

   Figure 13. *SangerRead* page - Trimming Method 2 (M2): ‘Trimmomatics Sliding Window Trimming’ UI.

:ref:`Figure 14<SangerAlignment_ShinyApp_11>` shows the quality report before and after trimming. After clicking the “Apply Trimming Parameters” button, the values of these information boxes will be updated to the latest values.

.. _SangerAlignment_ShinyApp_11:
.. figure::  ../image/SangerAlignment_ShinyApp_11.png
   :align:   center
   :scale:   45 %

   Figure 14. *SangerRead* page - read quality report before / after trimming.

In :ref:`Figure 15<SangerAlignment_ShinyApp_13>`, the x-axis is the index of the base pairs; the y-axis is the Phred quality score. The green horizontal bar at the top of the plot is the raw read region and the orange horizontal bar represents the trimmed read region. Both :ref:`Figure 15<SangerAlignment_ShinyApp_13>` trimming plot and :ref:`Figure 16<SangerAlignment_ShinyApp_14>` chromatogram will be updated once users change the quality trimming parameters and click the “Apply Trimming Parameters" button in :ref:`Figure 16<SangerAlignment_ShinyApp_14>`.

.. _SangerAlignment_ShinyApp_13:
.. figure::  ../image/SangerAlignment_ShinyApp_13.png
   :align:   center
   :scale:   50 %

   Figure 15. *SangerRead* page - quality trimming plot.

If we only see primary and secondary sequences in the table, we will loose some variations. Chromatogram is very helpful to check the peak resolution. :ref:`Figure 16<SangerAlignment_ShinyApp_14>` shows the panel of plotting chromatogram. Users can change four parameters: :code:`Base Number Per Row`, :code:`Height Per Row`, :code:`Signal Ratio Cutoff`, and :code:`Show Trimmed Region`. Among them, :code:`Signal Ratio Cutoff` is the key parameter. If its value is default value 0.33, it indicates that the lower peak should be at least 1/3rd as high as the higher peak for it count as a secondary peak.

.. _SangerAlignment_ShinyApp_14:
.. figure::  ../image/SangerAlignment_ShinyApp_14.png
   :align:   center
   :scale:   45 %

   Figure 16. *SangerRead* page - chromatogram panel.

Here is an example of applying new chromatogram parameters. We click “Show Trimmed Region” to set its value from FALSE to TRUE. :ref:`Figure 17<SangerAlignment_ShinyApp_15>` shows the loading notification popup during base calling and chromatogram plotting.

.. _SangerAlignment_ShinyApp_15:
.. figure::  ../image/SangerAlignment_ShinyApp_15.png
   :align:   center
   :scale:   45 %

   Figure 17. *SangerRead* page - loading notification popup during replotting chromatogram.

After replotting the chromatogram, trimmed region is showed in red striped region. :ref:`Figure 18<SangerAlignment_ShinyApp_16>` shows part of the the chromatogram (1 bp ~ 240 bp). Moreover, chromatogram will be replotted when trimmed positions or chromatogram parameters are updated.

.. _SangerAlignment_ShinyApp_16:
.. figure::  ../image/SangerAlignment_ShinyApp_16.png
   :align:   center
   :scale:   50 %

   Figure 18. *SangerRead* page - chromatogram with trimmed region showed.

To let users browse the trimmed primary/secondary sequences without finding “Trimming Start Point” and “Trimming End Point” by themselves, we provide the final trimmed primary/secondary sequences that will be used for reads alignment in table format with quality scores in :ref:`Figure 19<SangerAlignment_ShinyApp_17>`. Frameshift amino acid sequences are also provided.

.. _SangerAlignment_ShinyApp_17:
.. figure::  ../image/SangerAlignment_ShinyApp_17.png
   :align:   center
   :scale:   45 %

   Figure 19. *SangerRead* page - trimmed primary/secondary sequences and Phred quality score in table format.

We have updated the trimming and chromatogram parameters for each read. Now, we need to click “Re-calculate contig” button to do alignment again. Last but not least, we can save all data into a new ‘SangerContig’ S4 instance by clicking “Save S4 instance button”. New S4 instance will be saved in **Rda** format. Users can run :code:`readRDS` function to load it into current R environment. :ref:`Figure 20<SangerAlignment_ShinyApp_18>` shows some hints in the save notification popup.

.. _SangerAlignment_ShinyApp_18:
.. figure::  ../image/SangerAlignment_ShinyApp_18.png
   :align:   center
   :scale:   40 %

   Figure 20. *SangerRead* page - saving notification popup.


|


Writing *SangerAlignment* FASTA files :sub:`(AB1)`
+++++++++++++++++++++++++++++++++++++++++++++++++++
Users can write the *SangerAlignment* instance, :code:`my_sangerAlignment`, to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`contigs_unalignment`: Writing contigs into a single **FASTA** file.
* :code:`contigs_alignment`: Writing contigs alignment and contigs consensus read to a single **FASTA** file.
* :code:`all_reads`: Writing all reads to a single **FASTA** file.
* :code:`all`: Writing contigs, contigs alignment, and all reads into three different files.

Below is the oneliner for writing out **FASTA** files. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFasta` function.

.. code-block:: R

   writeFasta(my_sangerAlignment,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA,
              selection         = "all")

Users can download the output FASTA file of this example through the following three links:

(1) :download:`Sanger_contigs_unalignment.fa <../files/SangerAlignment_ab1/Sanger_contigs_unalignment.fa>`
(2) :download:`Sanger_contigs_alignment.fa <../files/SangerAlignment_ab1/Sanger_contigs_alignment.fa>`
(3) :download:`Sanger_all_trimmed_reads.fa <../files/SangerAlignment_ab1/Sanger_all_trimmed_reads.fa>`

|

Generating *SangerAlignment* report :sub:`(AB1)`
+++++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save *SangerAlignment* instance, :code:`my_sangerAlignment`, into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.

Users can set :code:`includeSangerContig` and :code:`includeSangerRead` parameters to decide to which level the *SangerAlignment* report will go. Moreover, after the reports are generated, users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it will take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` and :code:`includeSangerContig` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReport(my_sangerAlignment,
                  outputDir           = tempdir(),
                  includeSangerRead   = FALSE,
                  includeSangerContig = FALSE)


Here is the generated `SangerAlignment html report of this example (ABIF) <https://kuanhao-chao.github.io/sangeranalyseR_report/SangerAlignment/AB1/SangerAlignment/SangerAlignment_Report.html>`_. Users can access to '*Basic Information*', '*Contigs Consensus*', '*Contigs Alignment*', '*Contigs Tree*', and '*Contig Reports*' sections inside it. Furthermore, users can also navigate through html reports of all contigs and forward and reverse *SangerRead* in this *SangerAlignment* report.


-----

|
|


Code summary (*SangerAlignment*, **AB1**)
++++++++++++++++++++++++++++++++++++++++++++++++++


(1) Preparing *SangerAlignment* **AB1** inputs
-----------------------------------------------

.. code-block:: R

      rawDataDir <- system.file("extdata", package = "sangeranalyseR")
      parentDir <- file.path(rawDataDir, 'Allolobophora_chlorotica')

|

(2) Creating *SangerAlignment* instance from **AB1**
------------------------------------------------------

(2.1) "Regular Expression Method" *SangerAlignment* creation (**AB1**)
***********************************************************************

.. code-block:: R

   # using `constructor` function to create SangerAlignment instance
   my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                         processMethod        = "REGEX",
                                         ABIF_Directory       = parentDir,
                                         REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                                         REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                                         refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


   # using `new` method to create SangerAlignment instance
   my_sangerAlignment <- new("SangerAlignment",
                             inputSource          = "ABIF",
                             processMethod        = "REGEX",
                             ABIF_Directory       = parentDir,
                             REGEX_SuffixForward  = "_[0-9]*_F.ab1$",
                             REGEX_SuffixReverse  = "_[0-9]*_R.ab1$",
                             refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")

.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

      INFO [2021-14-07 03:10:24] #################################################
      INFO [2021-14-07 03:10:24] #### Start creating SangerAlignment instance ####
      INFO [2021-14-07 03:10:24] #################################################
      INFO [2021-14-07 03:10:24]   >> You are using Regular Expression Method to group AB1 files!
      INFO [2021-14-07 03:10:24] ========================================================
      INFO [2021-14-07 03:10:24] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:24] ========================================================
      INFO [2021-14-07 03:10:24]   >> Contig Name: 'Achl_ACHLO006-09'
      SUCCESS [2021-14-07 03:10:24] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:24] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:24] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:24]    >> 'Achl_ACHLO006-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:24] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:24] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:24] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:24]    >> 'Achl_ACHLO006-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:24]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:10:26] ==========================================================
      SUCCESS [2021-14-07 03:10:26] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:26] ==========================================================
      INFO [2021-14-07 03:10:26]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:26]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:26]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:26] ========================================================
      INFO [2021-14-07 03:10:26] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:26] ========================================================
      INFO [2021-14-07 03:10:26]   >> Contig Name: 'Achl_ACHLO007-09'
      SUCCESS [2021-14-07 03:10:27] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:27] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:27] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:27]    >> 'Achl_ACHLO007-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:27] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:27] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:27] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:27]    >> 'Achl_ACHLO007-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:27]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:10:29] ==========================================================
      SUCCESS [2021-14-07 03:10:29] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:29] ==========================================================
      INFO [2021-14-07 03:10:29]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:29]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:29]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:29] ========================================================
      INFO [2021-14-07 03:10:29] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:29] ========================================================
      INFO [2021-14-07 03:10:29]   >> Contig Name: 'Achl_ACHLO040-09'
      SUCCESS [2021-14-07 03:10:29] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:29] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:29] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:29]    >> 'Achl_ACHLO040-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:30] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:30] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:30] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:30]    >> 'Achl_ACHLO040-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:30]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:10:31] ==========================================================
      SUCCESS [2021-14-07 03:10:31] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:31] ==========================================================
      INFO [2021-14-07 03:10:31]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:31]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:31]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:31] ========================================================
      INFO [2021-14-07 03:10:31] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:31] ========================================================
      INFO [2021-14-07 03:10:31]   >> Contig Name: 'Achl_ACHLO041-09'
      SUCCESS [2021-14-07 03:10:31] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:31] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:31] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:32]    >> 'Achl_ACHLO041-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:32] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:32] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:32] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:32]    >> 'Achl_ACHLO041-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:32]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:10:33] ==========================================================
      SUCCESS [2021-14-07 03:10:33] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:33] ==========================================================
      INFO [2021-14-07 03:10:33]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:33]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:33]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:33] ========================================================
      INFO [2021-14-07 03:10:33] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:33] ========================================================
      INFO [2021-14-07 03:10:33]   >> Contig Name: 'Achl_RBNII384-13'
      SUCCESS [2021-14-07 03:10:34] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:34] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:34] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:34]    >> 'Achl_RBNII384-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:34] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:34] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:34] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:34]    >> 'Achl_RBNII384-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:34]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:10:35] ==========================================================
      SUCCESS [2021-14-07 03:10:35] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:35] ==========================================================
      INFO [2021-14-07 03:10:35]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:35]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:35]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:35] ========================================================
      INFO [2021-14-07 03:10:35] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:35] ========================================================
      INFO [2021-14-07 03:10:35]   >> Contig Name: 'Achl_RBNII395-13'
      SUCCESS [2021-14-07 03:10:36] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:36] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:36] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:36]    >> 'Achl_RBNII395-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:36] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:36] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:36] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:36]    >> 'Achl_RBNII395-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:36]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:10:37] ==========================================================
      SUCCESS [2021-14-07 03:10:37] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:37] ==========================================================
      INFO [2021-14-07 03:10:37]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:37]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:37]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:37] ========================================================
      INFO [2021-14-07 03:10:37] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:37] ========================================================
      INFO [2021-14-07 03:10:37]   >> Contig Name: 'Achl_RBNII396-13'
      SUCCESS [2021-14-07 03:10:38] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:38] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:38] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:38]    >> 'Achl_RBNII396-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:38] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:38] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:38] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:38]    >> 'Achl_RBNII396-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:38]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:10:40] ==========================================================
      SUCCESS [2021-14-07 03:10:40] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:40] ==========================================================
      INFO [2021-14-07 03:10:40]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:40]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:40]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:40] ========================================================
      INFO [2021-14-07 03:10:40] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:10:40] ========================================================
      INFO [2021-14-07 03:10:40]   >> Contig Name: 'Achl_RBNII397-13'
      SUCCESS [2021-14-07 03:10:40] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:40] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:40] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:40]    >> 'Achl_RBNII397-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:10:40] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:40] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:10:40] --------------------------------------------------------
      SUCCESS [2021-14-07 03:10:40]    >> 'Achl_RBNII397-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:10:40]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:10:42] ==========================================================
      SUCCESS [2021-14-07 03:10:42] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:10:42] ==========================================================
      INFO [2021-14-07 03:10:42]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:10:42]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:42]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      INFO [2021-14-07 03:10:42] Aligning consensus reads ... 
      SUCCESS [2021-14-07 03:10:43] #############################################################
      SUCCESS [2021-14-07 03:10:43] ######## 'SangerAlignment' S4 instance is created !! ########
      SUCCESS [2021-14-07 03:10:43] #############################################################
      INFO [2021-14-07 03:10:43]   >> 8 contigs detected from 'regular expression'.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_ACHLO006-09':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_ACHLO007-09':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_ACHLO040-09':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_ACHLO041-09':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_RBNII384-13':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_RBNII395-13':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_RBNII396-13':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]       >> Contig 'Achl_RBNII397-13':
      INFO [2021-14-07 03:10:43]           >> 1 forward reads.
      INFO [2021-14-07 03:10:43]           >> 1 reverse reads.
      INFO [2021-14-07 03:10:43]   >> 16 reads created from ABIF file.
      INFO [2021-14-07 03:10:43]   >> Reads are trimmed by 'M1 - Mott’s trimming algorithm'.
      DEBUG [2021-14-07 03:10:43]    >> For more information, please run 'object'.
      DEBUG [2021-14-07 03:10:43]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads

|

(2.2) "CSV file matching" *SangerAlignment* creation (**AB1**)
*****************************************************************


.. code-block:: R

   csv_namesConversion <- file.path(rawDataDir, "ab1", "SangerAlignment", "names_conversion_all.csv")

   # using `constructor` function to create SangerAlignment instance
   my_sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                         processMethod        = "CSV",
                                         ABIF_Directory       = parentDir,
                                         CSV_NamesConversion  = csv_namesConversion,
                                         refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


   # using `new` method to create SangerAlignment instance
   my_sangerAlignment <- new("SangerAlignment",
                             processMethod        = "CSV",
                             ABIF_Directory       = parentDir,
                             CSV_NamesConversion  = csv_namesConversion,
                             refAminoAcidSeq      = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

      WARN [2021-14-07 03:11:43] 'Achl_ACHLO006-09_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO006-09_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO007-09_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO007-09_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO040-09_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO040-09_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO041-09_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_ACHLO041-09_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII384-13_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII384-13_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII395-13_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII395-13_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII396-13_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII396-13_2_R.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII397-13_1_F.ab1' is not in the parent directory.
      WARN [2021-14-07 03:11:43] 'Achl_RBNII397-13_2_R.ab1' is not in the parent directory.
      INFO [2021-14-07 03:11:43] #################################################
      INFO [2021-14-07 03:11:43] #### Start creating SangerAlignment instance ####
      INFO [2021-14-07 03:11:43] #################################################
      INFO [2021-14-07 03:11:43] **** You are using CSV Name Conversion Method to group AB1 files!
      INFO [2021-14-07 03:11:43] **** Contig number in your Csv file is 8
      INFO [2021-14-07 03:11:43] ========================================================
      INFO [2021-14-07 03:11:43] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:43] ========================================================
      INFO [2021-14-07 03:11:43]   >> Contig Name: 'Achl_ACHLO006-09'
      SUCCESS [2021-14-07 03:11:43] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:43] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:43] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:43]    >> 'Achl_ACHLO006-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:44] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:44] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:44] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:44]    >> 'Achl_ACHLO006-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:44]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:11:45] ==========================================================
      SUCCESS [2021-14-07 03:11:45] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:45] ==========================================================
      INFO [2021-14-07 03:11:45]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:45]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:45]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:45] ========================================================
      INFO [2021-14-07 03:11:45] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:45] ========================================================
      INFO [2021-14-07 03:11:45]   >> Contig Name: 'Achl_ACHLO007-09'
      SUCCESS [2021-14-07 03:11:45] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:45] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:45] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:45]    >> 'Achl_ACHLO007-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:46] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:46] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:46] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:46]    >> 'Achl_ACHLO007-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:46]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:11:47] ==========================================================
      SUCCESS [2021-14-07 03:11:47] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:47] ==========================================================
      INFO [2021-14-07 03:11:47]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:47]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:47]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:47] ========================================================
      INFO [2021-14-07 03:11:47] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:47] ========================================================
      INFO [2021-14-07 03:11:47]   >> Contig Name: 'Achl_ACHLO040-09'
      SUCCESS [2021-14-07 03:11:47] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:47] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:47] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:47]    >> 'Achl_ACHLO040-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:48] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:48] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:48] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:48]    >> 'Achl_ACHLO040-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:48]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.04 secs
      SUCCESS [2021-14-07 03:11:49] ==========================================================
      SUCCESS [2021-14-07 03:11:49] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:49] ==========================================================
      INFO [2021-14-07 03:11:49]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:49]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:49]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:49] ========================================================
      INFO [2021-14-07 03:11:50] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:50] ========================================================
      INFO [2021-14-07 03:11:50]   >> Contig Name: 'Achl_ACHLO041-09'
      SUCCESS [2021-14-07 03:11:50] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:50] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:50] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:50]    >> 'Achl_ACHLO041-09_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:51] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:51] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:51] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:51]    >> 'Achl_ACHLO041-09_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:51]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.04 secs
      SUCCESS [2021-14-07 03:11:52] ==========================================================
      SUCCESS [2021-14-07 03:11:52] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:52] ==========================================================
      INFO [2021-14-07 03:11:52]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:52]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:52]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:52] ========================================================
      INFO [2021-14-07 03:11:52] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:52] ========================================================
      INFO [2021-14-07 03:11:52]   >> Contig Name: 'Achl_RBNII384-13'
      SUCCESS [2021-14-07 03:11:53] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:53] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:53] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:53]    >> 'Achl_RBNII384-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:53] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:53] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:53] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:53]    >> 'Achl_RBNII384-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:53]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.15 secs
      SUCCESS [2021-14-07 03:11:56] ==========================================================
      SUCCESS [2021-14-07 03:11:56] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:56] ==========================================================
      INFO [2021-14-07 03:11:56]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:56]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:56]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:56] ========================================================
      INFO [2021-14-07 03:11:56] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:56] ========================================================
      INFO [2021-14-07 03:11:56]   >> Contig Name: 'Achl_RBNII395-13'
      SUCCESS [2021-14-07 03:11:56] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:56] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:56] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:56]    >> 'Achl_RBNII395-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:11:57] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:57] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:57] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:57]    >> 'Achl_RBNII395-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:11:57]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.08 secs
      SUCCESS [2021-14-07 03:11:59] ==========================================================
      SUCCESS [2021-14-07 03:11:59] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:11:59] ==========================================================
      INFO [2021-14-07 03:11:59]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:11:59]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:11:59]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:11:59] ========================================================
      INFO [2021-14-07 03:11:59] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:11:59] ========================================================
      INFO [2021-14-07 03:11:59]   >> Contig Name: 'Achl_RBNII396-13'
      SUCCESS [2021-14-07 03:11:59] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:59] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:11:59] --------------------------------------------------------
      SUCCESS [2021-14-07 03:11:59]    >> 'Achl_RBNII396-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:12:00] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:00] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:12:00] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:00]    >> 'Achl_RBNII396-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:12:00]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.03 secs
      SUCCESS [2021-14-07 03:12:01] ==========================================================
      SUCCESS [2021-14-07 03:12:01] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:12:01] ==========================================================
      INFO [2021-14-07 03:12:01]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:12:01]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:12:01]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:12:01] ========================================================
      INFO [2021-14-07 03:12:01] ================ Creating 'SangerContig' ===============
      INFO [2021-14-07 03:12:01] ========================================================
      INFO [2021-14-07 03:12:01]   >> Contig Name: 'Achl_RBNII397-13'
      SUCCESS [2021-14-07 03:12:01] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:01] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:12:01] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:01]    >> 'Achl_RBNII397-13_1_F.ab1' is created (Forward Read; ABIF).
      SUCCESS [2021-14-07 03:12:02] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:02] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-14-07 03:12:02] --------------------------------------------------------
      SUCCESS [2021-14-07 03:12:02]    >> 'Achl_RBNII397-13_2_R.ab1' is created (Reverse Read; ABIF).
      INFO [2021-14-07 03:12:02]    >> The number of reads detected: 2
      Assessing frameshifts in nucleotide sequences:
      |=================================================================================================| 100%

      Time difference of 0.02 secs
      SUCCESS [2021-14-07 03:12:03] ==========================================================
      SUCCESS [2021-14-07 03:12:03] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-14-07 03:12:03] ==========================================================
      INFO [2021-14-07 03:12:03]    >> 2 read(s) created from ABIF file.
      INFO [2021-14-07 03:12:03]      >> 1 reads assigned to 'forward reads' according to 'csv file'.
      INFO [2021-14-07 03:12:03]      >> 1 reads assigned to 'reverse reads' according to 'csv file'.
      INFO [2021-14-07 03:12:03] Aligning consensus reads ... 
      SUCCESS [2021-14-07 03:12:04] #############################################################
      SUCCESS [2021-14-07 03:12:04] ######## 'SangerAlignment' S4 instance is created !! ########
      SUCCESS [2021-14-07 03:12:04] #############################################################
      INFO [2021-14-07 03:12:04]   >> 8 contigs detected from 'csv file'.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_ACHLO006-09':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_ACHLO007-09':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_ACHLO040-09':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_ACHLO041-09':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_RBNII384-13':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_RBNII395-13':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_RBNII396-13':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]       >> Contig 'Achl_RBNII397-13':
      INFO [2021-14-07 03:12:04]           >> 1 forward reads.
      INFO [2021-14-07 03:12:04]           >> 1 reverse reads.
      INFO [2021-14-07 03:12:04]   >> 16 reads created from ABIF file.
      INFO [2021-14-07 03:12:04]   >> Reads are trimmed by 'M1 - Mott’s trimming algorithm'.
      DEBUG [2021-14-07 03:12:04]    >> For more information, please run 'object'.
      DEBUG [2021-14-07 03:12:04]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads

|


(3) Updating *SangerAlignment* quality trimming parameters :sub:`(AB1)`
------------------------------------------------------------------------

.. code-block:: R

   newSangerAlignment <- updateQualityParam(my_sangerAlignment,
                                            TrimmingMethod       = "M2",
                                            M1TrimmingCutoff     = NULL,
                                            M2CutoffQualityScore = 29,
                                            M2SlidingWindowSize  = 15)

|


(4) Launching *SangerAlignment* Shiny app :sub:`(AB1)`
-------------------------------------------------------

.. code-block:: R

   launchApp(my_sangerAlignment)

|


(5) Writing *SangerAlignment* FASTA files :sub:`(AB1)`
-------------------------------------------------------

.. code-block:: R

   writeFasta(my_sangerAlignment)


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

(1) :download:`Sanger_contigs_unalignment.fa <../files/SangerAlignment_ab1/Sanger_contigs_unalignment.fa>`
(2) :download:`Sanger_contigs_alignment.fa <../files/SangerAlignment_ab1/Sanger_contigs_alignment.fa>`
(3) :download:`Sanger_all_trimmed_reads.fa <../files/SangerAlignment_ab1/Sanger_all_trimmed_reads.fa>`

|

(6) Generating *SangerAlignment* report :sub:`(AB1)`
-----------------------------------------------------

.. code-block:: R

   generateReport(my_sangerAlignment)

You can check the html report of `this SangerAlignment example (ABIF) <https://kuanhao-chao.github.io/sangeranalyseR_report/SangerAlignment/AB1/SangerAlignment/SangerAlignment_Report.html>`_.

-----

|
|
