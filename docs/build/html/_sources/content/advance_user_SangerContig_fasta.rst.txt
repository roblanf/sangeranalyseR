Advanced User Guide - *SangerContig* (**FASTA**)
================================================

*SangerContig* is in the intermediate level of sangeranalyseR (:ref:`Figure_1<SangerContig_hierarchy_fasta>`), and each *SangerContig* instance corresponds to a contig in a Sanger sequencing experiment. Among its slots, there are two lists, forward and reverse read list, storing *SangerRead* in the corresponding direction. 

In this section, we are going to go through details about a reproducible *SangerContig* analysis example with the **FASTA** file input in sangeranalyseR. By running the following example codes, you will get an end-to-end SangerContig analysis result.

.. _SangerContig_hierarchy_fasta:
.. figure::  ../image/SangerContig_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerContig* level.

|

Preparing *SangerContig* **FASTA** input
+++++++++++++++++++++++++++++++++++++++++

In :ref:`Advanced User Guide - *SangerContig* (**AB1**)`, we demonstrated how to use **AB1** input files to create *SangerContig* instance. Here, we explain another input format - the **FASTA** input. Before starting the analysis, users need to prepare one **FASTA** file, which must end with **.fa** or **.fasta**, containing sequences of all reads. In this example, the **FASTA** file is in the sangeranalyseR package, and you can simply get its path by running the following codes:

.. code-block:: R

      rawDataDir <- system.file("extdata", package = "sangeranalyseR")
      fastaFN <- file.path(rawDataDir, "fasta", "SangerContig", "Achl_ACHLO006-09.fa")


The value of :code:`fastaFN` is where the **FASTA** file is placed. If your operating system is macOS, then its value should look like this:

.. code-block:: 

      /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/fasta/SangerContig/Achl_ACHLO006-09.fa

And we showed the reads in :code:`fastaFN` in Figure_2 (:download:`example FASTA file <../files/SangerContig_fasta/Achl_ACHLO006-09.fa>`):

.. _SangerContig_fasta_input:
.. figure::  ../image/SangerContig_fasta_input.png
   :align:   center
   :scale:   45 %

   Figure 2. *SangerContig* **FASTA** input file.


Inside the **FASTA** file (:ref:`Figure_2<SangerContig_fasta_input>`; :download:`Achl_ACHLO006-09.fa <../files/SangerContig_fasta/Achl_ACHLO006-09.fa>`), the strings starting with ">" before each read are the read names. There are two ways of grouping reads which are **"regular expression matching"** and **"CSV file matching"**, and following are instructions of how to prepare your **FASTA** input file.



(1) "regular expression matching" *SangerContig* inputs (**FASTA**)
---------------------------------------------------------------------


For regular expression matching method, sangeranalyseR will group reads based on their contig name and read direction in their names automatically; therefore, users have to follow the read-naming regulations below:

.. note::

    *  All reads in the same contig group must include the same contig name in their read names.
    *  Forward or reverse direction also has to be specified in their read names.

There are four parameters, :code:`FASTA_File`, :code:`contigName`, :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`, that define the grouping rule to let sangeranalyseR automatically match correct reads in **FASTA** file and divide them into forward and reverse directions.

.. note::

  * :code:`FASTA_File`: this is the path to **FASTA** file that contains all sequences of reads, and it can be either an absolute or relative path. We suggest users to include only target reads inside this **FASTA** file and do not include any other unrelated reads.
  * :code:`contigName`: this is a regular expression that matches read names that are going to be included in the *SangerContig* analysis. :code:`grepl` function in R is used.
  * :code:`REGEX_SuffixForward`: this is a regular expression that matches all read names in forward direction. :code:`grepl` function in R is used.
  * :code:`REGEX_SuffixReverse`: this is a regular expression that matches all read names in reverse direction. :code:`grepl` function in R is used.

If you don't know what regular expression is, don't panic - it's just a way of recognising text. Please refer to :ref:`What is a regular expression?` for more details. Here is an example of how it works in sangeranalseR:


So how sangeranalyseR works is that it first matches the :code:`contigName` to exclude unrelated files and then separate the forward and reverse reads by matching :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. Therefore, it is important to make sure that all target reads in the **FASTA** file share the same :code:`contigName` and carefully select your :code:`REGEX_SuffixForward` and :code:`REGEX_SuffixReverse`. The bad file-naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate wrong results. Therefore, it is important to have a consistent naming strategy. So, how should we systematically name the reads? We suggest users to follow the file-naming regulation in :ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig_fasta>`. 


.. _sangeranalyseR_filename_convention_SangerContig_fasta:
.. figure::  ../image/sangeranalyseR_filename_convention_fasta.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested read naming regulation in **FASTA** file - *SangerContig*.

As you can see, the first part of the regulation is a consensus read name (or contig name), which narrows down the scope of reads to those we are going to examine. The second part of the regulation is an index. Since there might be more than one read that is in the forward or reverse direction, we recommend you to number your reads in the same contig group. The last part is a direction which is either 'F' (forward) or 'R' (reverse). 

To make it more specific, let's go back to the true example. In :ref:`Figure_2<SangerContig_fasta_input>`, there are two reads in the :code:`FASTA` file (:code:`fasta_FN`). First, we set :code:`contigName` to :code:`"Achl_ACHLO006-09"` to confirm that two of them, :code:`Achl_ACHLO006-09_1_F` and :code:`Achl_ACHLO006-09_2_R`, contain our target :code:`contigName` and should be included. Then, we set :code:`REGEX_SuffixForward` to :code:`"_[0-9]*_F$"` and :code:`REGEX_SuffixReverse` to :code:`"_[0-9]*_R$"` to let sangeranalyseR match and group forward and reverse reads automatically. By the regular expression rule, :code:`Achl_ACHLO006-09_1_F` and :code:`Achl_ACHLO006-09_2_R` will be categorized into "forward read list" and "reverse read list" respectively. The reason why we strongly recommend you to follow this file-naming regulation is that by doing so, you can directly adopt the example regular expression matching values, :code:`"_[0-9]*_F$"` and :code:`"_[0-9]*_R$"`, to group reads and reduce chances of error. 

After understanding how parameters work, please refer to :ref:`Creating *SangerContig* instance from **FASTA**` below to see how to create 'Achl_ACHLO006-09' *SangerContig* instance.


(2) "CSV file matching" *SangerContig* inputs (**FASTA**)
----------------------------------------------------------

No doubt that read names in the original **FASTA** file do not follow the naming regulation, and you do not want to change the original **FASTA** file; thus, we provide a second grouping approach, CSV file matching method. sangeranalyseR will group reads in the **FASTA** file based on the information in a CSV file automatically, and users do not need to alter the read names in the **FASTA** file; therefore, users have to follow the regulations below:

.. note::

    Here is an :download:`example CSV file <../files/SangerContig_fasta/names_conversion_1.csv>` (:ref:`Figure_4<sangeranalyseR_csv_file_SangerContig_fasta>`)

      .. _sangeranalyseR_csv_file_SangerContig_fasta:
      .. figure::  ../image/sangeranalyseR_csv_file_sangercontig_fasta.png
         :align:   center
         :scale:   45 %

         Figure 4. Example CSV file for *SangerContig* instance creation.  

    *  There must be three columns, "**reads**", "**direction**", and "**contig**", in the CSV file.
    *  The "**reads**" column stores the read names in the **FASTA** file that are going to be included in the analysis.
    *  The "**direction**" column stores the direction of the reads. It must be "F" (forward) or "R" (reverse).
    *  The "**contig**" column stores the contig name that each read blongs. Reads in the same contig have to have the same contig name, and they will be grouped into the same *SangerContig* instance.

There are three parameters, :code:`FASTA_File`, :code:`contigName`, and :code:`CSV_NamesConversion`,that define the grouping rule to help sangeranalseR to automatically match correct reads in a **FASTA** file and divide them into forward and reverse directions.

.. note::

  * :code:`FASTA_File`: this is the path to **FASTA** file that contains all sequences of reads, and it can be either an absolute or relative path. We suggest users to include only target reads inside this **FASTA** file and do not include any other unrelated reads.
  * :code:`contigName`: this is a regular expression that matches read names that are going to be included in the *SangerContig* analysis. :code:`grepl` function in R is used.
  * :code:`CSV_NamesConversion`: this is the path to the **CSV** file. It can be either an absolute or relative path.


The main difference between "CSV file matching" and "regular expression matching" is where the grouping rule is written. For "regular expression matching", rules are writtein in read names, and thus there are more naming requirements for users to follow. In contrast, rules of "CSV file matching" are written in an additional **CSV** file so it is more flexible on naming reads.


So how sangeranalyseR works is that it first reads in the **CSV** file (with *"reads"*, *"direction"*, and *"contig"* columns), filter out rows whose *"contig"* is not the value of :code:`contigName` parameter, find the read names in the **FASTA** file listed in *"reads"*, and assign directions to them based on *"direction"*.

To make it more specific, let's go back to the true example. First, we prepare a :download:`CSV file <../files/SangerContig_fasta/names_conversion_1.csv>` (:code:`CSV_NamesConversion`) and a :download:`FASTA file <../files/SangerContig_fasta/Achl_ACHLO006-09.fa>` (:code:`FASTA_File`). In the **CSV** file, both rows have the contig name :code:`"Achl_ACHLO006-09"`, which is what we need to assign to the :code:`contigName` parameter. sangeranalyseR then checks and matches *"reads"* of these two rows, :code:`"Achl_ACHLO006-09_1_F"` and :code:`"Achl_ACHLO006-09_2_R"`. Last, these two reads are assigned into "forward read list" and "reverse read list" respectively by the *"direction"* column.

After understanding how parameters work, please refer to :ref:`Creating *SangerContig* instance from **FASTA**` below to see how to create 'Achl_ACHLO006-09' *SangerContig* instance.



|

Creating *SangerContig* instance from **FASTA**
++++++++++++++++++++++++++++++++++++++++++++++++++

After preparing the input directory, we can create a *SangerContig* instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Their input parameters are same, and all of them have their default values. For more details about *SangerContig* inputs and slots definition, please refer to `sangeranalyseR reference manual <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. We will explain two *SangerContig* instance creation methods, "regular expression matching" and "CSV file matching".


(1) "regular expression matching" *SangerContig* creation (**FASTA**)
---------------------------------------------------------------------

The consturctor function and :code:`new` method below contain four parameters, :code:`FASTA_File`, :code:`contigName`, :code:`REGEX_SuffixForward`, and :code:`REGEX_SuffixReverse`, that we mentioned in the previous section. In contrast to **AB1** input method, it does not include quality trimming and chromatogram visualization parameters. Run the following code and create :code:`my_sangerContigFa` instance. 

.. code-block:: R

   # using `constructor` function to create SangerRead instance
   my_sangerContigFa <- SangerContig(inputSource           = "FASTA",
                                     processMethod         = "REGEX",
                                     FASTA_File            = fastaFN,
                                     contigName            =  "Achl_ACHLO006-09",
                                     REGEX_SuffixForward   = "_[0-9]*_F$",
                                     REGEX_SuffixReverse   = "_[0-9]*_R$",
                                     refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                     minReadsNum           = 2,
                                     minReadLength         = 20,
                                     minFractionCall       = 0.5,
                                     maxFractionLost       = 0.5,
                                     geneticCode           = GENETIC_CODE,
                                     acceptStopCodons      = TRUE,
                                     readingFrame          = 1,
                                     processorsNum         = 1)

   # using `new` method to create SangerRead instance
   my_sangerContigFa <- new("SangerContig",
                            inputSource           = "FASTA",
                            processMethod         = "REGEX",
                            FASTA_File            = fastaFN,
                            contigName            = "Achl_ACHLO006-09",
                            REGEX_SuffixForward   = "_[0-9]*_F$",
                            REGEX_SuffixReverse   = "_[0-9]*_R$",
                            refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                            minReadsNum           = 2,
                            minReadLength         = 20,
                            minFractionCall       = 0.5,
                            maxFractionLost       = 0.5,
                            geneticCode           = GENETIC_CODE,
                            acceptStopCodons      = TRUE,
                            readingFrame          = 1,
                            processorsNum         = 1)

In this example, :code:`contigName` is set to :code:`Achl_ACHLO006-09`, so :code:`Achl_ACHLO006-09_1_F` and :code:`Achl_ACHLO006-09_2_R` are matched and selected. Moreover, by regular expression pattern matching, :code:`Achl_ACHLO006-09_1_F` is categorized into the forward list, and :code:`Achl_ACHLO006-09_2_R` is categorized into the reverse read. Both reads are aligned into a contig, :code:`my_sangerContigFa`, and it will be used as the input for the following functions.

Inside the R shell, you can run :code:`my_sangerContigFa` to get basic information of the instance or run :code:`my_sangerContigFa@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerContigFa` is successfully created.

Here is the output of :code:`my_sangerContigFa`::

   SangerContig S4 instance
            Input Source :  FASTA 
            Process Method :  REGEX 
         Fasta File Name :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/fasta/SangerContig/Achl_ACHLO006-09.fa 
      REGEX Suffix Forward :  _[0-9]*_F$ 
      REGEX Suffix Reverse :  _[0-9]*_R$ 
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
   SUCCESS [2021-13-07 11:52:40] 'Achl_ACHLO006-09' is successfully created!

Here is the output of :code:`my_sangerContigFa@objectResults@readResultTable`::


               readName creationResult errorType errorMessage inputSource    direction
   1 Achl_ACHLO006-09_1_F           TRUE      None         None       FASTA Forward Read
   2 Achl_ACHLO006-09_2_R           TRUE      None         None       FASTA Reverse Read



(2) "CSV file matching" *SangerContig* creation (**FASTA**)
------------------------------------------------------------
The consturctor function and :code:`new` method below contain three parameters, :code:`FASTA_File`, :code:`contigName`, and :code:`CSV_NamesConversion`, that we mentioned in the previous section. Run the following code and create :code:`my_sangerContigFa` instance. 

.. code-block:: R

   csv_namesConversion <- file.path(rawDataDir, "fasta", "SangerContig", "names_conversion_1.csv")

   # using `constructor` function to create SangerRead instance
   my_sangerContigFa <- SangerContig(inputSource           = "FASTA",
                                     processMethod         = "CSV",
                                     FASTA_File            = fastaFN,
                                     contigName            = "Achl_ACHLO006-09",
                                     CSV_NamesConversion   = csv_namesConversion,
                                     refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                     minReadsNum           = 2,
                                     minReadLength         = 20,
                                     minFractionCall       = 0.5,
                                     maxFractionLost       = 0.5,
                                     geneticCode           = GENETIC_CODE,
                                     acceptStopCodons      = TRUE,
                                     readingFrame          = 1,
                                     processorsNum         = 1)

   # using `new` method to create SangerRead instance
   my_sangerContigFa <- new("SangerContig",
                            inputSource           = "FASTA",
                            processMethod         = "CSV",
                            FASTA_File            = fastaFN,
                            contigName            = "Achl_ACHLO006-09",
                            CSV_NamesConversion   = csv_namesConversion,
                            refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                            minReadsNum           = 2,
                            minReadLength         = 20,
                            minFractionCall       = 0.5,
                            maxFractionLost       = 0.5,
                            geneticCode           = GENETIC_CODE,
                            acceptStopCodons      = TRUE,
                            readingFrame          = 1,
                            processorsNum         = 1)

First, you need to load the **CSV** file into the R environment. If you are still don't know how to prepare it, please check :ref:`(2) "CSV file matching" *SangerContig* inputs (**FASTA**)`. Then, it will follow rules in the CSV file and create :code:`my_sangerContigFa`. After it's created, inside the R shell, you can run :code:`my_sangerContigFa` to get basic information of the instance or run :code:`my_sangerContigFa@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`my_sangerContigFa` is successfully created.

Here is the output of :code:`my_sangerContigFa`::

   SangerContig S4 instance
            Input Source :  FASTA 
            Process Method :  CSV 
         Fasta File Name :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/fasta/SangerContig/Achl_ACHLO006-09.fa 
      CSV Names Conversion :  /Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata/fasta/SangerContig/names_conversion_1.csv 
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
   SUCCESS [2021-13-07 12:01:57] 'Achl_ACHLO006-09' is successfully created!


Here is the output of :code:`my_sangerContigFa@objectResults@readResultTable`::

               readName creationResult errorType errorMessage inputSource    direction
   1 Achl_ACHLO006-09_1_F           TRUE      None         None       FASTA Forward Read
   2 Achl_ACHLO006-09_2_R           TRUE      None         None       FASTA Reverse Read


|


Writing *SangerContig* FASTA files :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++++++
Users can write the *SangerContig* instance, :code:`my_sangerContigFa`, to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`reads_unalignment`: Writing reads into a single **FASTA** file (only trimmed without alignment).
* :code:`reads_alignment`: Writing reads alignment and contig read to a single **FASTA** file.
* :code:`contig`: Writing the contig to a single **FASTA** file.
* :code:`all`: Writing reads, reads alignment, and the contig into three different files.

Below is the oneliner for writing out **FASTA** files. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFasta` function.

.. code-block:: R

   writeFasta(my_sangerContigFa,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA,
              selection         = "all")

Users can download the output FASTA file of this example through the following three links:

(1) :download:`Achl_ACHLO006-09_reads_unalignment.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_reads_unalignment.fa>`
(2) :download:`Achl_ACHLO006-09_reads_alignment.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_reads_alignment.fa>`
(3) :download:`Achl_ACHLO006-09_contig.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_contig.fa>`


|

Generating *SangerContig* report :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save *SangerContig* instance, :code:`my_sangerContigFa`, into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.

Users can set :code:`includeSangerRead` parameter to decide to which level the *SangerContig* report will go. Moreover, after the reports are generated,
users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it will take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReport(my_sangerContigFa,
                  outputDir           = tempdir(),
                  includeSangerRead   = FALSE)

Users can access to '*Basic Information*', '*SangerContig Input Parameters*', '*Contig Sequence*' and '*Contig Results*' sections inside the generated `SangerContig html report of this example <https://howardchao.github.io/sangeranalyseR_report/SangerContig/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]/SangerContig_Report.html>`_. Furthermore, users can also navigate through html reports of all forward and reverse *SangerRead* in this *SangerContig* report.



-----

|
|

Code summary (*SangerContig*, **FASTA**)
++++++++++++++++++++++++++++++++++++++++++++++++++


1. Preparing *SangerContig* **FASTA** input
---------------------------------------------

.. code-block:: R

      rawDataDir <- system.file("extdata", package = "sangeranalyseR")
      fastaFN <- file.path(rawDataDir, "fasta", "SangerContig", "Achl_ACHLO006-09.fa")

|

2. Creating *SangerContig* instance from **FASTA**
----------------------------------------------------

.. code-block:: R

   # using `constructor` function to create SangerRead instance
   my_sangerContigFa <- SangerContig(inputSource           = "FASTA",
                                     processMethod         = "REGEX",
                                     FASTA_File            = fastaFN,
                                     contigName            =  "Achl_ACHLO006-09",
                                     REGEX_SuffixForward   = "_[0-9]*_F$",
                                     REGEX_SuffixReverse   = "_[0-9]*_R$",
                                     refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")

   # using `new` method to create SangerRead instance
   my_sangerContigFa <- new("SangerContig",
                            inputSource           = "FASTA",
                            processMethod         = "REGEX",
                            FASTA_File            = fastaFN,
                            contigName            = "Achl_ACHLO006-09",
                            REGEX_SuffixForward   = "_[0-9]*_F$",
                            REGEX_SuffixReverse   = "_[0-9]*_R$",
                            refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")

.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

      INFO [2021-13-07 12:18:32] ========================================================
      INFO [2021-13-07 12:18:32] ================ Creating 'SangerContig' ===============
      INFO [2021-13-07 12:18:32] ========================================================
      INFO [2021-13-07 12:18:32]   >> Contig Name: 'Achl_ACHLO006-09'
      INFO [2021-13-07 12:18:32]   >> You are using Regular Expression Method to group reads in FASTA file (No CSV file)!
      INFO [2021-13-07 12:18:32] >> Your contig name is Achl_ACHLO006-09
      SUCCESS [2021-13-07 12:18:32] --------------------------------------------------------
      SUCCESS [2021-13-07 12:18:32] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-13-07 12:18:32] --------------------------------------------------------
      SUCCESS [2021-13-07 12:18:32]    >> 'Achl_ACHLO006-09_1_F' is created (Forward Read; FASTA).
      SUCCESS [2021-13-07 12:18:32] --------------------------------------------------------
      SUCCESS [2021-13-07 12:18:32] -------- 'SangerRead' S4 instance is created !! --------
      SUCCESS [2021-13-07 12:18:32] --------------------------------------------------------
      SUCCESS [2021-13-07 12:18:32]    >> 'Achl_ACHLO006-09_2_R' is created (Reverse Read; FASTA).
      INFO [2021-13-07 12:18:32]    >> The number of reads detected: 2
      INFO [2021-13-07 12:18:32] Correcting frameshifts in reads using amino acidreference sequence
      Assessing frameshifts in nucleotide sequences:
      |=====================================================================================| 100%

      Time difference of 0.05 secs
      SUCCESS [2021-13-07 12:18:34] ==========================================================
      SUCCESS [2021-13-07 12:18:34] ======== 'SangerContig' S4 instance is created !! ========
      SUCCESS [2021-13-07 12:18:34] ==========================================================
      INFO [2021-13-07 12:18:34]    >> 2 read(s) created from FASTA file.
      INFO [2021-13-07 12:18:34]      >> 1 reads assigned to 'forward reads' according to 'regular expression'.
      INFO [2021-13-07 12:18:34]      >> 1 reads assigned to 'reverse reads' according to 'regular expression'.
      DEBUG [2021-13-07 12:18:34]    >> For more information, please run 'object'
      DEBUG [2021-13-07 12:18:34]    >> Run 'object@objectResults@readResultTable' to check the results of each Sanger reads

|


3. Writing *SangerContig* FASTA files :sub:`(FASTA)`
----------------------------------------------------

.. code-block:: R

   writeFasta(my_sangerContigFa)


.. container:: toggle

     .. container:: header

        Following is the R shell output that you will get.

     .. code-block::

         INFO [2021-13-07 12:19:16] Your input is 'SangerContig' S4 instance
         INFO [2021-13-07 12:19:16] >>> outputDir : /private/var/folders/33/7v38jdjd2874jcxb6l71m00h0000gn/T/RtmpPlQHzO
         INFO [2021-13-07 12:19:16] Start to write 'Achl_ACHLO006-09' to FASTA format ...
         INFO [2021-13-07 12:19:16] >> Writing alignment to FASTA ...
         INFO [2021-13-07 12:19:16] >> Writing all single reads to FASTA ...
         INFO [2021-13-07 12:19:16] >> Writing consensus read to FASTA ...
         INFO [2021-13-07 12:19:16] Finish writing 'Achl_ACHLO006-09' to FASTA format

|

And you will get three FASTA files:

(1) :download:`Achl_ACHLO006-09_reads_unalignment.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_reads_unalignment.fa>`
(2) :download:`Achl_ACHLO006-09_reads_alignment.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_reads_alignment.fa>`
(3) :download:`Achl_ACHLO006-09_contig.fa <../files/SangerContig_fasta/Achl_ACHLO006-09_contig.fa>`

|

4. Generating *SangerContig* report :sub:`(FASTA)`
---------------------------------------------------

.. code-block:: R

   generateReport(my_sangerContigFa)

-----

|
|
