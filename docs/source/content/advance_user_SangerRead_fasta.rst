Advanced User Guide - *SangerRead* (**FASTA**)
==============================================

*SangerRead* is in the bottommost level of sangeranalyseR (:ref:`Figure_1<SangerRead_hierarchy_fasta>`), and each *SangerRead* object corresponds to a single read in Sanger sequencing. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerRead level* with **FASTA** file input.


.. _SangerRead_hierarchy_fasta:
.. figure::  ../image/SangerRead_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerRead* level.

|

Preparing *SangerRead* **FASTA** input
+++++++++++++++++++++++++++++++++++++++

The **FASTA** input method is designed for those who do not want to do quality trimming and base calling on their Sanger sequencing data; therefore, no quality trimming and chromatogram input parameters are needed. Before starting the analysis, users need to prepare a **FASTA** file, and in this example, it is in the sangeranalyseR package; thus, you can simply get its path by running the following codes:

.. code-block:: R

   inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
   A_chloroticaFFNfa <- file.path(inputFilesPath,
                                  "fasta",
                                  "SangerRead",
                                  "Achl_ACHLO006-09_1_F.fa")






The only hard regulation of the filename, :code:`Achl_ACHLO006-09_1_F.fa` in this example, is that file extension must be **.fasta** or **.fa**.


|

Creating *SangerRead* instance from **FASTA**
++++++++++++++++++++++++++++++++++++++++++++++

After preparing an input **FASTA** file, the next step is to create a *SangerRead* instance by running :code:`SangerRead` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method which makes instance creation more intuitive. All of the input parameters have their default values. We list important parameters in the two *SangerRead* creation methods below. :code:`readFileName` stores the **FASTA** filename, and inside it, the string in the first line after ">" is the name of the read. Users need to assign the name of the read to :code:`fastaReadName` which is used for read-matching. :ref:`Figure 2<SangerRead_fasta_input_file>` is a valid **FASTA** file, :code:`Achl_ACHLO006-09_1_F.fa` (:download:`example FASTA file <../files/SangerRead_fasta/Achl_ACHLO006-09_1_F.fa>`), and the value of :code:`fastaReadName` is :code:`Achl_ACHLO006-09_1_F`.

.. _SangerRead_fasta_input_file:
.. figure::  ../image/SangerRead_fasta_input_file.png
   :align:   center
   :scale:   40 %

   Figure 2. *SangerRead* **FASTA** input file.


.. code-block:: R

   # using `constructor` function to create SangerRead instance
   sangerReadFfa <- SangerRead(inputSource        = "FASTA",
                               readFeature        = "Forward Read",
                               readFileName       = A_chloroticaFFNfa,
                               fastaReadName      = "Achl_ACHLO006-09_1_F",
                               geneticCode        = GENETIC_CODE)


   # using `new` method to create SangerRead instance
   sangerReadFfa <- new("SangerRead",
                        inputSource        = "FASTA",
                        readFeature        = "Forward Read",
                        readFileName       = A_chloroticaFFNfa,
                        fastaReadName      = "Achl_ACHLO006-09_1_F", 
                        geneticCode        = GENETIC_CODE)



The inputs of :code:`SangerRead` constructor function and :code:`new` method are the same. For more details about *SangerRead* inputs and slots definition, please refer to `sangeranalyseR reference manual <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_.


Inside the R shell, you can run :code:`sangerReadFfa` to get basic information of the instance or run :code:`sangerReadFfa@objectResults@readResultTable` to check the creation result of every Sanger read after :code:`sangerReadFfa` is successfully created.

Here is the output of :code:`sangerReadFfa`::

   SangerRead S4 instance
            Input Source :  FASTA 
            Read Feature :  Forward Read 
            Read FileName :  Achl_ACHLO006-09_1_F.fa 
         Fasta Read Name :  Achl_ACHLO006-09_1_F 
         Primary Sequence :  CTGGGCGTCTGAGCAGGAATGGTTGGAGCCGGTATAAGACTTCTAATTCGAATCGAGCTAAGACAACCAGGAGCGTTCCTGGGCAGAGACCAACTATACAATACTATCGTTACTGCACACGCATTTGTAATAATCTTCTTTCTAGTAATGCCTGTATTCATCGGGGGATTCGGAAACTGGCTTTTACCTTTAATACTTGGAGCCCCCGATATAGCATTCCCTCGACTCAACAACATGAGATTCTGACTACTTCCCCCATCACTGATCCTTTTAGTGTCCTCTGCGGCGGTAGAAAAAGGCGCTGGTACGGGGTGAACTGTTTATCCGCCTCTAGCAAGAAATCTTGCCCACGCAGGCCCGTCTGTAGATTTAGCCATCTTTTCCCTTCATTTAGCGGGTGCGTCTTCTATTCTAGGGGCTATTAATTTTATCACCACAGTTATTAATATGCGTTGAAGAGG 
   SUCCESS [2021-12-07 23:37:43] 'Achl_ACHLO006-09_1_F.fa' is successfully created!


Here is the output of :code:`sangerReadFfa@objectResults@readResultTable`::
   
               readName creationResult errorType errorMessage inputSource    direction
   1 Achl_ACHLO006-09_1_F           TRUE      None         None       FASTA Forward Read



|


Writing *SangerRead* FASTA files :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++++
Users can write :code:`sangerReadFfa` to a **FASTA** file. Because the **FASTA** input method does not support quality trimming or base calling, in this example, the sequence of the output **FASTA** file will be the same as the input **FASTA** file. Moreover, users can set the compression level through the one-liner, :code:`writeFasta`, which mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package.


.. code-block:: R

   writeFasta(sangerReadFfa,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA)

Users can download the :download:`Achl_ACHLO006-09_1_F.fa <../files/SangerRead_fasta/Achl_ACHLO006-09_1_F.fa>` of this example.


|

Generating *SangerRead* report :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save :code:`sangerReadFfa` into a static **HTML** report by knitting **Rmd** files. In this example, :code:`tempdir` function will generate a random path.


.. code-block:: R

   generateReport(sangerReadFfa,
                  outputDir = tempdir())

.. `SangerRead_Report_fasta.html <https://howardchao.github.io/sangeranalyseR_report/SangerRead/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F/SangerRead_Report_fasta.html>`_ is the generated *SangerRead* report html of this example. Users can access to '*Basic Information*', '*DNA Sequence*' and '*Amino Acids Sequence*' sections inside this report.

.. It might seem a bit confusing why we go through all troubles to create an R class only for the purpose of  storing the **FASTA** file. That is because when it comes to the higher levels, *SangerContig* and *SangerAlignment*, we can build upon *SangerRead* class and do further analyses like sequence alignment, frameshifts correcting, and so on. Please refer to :ref:`Advanced User Guide - *SangerContig* (**FASTA**)` and :ref:`Advanced User Guide - *SangerAlignment* (**FASTA**)` to see how to start the Sanger sequencing analysis with **FASTA** files in a higher level.

-----

|
|








Code summary (*SangerRead*, **fasta**)
+++++++++++++++++++++++++++++++++++++++++++++++++


(1) Preparing *SangerRead* **FASTA** input
------------------------------------------

.. code-block:: R

   inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
   A_chloroticaFFNfa <- file.path(inputFilesPath,
                                  "fasta",
                                  "SangerRead",
                                  "Achl_ACHLO006-09_1_F.fa")

|

(2) Creating *SangerRead* instance from **FASTA**
-------------------------------------------------


.. code-block:: R

   # using `constructor` function to create SangerRead instance
   sangerReadFfa <- SangerRead(inputSource        = "FASTA",
                               readFeature        = "Forward Read",
                               readFileName       = A_chloroticaFFNfa,
                               fastaReadName      = "Achl_ACHLO006-09_1_F")

   # using `new` method to create SangerRead instance
   sangerReadFfa <- new("SangerRead",
                        inputSource        = "FASTA",
                        readFeature        = "Forward Read",
                        readFileName       = A_chloroticaFFNfa,
                        fastaReadName      = "Achl_ACHLO006-09_1_F")


.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

         INFO [2021-29-06 17:07:40] ------------------------------------------------
         INFO [2021-29-06 17:07:40] -------- Creating 'SangerRead' instance --------
         INFO [2021-29-06 17:07:40] ------------------------------------------------
         INFO [2021-29-06 17:07:40] Forward Read: Creating SangerRead from FASTA ...
         SUCCESS [2021-29-06 17:07:41] --------------------------------------------------------
         SUCCESS [2021-29-06 17:07:41] -------- 'SangerRead' S4 instance is created !! --------
         SUCCESS [2021-29-06 17:07:41] --------------------------------------------------------
         SUCCESS [2021-29-06 17:07:41]    >> 'Achl_ACHLO006-09_1_F' is created (Forward Read; FASTA).
         INFO [2021-29-06 17:07:41]    >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.
         DEBUG [2021-29-06 17:07:41]    >> For more information, please run 'object'.
         DEBUG [2021-29-06 17:07:41]    >> Run 'object@objectResults@readResultTable' to check the result of the Sanger read

|


(3) Writing *SangerRead* FASTA files :sub:`(FASTA)`
---------------------------------------------------


.. code-block:: R

   writeFasta(sangerReadFfa)


.. container:: toggle

     .. container:: header

        Following is the R shell output that you will get.

     .. code-block::

         INFO [2021-29-06 16:30:17] Your input is 'SangerRead' S4 instance
         INFO [2021-29-06 16:30:17] >>> outputDir : /private/var/folders/33/7v38jdjd2874jcxb6l71m00h0000gn/T/RtmpRAPaMV
         INFO [2021-29-06 16:30:17] Start writing '/Library/Frameworks/R.framework/Versions/4.0/Resources/library/sangeranalyseR/extdata//Allolobophora_chlorotica/ACHLO/Achl_ACHLO006-09_1_F.ab1' to FASTA format ...
         INFO [2021-29-06 16:30:17] >> '/private/var/folders/33/7v38jdjd2874jcxb6l71m00h0000gn/T/RtmpRAPaMV/Achl_ACHLO006-09_1_F.fa' is written
         [1] "/private/var/folders/33/7v38jdjd2874jcxb6l71m00h0000gn/T/RtmpRAPaMV/Achl_ACHLO006-09_1_F.fa"

|

And you will get one FASTA file:

(1) :download:`Achl_ACHLO006-09_1_F.fa <../files/SangerRead_fasta/Achl_ACHLO006-09_1_F.fa>`


|

(4) Generating *SangerRead* report :sub:`(FASTA)`
-------------------------------------------------


.. code-block:: R

   generateReport(sangerReadFfa)


-----

|
|
