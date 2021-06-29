Advanced User Guide - *SangerRead* (**FASTA**)
==============================================

*SangerRead* is the bottommost level in sangeranalyseR (:ref:`Figure_1<SangerRead_hierarchy_fasta>`) and each *SangerRead* object corresponds to a single read in Sanger sequencing. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerRead level* with **FASTA** file input.

.. _SangerRead_hierarchy_fasta:
.. figure::  ../image/SangerRead_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerRead* level.

|

Preparing *SangerRead* **FASTA** input
+++++++++++++++++++++++++++++++++++++++
We design the **FASTA** file input for those who do not want to do quality trimming and base calling for their *SangerRead*; therefore, it does not contain quality trimming and chromatogram input parameters and results in its slots. Before starting the analysis, users need to prepare one target **FASTA** file. The only hard regulation of the filename is that file extension must be **.fasta** or **.fa**.

|

Creating *SangerRead* instance from **FASTA**
++++++++++++++++++++++++++++++++++++++++++++++

After preparing the *SangerRead* input **FASTA** file, the next step is to create the *SangerRead* S4 instance by running :code:`SangerRead` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method which makes instance creation more intuitive. Most of the input parameters have their own default values. In the constructor below, we list important parameters. The filename is assigned to :code:`readFileName`. Inside **FASTA** file, the string in the first line after ">" is the name of the read. Users need to assign the read name to :code:`fastaReadName` which is used to match the target read in **FASTA** input file. :ref:`Figure 2<SangerRead_fasta_input_file>` is a valid **FASTA** file and the value of :code:`fastaReadName` is :code:`Achl_ACHLO006-09_1_Forward`.

.. code-block:: R

   sangerReadFfa <- new("SangerRead",
                        inputSource          = "FASTA",
                        readFeature          = "Forward Read",
                        readFileName         = "ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F.fa",
                        fastaReadName        = "Achl_ACHLO006-09_1_Forward",
                        geneticCode          = GENETIC_CODE)

.. _SangerRead_fasta_input_file:
.. figure::  ../image/SangerRead_fasta_input_file.png
   :align:   center
   :scale:   40 %

   Figure 2. *SangerRead* **FASTA** input file.

The inputs of :code:`SangerRead` constructor function and :code:`new` method are same. For more details about *SangerRead* inputs and slots definition, please refer to `sangeranalyseR reference manual (need update after upload function manul) <http://packages.python.org/an_example_pypi_project/>`_.

|


Writing *SangerRead* FASTA files :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++++
Users can write the *SangerRead* instance to **FASTA** files. Because the **FASTA** input does not support quality trimming and base calling, in this example, the sequence of the written **FASTA** file will be same as the input **FASTA** file. Moreover, users can set the compression level through the one-line function :code:`writeFasta` which mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package.

.. code-block:: R

   writeFasta(sangerReadFfa,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA)

Users can download the `output FASTA file <https://howardchao.github.io/sangeranalyseR_report/SangerRead/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F.fa>`_ of this example.


|

Generating *SangerRead* report :sub:`(FASTA)`
++++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save *SangerRead* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files. The results in the report are static.

.. code-block:: R

   generateReport(sangerReadFfa,
                  outputDir = tempdir())

`SangerRead_Report_fasta.html <https://howardchao.github.io/sangeranalyseR_report/SangerRead/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F/SangerRead_Report_fasta.html>`_ is the generated *SangerRead* report html of this example. Users can access to '*Basic Information*', '*DNA Sequence*' and '*Amino Acids Sequence*' sections inside this report.

-----

|
|










A Reproducible Example (*SangerRead*, **fasta**)
+++++++++++++++++++++++++++++++++++++++++++++++++


1. Preparing *SangerRead* **FASTA** input
------------------------------------------
The data of this example is in the sangeranalyseR package; thus, you can simply get its path from the library.

.. code-block:: R

   inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
   A_chloroticaFFNfa <- file.path(inputFilesPath,
                                  "fasta",
                                  "SangerRead",
                                  "Achl_ACHLO006-09_1_F.fa")
   readNameFfa <- "Achl_ACHLO006-09_1_F"

|

2. Creating *SangerRead* instance from **FASTA**
-------------------------------------------------
Run the following on-liner to create the *SangerRead* object.


.. code-block:: R

   # using `constructor` function to create SangerRead instance
   sangerReadFfa <- SangerRead(inputSource        = "FASTA",
                               readFeature        = "Forward Read",
                               readFileName       = A_chloroticaFFNfa,
                               fastaReadName      = readNameFfa)
   
   # using `new` method to create SangerRead instance
   sangerReadFfa <- new("SangerRead",
                        inputSource        = "FASTA",
                        readFeature        = "Forward Read",
                        readFileName       = A_chloroticaFFNfa,
                        fastaReadName      = readNameFfa)


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


3. Writing *SangerRead* FASTA files :sub:`(FASTA)`
---------------------------------------------------

Write the read into a FASTA file.

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

(1) :download:`Sanger_all_trimmed_reads.fa <../files/SangerRead_fasta/Achl_ACHLO006-09_1_F.fa>`

|

4. Generating *SangerRead* report :sub:`(FASTA)`
-------------------------------------------------

Last but not least, generate an Rmarkdown report to store all the sequence information.

.. code-block:: R

   generateReport(sangerReadFfa)


-----

|
|
