Advanced User Guide - *SangerRead* (**AB1**)
============================================

*SangerRead* is in the bottommost level of sangeranalyseR (:ref:`Figure_1<SangerRead_hierarchy>`), and each *SangerRead* object corresponds to a single read (one **AB1** file) in a Sanger sequencing experiment. *SangerRead* class extends *sangerseq* class from `sangerseqR <https://www.bioconductor.org/packages/release/bioc/html/sangerseqR.html>`_ package and contains input parameters and results of quality trimming and chromatogram. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerRead level* with **AB1** file input.

.. _SangerRead_hierarchy:
.. figure::  ../image/SangerRead_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerRead* level.


Preparing *SangerRead* **AB1** input
+++++++++++++++++++++++++++++++++++++
The main input file format to create *SangerRead* instance is **AB1**. Before starting the analysis, users need to prepare one target **AB1** file. The only hard regulation of the filename is that the input file must have **.ab1** as its file extension. There are some suggestions about the filename in the note below:

.. note::

    * **AB1** file should be indexed for better consistency with file-naming regulation for *SangerContig* and *SangerAlignment*.
    * Forward or reverse direction should be specified in the filename.

:ref:`Figure_2<SangerRead_file_structure>` shows the suggested file-naming strategy. The filename should contain four main parts: **"Contig name"**, **"Index number"**, **"Direction"** and **"ab1 file extension"**.

* **"Contig name"** :  :code:`Achl_RBNII397-13`
* **"Index number"** :  :code:`1`
* **"Direction"** :  :code:`F`
* **"ab1 file extension"** :  :code:`.ab1`

.. _SangerRead_file_structure:
.. figure::  ../image/SangerRead_file_structure.png
   :align:   center
   :scale:   80 %

   Figure 2. *SangerRead* filename regulation.

In *SangerRead* section, it is not compulsory to follow the file-naming regulation because users can directly specify the filename in input (see :ref:`Creating *SangerRead* instance from **AB1**`); however, in the *SangerContig* and *SangerAlignment*, sangeranalyseR will automatically group files, so it is compulsory to have systematic file-naming strategy. For more details, please read :ref:`Advanced User Guide - *SangerContig* (**AB1**)` and :ref:`Advanced User Guide - *SangerAlignment* (**AB1**)`. :ref:`Figure_3<sangeranalyseR_filename_convention_SangerRead>` shows the suggested **AB1** file-naming regulation.


.. _sangeranalyseR_filename_convention_SangerRead:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file-naming regulation - *SangerRead*.


|

Creating *SangerRead* instance from **AB1**
++++++++++++++++++++++++++++++++++++++++++++

After preparing the *SangerRead* input **AB1** file, the next step is to create the *SangerRead* instance by running :code:`SangerRead` constructor function or :code:`new` method. The constructor function is a wrapper for the :code:`new` method which makes instance creation more intuitive. The inputs include **Basic Parameters**, **Trimming Parameters**, and **Chromatogram Parameters**, and all of them have default values. In the example below, we show both *SangerRead* creation methods with important parameters.

.. code-block:: R

   # using `constructor` function to create SangerRead instance
   sangerReadF <- SangerRead(inputSource           = "ABIF",
                             readFeature           = "Forward Read",
                             readFileName          = "Achl_RBNII397-13_1_F.ab1",
                             geneticCode           = GENETIC_CODE,
                             TrimmingMethod        = "M1",
                             M1TrimmingCutoff      = 0.0001,
                             M2CutoffQualityScore  = NULL,
                             M2SlidingWindowSize   = NULL,
                             baseNumPerRow         = 100,
                             heightPerRow          = 200,
                             signalRatioCutoff     = 0.33,
                             showTrimmed           = TRUE)

   # using `new` method to create SangerRead instance
   sangerReadF <- new("SangerRead",
                       inputSource           = "ABIF",
                       readFeature           = "Forward Read",
                       readFileName          = "Achl_RBNII397-13_1_F.ab1",
                       geneticCode           = GENETIC_CODE,
                       TrimmingMethod        = "M1",
                       M1TrimmingCutoff      = 0.0001,
                       M2CutoffQualityScore  = NULL,
                       M2SlidingWindowSize   = NULL,
                       baseNumPerRow         = 100,
                       heightPerRow          = 200,
                       signalRatioCutoff     = 0.33,
                       showTrimmed           = TRUE)


The inputs of :code:`SangerRead` constructor function and :code:`new` method are the same. For more details about *SangerRead* inputs and slots definition, please refer to the `sangeranalyseR reference manual <https://bioconductor.org/packages/release/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. The created *SangerRead* instance, :code:`sangerReadF`, is used as the input for the following functions.

|

Visualizing *SangerRead* trimmed read
++++++++++++++++++++++++++++++++++++++
Before going to :ref:`Writing *SangerRead* FASTA file :sub:\`(AB1)\`` and :ref:`Generating *SangerRead* report :sub:\`(AB1)\`` pages, it is suggested to visualize the trimmed *SangerRead*. Run the :code:`qualityBasePlot` function to get the result in :ref:`Figure_4 <SangerRead_qualityBasePlot>`. It shows the quality score for each base pairs and the trimming start/end points of the sequence.


.. _SangerRead_qualityBasePlot:
.. figure::  ../image/SangerRead_qualityBasePlot.png
   :align:   center
   :scale:   30 %

   Figure 4. *SangerRead* trimmed read visualization.

.. code-block:: R

   qualityBasePlot(sangerReadF)

|

Updating *SangerRead* quality trimming parameters
++++++++++++++++++++++++++++++++++++++++++++++++++
In the previous :ref:`Creating *SangerRead* instance from **AB1**` part, the constructor function applies the quality trimming parameters to the read. These parameters are not fixed. After instance creation, users can run :code:`updateQualityParam` function which will change the *QualityReport* instance inside the *SangerRead* and update frameshift amino acid sequences.

.. code-block:: R

   newSangerRead <- updateQualityParam(sangerReadF,
                                       TrimmingMethod       = "M2",
                                       M1TrimmingCutoff     = NULL,
                                       M2CutoffQualityScore = 29,
                                       M2SlidingWindowSize  = 15)

|



Writing *SangerRead* FASTA file :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++++

After quality trimming, users can write :code:`sangerReadF` into a **FASTA** file. Below is the one-liner that needs to be run. This function, :code:`writeFasta`, mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can further set the compression level through it.

.. code-block:: R

   writeFasta(newSangerRead,
              outputDir         = tempdir(),
              compress          = FALSE,
              compression_level = NA)

Users can download the :download:`output FASTA file <../files/SangerRead_ab1/Achl_RBNII384-13_1_F.fa>` of this example.

|


Generating *SangerRead* report :sub:`(AB1)`
++++++++++++++++++++++++++++++++++++++++++++
Last but not least, users can save :code:`sangerReadF` into a static **HTML** report by knitting **Rmd** files. In this example, :code:`tempdir` function will generate a random path.

.. code-block:: R

   generateReport(newSangerRead,
                  outputDir = tempdir())

.. `SangerRead_Report_ab1.html <https://howardchao.github.io/sangeranalyseR_report/SangerRead/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F/SangerRead_Report_ab1.html>`_ is the generated *SangerRead* report html of this example. Users can access to '*Basic Information*', '*DNA Sequence*', '*Amino Acids Sequence*', '*Quality Trimming*' and '*Chromatogram*' sections inside this report.

-----

|
|




















A Reproducible Example (*SangerRead*, **ab1**)
+++++++++++++++++++++++++++++++++++++++++++++++


(1) Preparing *SangerRead* **AB1** input
----------------------------------------
The data of this example is in the sangeranalyseR package; thus, you can simply get its path from the library.

.. code-block:: R

   inputFilesPath <- system.file("extdata/", package = "sangeranalyseR")
   A_chloroticaFFN <- file.path(inputFilesPath,
                                "Allolobophora_chlorotica",
                                "ACHLO",
                                "Achl_ACHLO006-09_1_F.ab1")

|

(2) Creating *SangerRead* instance from **AB1**
-----------------------------------------------
Run the following on-liner, SangerRead :code:`constructor` or :code:`new` method, to create the *SangerRead* object.


.. code-block:: R

   # using `constructor` function to create SangerRead instance
   sangerReadF <- SangerRead(readFeature           = "Forward Read",
                             readFileName          = A_chloroticaFFN)

   # using `new` method to create SangerRead instance
   sangerReadF <- new("SangerRead",
                      readFeature           = "Forward Read",
                      readFileName          = A_chloroticaFFN)


.. container:: toggle

    .. container:: header

        Following is the R shell output that you will get.
    .. code-block::

         INFO [2021-29-06 16:28:39] ------------------------------------------------
         INFO [2021-29-06 16:28:39] -------- Creating 'SangerRead' instance --------
         INFO [2021-29-06 16:28:39] ------------------------------------------------
         INFO [2021-29-06 16:28:39] >> Forward Read: Creating abif & sangerseq ...
         INFO [2021-29-06 16:28:39]     >> Creating Forward Read raw abif ...
         INFO [2021-29-06 16:28:39]     >> Creating Forward Read raw sangerseq ...
         INFO [2021-29-06 16:28:39]           * Making basecall !!
         INFO [2021-29-06 16:28:40]           * Updating slots in 'SangerRead' instance !!
         SUCCESS [2021-29-06 16:28:40] --------------------------------------------------------
         SUCCESS [2021-29-06 16:28:40] -------- 'SangerRead' S4 instance is created !! --------
         SUCCESS [2021-29-06 16:28:40] --------------------------------------------------------
         SUCCESS [2021-29-06 16:28:40]    >> 'Achl_ACHLO006-09_1_F.ab1' is created (Forward Read; ABIF).
         INFO [2021-29-06 16:28:40]    >> Read is trimmed by 'M1 - Mott’s trimming algorithm'.
         DEBUG [2021-29-06 16:28:40]    >> For more information, please run 'object'.
         DEBUG [2021-29-06 16:28:40]    >> Run 'object@objectResults@readResultTable' to check the result of the Sanger read

|

(3) Visualizing *SangerRead* trimmed read
-----------------------------------------

Launch an interactive plotly plot to check the trimmed read.

.. code-block:: R

   qualityBasePlot(sangerReadF)

|


(4) Writing *SangerRead* FASTA file :sub:`(AB1)`
-------------------------------------------------

Write the trimmed read into a FASTA file.

.. code-block:: R

   writeFasta(sangerReadF)


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

(1) :download:`Achl_ACHLO006-09_1_F.fa <../files/SangerRead_ab1/Achl_ACHLO006-09_1_F.fa>`

|

(5) Generating *SangerRead* report :sub:`(AB1)`
-----------------------------------------------

Last but not least, generate an Rmarkdown report to store all sequence information.

.. code-block:: R

   generateReport(sangerReadF)


-----

|
|
