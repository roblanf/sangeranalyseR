Advanced User Guide - *SangerRead* (**AB1**)
============================================

*SangerRead* is the lowest level in sangeranalyseR showed in :ref:`Figure_1<SangerRead_hierarchy>` which corresponds to a single read (one **AB1** file) in Sanger sequencing. It extends *sangerseq* S4 class from `sangerseqR <https://www.bioconductor.org/packages/release/bioc/html/sangerseqR.html>`_ package and contains quality trimming as well as chromatogram input parameters and results. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerRead level* with **AB1** file input.

.. _SangerRead_hierarchy:
.. figure::  ../image/SangerRead_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerRead* level.


Preparing *SangerRead* **AB1** input
------------------------------------
The main input file format to create *SangerRead* instance is **AB1**. Before starting the analysis, users need to prepare one target **AB1** file. The only hard regulation of the filename is that the input file must have **.ab1** as its file extension. There are some suggestions about the filename in the note below:

.. note::

    * **AB1** file should be indexed for better consistency with file naming regulation for *SangerContig* and *SangerAlignment*.
    * Forward or reverse direction should be specified in the filename.

:ref:`Figure_2<SangerRead_file_structure>` shows the suggested file naming strategy. The filename should contain four main parts: **"Contig name"**, **"Index number"**, **"Direction"** and **"ab1 file extension"**.

* **"Contig name"** :  :code:`RBNII397-13[C_LepFolF,C_lepFolR]`
* **"Index number"** :  :code:`1`
* **"Direction"** :  :code:`F`
* **"ab1 file extension"** :  :code:`.ab1`

.. _SangerRead_file_structure:
.. figure::  ../image/SangerRead_file_structure.png
   :align:   center
   :scale:   80 %

   Figure 2. *SangerRead* filename regulation.

In *SangerRead* section, it is not compulsory to follow the file naming regulation because users can directly specify the filename in input (see :ref:`Creating *SangerRead* instance from **AB1**`); however, in the *SangerContig* and *SangerAlignment*, sangeranalyseR will automatically group files, so it is compulsory to have systematic file naming strategy. For more details, please read :ref:`Advanced User Guide - *SangerContig* (**AB1**)` and :ref:`Advanced User Guide - *SangerAlignment* (**AB1**)`. :ref:`Figure_3<sangeranalyseR_filename_convention_SangerRead>` shows the suggested **AB1** file naming regulation.


.. _sangeranalyseR_filename_convention_SangerRead:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file naming regulation - *SangerRead*.


|

Creating *SangerRead* instance from **AB1**
-------------------------------------------

After preparing the *SangerRead* input **AB1** file, the next step is to create the *SangerRead* S4 instance by running :code:`SangerRead` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method which makes instance creation more intuitive. The inputs include **Basic Parameters**, **Trimming Parameters** and **Chromatogram Parameters** and most of them have their own default values. In the constructor below, we list important parameters.

.. code-block:: R

   sangerReadF <- SangerRead(inputSource           = "ABIF",
                             readFeature           = "Forward Read",
                             readFileName          = "RBNII397-13[C_LepFolF,C_LepFolR]_1_F.ab1",
                             geneticCode           = GENETIC_CODE,
                             TrimmingMethod        = "M1",
                             M1TrimmingCutoff      = 0.0001,
                             M2CutoffQualityScore  = NULL,
                             M2SlidingWindowSize   = NULL,
                             baseNumPerRow         = 100,
                             heightPerRow          = 200,
                             signalRatioCutoff     = 0.33,
                             showTrimmed           = TRUE)


The inputs of :code:`SangerRead` constructor function and :code:`new` method are same. For more details about *SangerRead* inputs & slots definition, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_. The created *SangerRead* instance, :code:`sangerRead`, is used as the input for the following functions.

|

Visualizing *SangerRead* trimmed read
-------------------------------------
Before going to :ref:`Writing *SangerRead* FASTA files :sub:\`(AB1)\`` and :ref:`Generating *SangerRead* report :sub:\`(AB1)\`` pages, it is suggested to visualize the trimmed *SangerRead*. Run the :code:`qualityBasePlot` function to get the result in :ref:`Figure_4 <SangerRead_qualityBasePlot>`. It shows the quality score for each base pairs and the trimming start/end points of the sequence.


.. _SangerRead_qualityBasePlot:
.. figure::  ../image/SangerRead_qualityBasePlot.png
   :align:   center
   :scale:   30 %

   Figure 4. *SangerRead* trimmed read visualization.

.. code-block:: R

   qualityBasePlot(sangerReadF)

|

Updating *SangerRead* quality trimming parameters
-------------------------------------------------
In the previous :ref:`Creating *SangerRead* instance from **AB1**` part, the constructor function applies the quality trimming parameters to the read. After creating the *SangerRead* S4 instance, users can change the trimming parameters by running :code:`updateQualityParam` function which will change the *QualityReport* instance inside the *SangerRead* and update frameshift amino acid sequences.

.. code-block:: R

   newSangerRead <- updateQualityParam(sangerReadF,
                                       TrimmingMethod       = "M2",
                                       M1TrimmingCutoff     = NULL,
                                       M2CutoffQualityScore = 29,
                                       M2SlidingWindowSize  = 15)

|



Writing *SangerRead* FASTA files :sub:`(AB1)`
-------------------------------------------------
Users can write the *SangerRead* instance to **FASTA** files. The trimmed read sequence will be written into a **FASTA** file. Below is the one-line function that users need to run. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFastaSA` function.

.. code-block:: R

   writeFastaSR(newSangerRead,
                outputDir         = tempdir(),
                compress          = FALSE,
                compression_level = NA)

Users can download the `output **FASTA** file <https://howardchao.github.io/sangeranalyseR_report/SangerRead/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F.fa>`_ of this example.

|


Generating *SangerRead* report :sub:`(AB1)`
------------------------------------------------
Last but not least, users can save *SangerRead* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files. The results in the report are static.

.. code-block:: R

   generateReportSR(newSangerRead,
                    outputDir = tempdir())

`SangerRead_Report_ab1.html <https://howardchao.github.io/sangeranalyseR_report/SangerRead/AB1/ACHLO006-09[LCO1490_t1,HCO2198_t1]_1_F/SangerRead_Report_ab1.html>`_ is the generated report html of this example. Users can access to '*Basic Information*', '*DNA Sequence*', '*Amino Acids Sequence*', '*Quality Trimming*' and '*Chromatogram*' sections inside this report.
