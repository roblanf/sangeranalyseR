Advanced User Guide - *SangerContig* (**FASTA**)
================================================

*SangerContig* is the second level in sangeranalyseR showed in :ref:`Figure_1<SangerContig_hierarchy_fasta>` which corresponds to a contig in Sanger sequencing. Among slots inside it, there are two lists, forward and reverse read list, storing *SangerRead* in the corresponding direction. In this section, we are going to go through details about sangeranalyseR data analysis in *SangerContig* level with **FASTA** file input.

.. _SangerContig_hierarchy_fasta:
.. figure::  ../image/SangerContig_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR, *SangerContig* level.

|

Preparing *SangerContig* **FASTA** input
----------------------------------------

We design the **FASTA** file input for those who do not want to do quality trimming and base calling for each *SangerRead* in *SangerContig*; therefore, it does not contain quality trimming and chromatogram input parameters and results in slots of *SangerRead*.


Before starting the analysis, users need to prepare one **FASTA** file containing sequence of all reads. Inside the **FASTA** file, there are strings starting with ">" before each read which are the read names. Because sangeranalyseR will group reads into "Forward Read List" and "Reverse Read List", users have to follow the naming regulations for the read names. Below are some regulations:

.. note::

    *  The same contig name must be included in all read names.
    *  Forward or reverse direction also has to be specified in the read names.

There are four parameters, :code:`fastaFileName`, :code:`contigName`, :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`,that users need to provide so that program can automatically match correct reads in **FASTA** file and divide them into forward and reverse direction.

.. note::

  * :code:`fastaFileName`: The **FASTA** file that contains sequence of all reads. The read names have to follow the naming regulation.
  * :code:`contigName`: The value of this parameter is a regular expression that matches read names that are going to be included in the *SangerContig* level analysis. :code:`grepl` function in R is used.
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all read names in forward direction. :code:`grepl` function in R is used to select forward reads from all **AB1** files.
  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all read names in reverse direction. :code:`grepl` function in R is used to select reverse reads from all **AB1** files.

namesConversionCSV

Here, we have an example:

.. _SangerContig_fasta_input:
.. figure::  ../image/SangerContig_fasta_input.png
   :align:   center
   :scale:   50 %

   Figure 2. *SangerContig* **FASTA** input file.

:ref:`Figure_2<SangerContig_fasta_input>` shows the **FASTA** input file and read names naming convention. sangeranalyseR will first match the :code:`contigName` to exclude unrelated reads and then separate the forward and reverse reads by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. Therefore, it is important to make sure all target reads share the same :code:`contigName` and carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1$"`, :code:`"_[0-9]*_R.ab1$"` for matching forward and reverse reads are highly suggested and are used as default. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerContig_fasta:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested read naming regulation - *SangerContig*.

:ref:`Figure_3<sangeranalyseR_filename_convention_SangerContig_fasta>` shows the suggested read naming convention. Users are strongly recommended to follow this file naming regulation and use the default :code:`suffixForwardRegExp` : ":code:`_[0-9]*_F.ab1$`" and :code:`suffixReverseRegExp` : ":code:`_[0-9]*_R.ab1$`" to reduce any chance of error.


|


Creating *SangerContig* instance from **FASTA**
-----------------------------------------------

After preparing the input directory, we can create the *SangerContig* S4 instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list important parameters.

.. code-block:: R

   sangerContigFa <- SangerContig(inputSource           = "FASTA",
                                  fastaFileName         = fastaFN,
                                  namesConversionCSV    = namesConversionCSV,
                                  contigName            = contigName,
                                  suffixForwardRegExp   = suffixForwardRegExpFa,
                                  suffixReverseRegExp   = suffixReverseRegExpFa,
                                  refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


In this example, :code:`contigName` is set to :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]"`, so only :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_1_F.ab1"` and :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_2_R.ab1"` reads will be selected from **FASTA** file to align to a contig.

The inputs of :code:`SangerContig` constructor function and :code:`new` method are same. For more details about *SangerContig* inputs & slots definition, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_.

|


Writing *SangerContig* FASTA files :sub:`(FASTA)`
-------------------------------------------------
Users can write the *SangerContig* instance to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`reads_unalignment`: Writing reads into a single **FASTA** file (only trimmed without alignment).
* :code:`reads_alignment`: Writing reads alignment and contig read to a single **FASTA** file.
* :code:`contig`: Writing the contig to a single **FASTA** file.
* :code:`all`: Writing reads, reads alignment, and the contig into three different files.

Below is the one-line function that users need to run. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFastaSA` function.

.. code-block:: R

   writeFastaSC(sangerContig,
                outputDir         = tempdir(),
                compress          = FALSE,
                compression_level = NA,
                selection         = "all")

|

Generating *SangerContig* report :sub:`(FASTA)`
-----------------------------------------------
Last but not least, users can save *SangerContig* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.


Users can set :code:`includeSangerRead` parameter to decide to which level the *SangerContig* report will go. Moreover, after the reports are generated, users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it would take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReportSC(sangerContig,
                    outputDir           = tempdir(),
                    includeSangerRead   = TRUE)
