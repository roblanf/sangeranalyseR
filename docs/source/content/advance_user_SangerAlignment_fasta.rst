Advanced User Guide - *SangerAlignment* (**FASTA**)
===================================================

*SangerAlignment* is the highest class level in sangeranalyseR showed in :ref:`Figure_1<SangerAlignment_hierachy_fasta>`. It contains *SangerContig* list and the contigs alignment result. Users can access to *SangerContig* and *SangerRead* instance inside *SangerAlignment* instance. In this section, we are going to go through detailed sangeranalyseR data analysis steps in *SangerAlignment* level.

.. _SangerAlignment_hierachy_fasta:
.. figure::  ../image/SangerAlignment_hierachy.png
   :align:   center
   :scale:   20 %

   Figure 1. Classes hierarchy in sangeranalyseR, *SangerAlignment* level.

|

Preparing *SangerAlignment* **FASTA** input
-------------------------------------------
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


.. _SangerAlignment_file_structure_complex_fasta:
.. figure::  ../image/SangerAlignment_file_structure_complex.png
   :align:   center
   :scale:   120 %

   Figure 2. Input ab1 files inside the parent directory, :code:`./tmp/`.


:ref:`Figure_2<SangerAlignment_file_structure_complex_fasta>` shows the file naming regulation and directory hierarchy. In this example, the parent directory is :code:`extdata` and the directories in first layer are :code:`Allolobophora_chlorotica` and :code:`Drosophila_melanogaster`. All target **AB1** files need to be inside parent directory but it is not necessary to put them in the same level of directory. sangeranalyseR will recursively search all files with **.ab1** file extension and automatically group reads with the same contig name. The direction of reads in each contig will be grouped by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp` with filenames. Therefore, it is important to carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, :code:`"_[0-9]*_F.ab1$"`, :code:`"_[0-9]*_R.ab1$"` for matching forward and reverse reads are highly suggested and are used as default. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerAlignment_fasta:
.. figure::  ../image/sangeranalyseR_filename_convention.png
   :align:   center
   :scale:   25 %

   Figure 3. Suggested **AB1** file naming regulation - *SangerAlignment*.

:ref:`Figure_3<sangeranalyseR_filename_convention_SangerAlignment_fasta>` shows the suggested **AB1** file naming regulation. Users are strongly recommended to follow this file naming regulation and use the default :code:`suffixForwardRegExp` : ":code:`_[0-9]*_F.ab1$`" and :code:`suffixReverseRegExp` : ":code:`_[0-9]*_R.ab1$`" to reduce any chance of error.

|

Creating *SangerAlignment* instance from **FASTA**
--------------------------------------------------
After preparing the input directory, we can create the *SangerAlignment* S4 instance by running :code:`SangerAlignment` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list important parameters. For a simpler command, please go to :ref:`Quick Commands`.

.. code-block:: R

   sangerAlignment <- SangerAlignment(inputSource           = "FASTA",
                                      fastaFileName         = fastaFN,
                                      namesConversionCSV    = namesConversionCSV,
                                      suffixForwardRegExp   = suffixForwardRegExpFa,
                                      suffixReverseRegExp   = suffixReverseRegExpFa,
                                      refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN"
                                      )


The inputs of :code:`SangerAlignment` constructor function and :code:`new` method are same. For more details about *SangerAlignment* inputs & slots definition, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_.

|

Writing *SangerAlignment* FASTA files :sub:`(FASTA)`
----------------------------------------------------
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

Generating *SangerAlignment* report :sub:`(FASTA)`
--------------------------------------------------
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
