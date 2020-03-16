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

We design the **FASTA** file input for those who do not want to do quality trimming and base calling for each *SangerRead* in *SangerContig*; therefore, it does not contain quality trimming and chromatogram input parameters and results in *SangerRead* slots. Before starting the analysis, users need to prepare one **FASTA** file containing sequence of all reads. Inside the **FASTA** file, the strings starting with ">" before each read are the read names. Because sangeranalyseR will automatically group reads into "Forward Read List" and "Reverse Read List", users have to follow the naming regulations. Below are some regulations:

.. note::

    *  The same contig name must be included in all read names.
    *  Forward or reverse direction also has to be specified in the read names.

There are four parameters, :code:`fastaFileName`, :code:`contigName`, :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`,that users need to provide so that program can automatically match correct reads in **FASTA** file and divide them into forward and reverse direction.

.. note::

  * :code:`fastaFileName`: The path of **FASTA** file that contains sequence of all reads. The read names have to follow the naming regulation.
  * :code:`contigName`: The value of this parameter is a regular expression that matches read names that are going to be included in the *SangerContig* level analysis. :code:`grepl` function in R is used.
  * :code:`suffixForwardRegExp`: The value of this parameter is a regular expression that matches all read names in forward direction. :code:`grepl` function in R is used to select forward reads from all read names in **FASTA** files.
  * :code:`suffixReverseRegExp`: The value of this parameter is a regular expression that matches all read names in reverse direction. :code:`grepl` function in R is used to select reverse reads from all read names in **FASTA** files.

No doubt read names in the original **FASTA** file will not follow the naming regulation; however, it is highly not recommended to change the name directly in the raw **FASTA** file. Therefore, we provide a feature to let users do read names mapping conversion by a **CSV** file showed in :ref:`Figure_2<SangerContig_read_names_conversion>`. The first column is "original_read_name" which are the read names in the raw **FASTA** file, and the second column is "analysis_read_name" which are the read names that follow the naming regulation. The read names will be mapped onto the names in "original_read_name" without changing the raw **FASTA** file. :code:`namesConversionCSV` is the parameter that stores the path to this **CSV** file.

.. _SangerContig_read_names_conversion:
.. figure::  ../image/SangerContig_read_names_conversion.png
   :align:   center
   :scale:   40 %

   Figure 2. *SangerContig* **CSV** file - read names conversion.


Here, we have an example:

.. _SangerContig_fasta_input:
.. figure::  ../image/SangerContig_fasta_input.png
   :align:   center
   :scale:   40 %

   Figure 3. *SangerContig* **FASTA** input file.

:ref:`Figure_3<SangerContig_fasta_input>` shows the **FASTA** input file and the read names in it will be mapped onto the **CSV** file showed in :ref:`Figure_2<SangerContig_read_names_conversion>`. sangeranalyseR will first match the :code:`contigName` to exclude unrelated reads and then separate the forward and reverse reads by matching :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. Therefore, it is important to make sure all target reads share the same :code:`contigName` and carefully select :code:`suffixForwardRegExp` and :code:`suffixReverseRegExp`. The bad file naming and wrong regex matching might accidentally include reverse reads into the forward read list or vice versa, which will make the program generate totally wrong results. Therefore, users should have a consistent naming strategy. In this example, ":code:`_[0-9]*_F$`", ":code:`_[0-9]*_R$`" for matching forward and reverse reads are highly suggested. It is a good habit to index your reads in the same contig group because there might be more than one read that are in the forward or reverse direction.

.. _sangeranalyseR_filename_convention_SangerContig_fasta:
.. figure::  ../image/sangeranalyseR_filename_convention_fasta.png
   :align:   center
   :scale:   25 %

   Figure 4. Suggested read naming regulation in **FASTA** file - *SangerContig*.

:ref:`Figure_4<sangeranalyseR_filename_convention_SangerContig_fasta>` shows the suggested read naming regulation which is used in the "analysis_read_name" column in **CSV** file (:ref:`Figure_2<SangerContig_read_names_conversion>`). Users are strongly recommended to follow this file naming regulation and use :code:`suffixForwardRegExp` : ":code:`_[0-9]*_F$`" and :code:`suffixReverseRegExp` : ":code:`_[0-9]*_R$`" to reduce any chance of error.

|

Creating *SangerContig* instance from **FASTA**
-----------------------------------------------

After preparing the input directory, we can create the *SangerContig* S4 instance by running :code:`SangerContig` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method and it makes instance creation more intuitive. Most parameters in the constructor have their own default values. In the constructor below, we list important parameters.

.. code-block:: R

   sangerContigFa <- SangerContig(inputSource           = "FASTA",
                                  fastaFileName         = "ACHLO006-09[LCO1490_t1,HCO2198_t1].fa",
                                  namesConversionCSV    = "names_conversion_1.csv",
                                  contigName            = "ACHLO006-09[LCO1490_t1,HCO2198_t1]",
                                  suffixForwardRegExp   = "_[0-9]*_F$",
                                  suffixReverseRegExp   = "_[0-9]*_R$",
                                  refAminoAcidSeq       = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN")


In this example, :code:`contigName` is set to :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]"`, so only :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_1_F.ab1"` and :code:`"ACHLO006-09[LCO1490_t1.HCO2198_t1]_2_R.ab1"` reads will be selected from **FASTA** file to align into a contig.

The inputs of :code:`SangerContig` constructor function and :code:`new` method are same. For more details about *SangerContig* inputs and slots definition, please refer to `sangeranalyseR reference manual (need update) <http://packages.python.org/an_example_pypi_project/>`_.

|


Writing *SangerContig* FASTA files :sub:`(FASTA)`
-------------------------------------------------
Users can write the *SangerContig* instance to **FASTA** files. There are four options for users to choose from in :code:`selection` parameter.

* :code:`reads_unalignment`: Writing reads into a single **FASTA** file.
* :code:`reads_alignment`: Writing reads alignment and the aligned contig to a single **FASTA** file.
* :code:`contig`: Writing the contig to a single **FASTA** file.
* :code:`all`: Executing the three options mentioned above and writing *SangerContig* instance into three different files.

Below is the one-line function that users need to run. This function mainly depends on :code:`writeXStringSet` function in `Biostrings <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>`_ R package. Users can set the compression level through :code:`writeFastaSA` function.

.. code-block:: R

   writeFastaSC(sangerContigFa,
                outputDir         = tempdir(),
                compress          = FALSE,
                compression_level = NA,
                selection         = "all")


Users can download the output **FASTA** file of this example through the following three links:

* `reads_unalignment FASTA file <https://howardchao.github.io/sangeranalyseR_report/SangerContig/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_reads_unalignment.fa>`_
* `reads_alignment FASTA file <https://howardchao.github.io/sangeranalyseR_report/SangerContig/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_reads_alignment.fa>`_
* `contig FASTA file <https://howardchao.github.io/sangeranalyseR_report/SangerContig/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]_contig.fa>`_

|

Generating *SangerContig* report :sub:`(FASTA)`
-----------------------------------------------
Last but not least, users can save *SangerContig* instance into a report after the analysis. The report will be generated in **HTML** by knitting **Rmd** files.


Users can set :code:`includeSangerRead` parameter to decide to which level the *SangerContig* report will go. Moreover, after the reports are generated, users can easily navigate through reports in different levels within the **HTML** file.

One thing to pay attention to is that if users have many reads, it would take quite a long time to write out all reports. If users only want to generate the contig result, remember to set :code:`includeSangerRead` to :code:`FALSE` in order to save time.

.. code-block:: R

   generateReportSC(sangerContigFa,
                    outputDir           = tempdir(),
                    includeSangerRead   = TRUE)

Users can access to '*Basic Information*', '*SangerContig Input Parameters*', '*Contig Sequence*' and '*Contig Results*' sections inside the generated `SangerContig html report of this example <https://howardchao.github.io/sangeranalyseR_report/SangerContig/FASTA/ACHLO006-09[LCO1490_t1,HCO2198_t1]/SangerContig_Report.html>`_. Furthermore, users can also navigate through html reports of all forward and reverse *SangerRead* in this *SangerContig* report.
