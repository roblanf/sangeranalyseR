Advanced User Guide - *SangerRead*
==================================

*SangerRead* is the highest level in sangeranalyseR. :ref:`Figure_1<SangerRead_hierarchy>`

.. _SangerRead_hierarchy:
.. figure::  ../image/SangerRead_hierarchy.png
   :align:   center
   :scale:   20 %

   Figure 1. Hierarchy of classes in sangeranalyseR.


   *SangerAlignment* creation
   --------------------------
   Now we can create the *SangerAlignment* S4 instance by :code:`SangerAlignment` constructor function or :code:`new` method. The constructor function is a wrapper for :code:`new` method.

   .. code-block:: R

      sangerAlignment <- SangerAlignment(inputSource          = "ABIF",
                                         parentDirectory      = "./tmp/",
                                         suffixForwardRegExp  = "_[0-9]*_F.ab1$",
                                         suffixReverseRegExp  = "_[0-9]*_R.ab1$",
                                         refAminoAcidSeq = "SRQWLFSTNHKDIGTLYFIFGAWAGMVGTSLSILIRAELGHPGALIGDDQIYNVIVTAHAFIMIFFMVMPIMIGGFGNWLVPLMLGAPDMAFPRMNNMSFWLLPPALSLLLVSSMVENGAGTGWTVYPPLSAGIAHGGASVDLAIFSLHLAGISSILGAVNFITTVINMRSTGISLDRMPLFVWSVVITALLLLLSLPVLAGAITMLLTDRNLNTSFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGMISHIISQESGKKETFGSLGMIYAMLAIGLLGFIVWAHHMFTVGMDVDTRAYFTSATMIIAVPTGIKIFSWLATLHGTQLSYSPAILWALGFVFLFTVGGLTGVVLANSSVDIILHDTYYVVAHFHYVLSMGAVFAIMAGFIHWYPLFTGLTLNNKWLKSHFIIMFIGVNLTFFPQHFLGLAGMPRRYSDYPDAYTTWNIVSTIGSTISLLGILFFFFIIWESLVSQRQVIYPIQLNSSIEWYQNTPPAEHSYSELPLLTN",
                                         TrimmingMethod        = "M1",
                                         M1TrimmingCutoff      = 0.0001,
                                         M2CutoffQualityScore  = NULL,
                                         M2SlidingWindowSize   = NULL,
                                         baseNumPerRow         = 100,
                                         heightPerRow          = 200,
                                         signalRatioCutoff     = 0.33,
                                         showTrimmed           = TRUE)

   Inputs Definition
   +++++++++++++++++
   The inputs of :code:`SangerAlignment` constructor function and :code:`new` method are same. Both of their input parameters can be divided into three main categories: :code:`Basic Parameters`, :code:`Quality Trimming Parameters`, and :code:`Chromatogram Parameters`. Here, we will only explain important parameters. For more details, please refer to `sangeranalyseR reference manul <http://packages.python.org/an_example_pypi_project/>`_.

   * Basic Parameters:

      * :code:`readFeature`: Specify whether the read is forward or reverse reads. The value must be “Forward Read” or “Reverse Read”.
      * :code:`parentDirectory`: see :ref:`Input files preparation`.
      * :code:`suffixForwardRegExp`: see :ref:`Input files preparation`.
      * :code:`suffixReverseRegExp`: see :ref:`Input files preparation`.
      * :code:`refAminoAcidSeq`: a reference amino acid sequence for the primary sequence as an AAString object in Biostring R packge.


   * Quality Trimming Parameters:

      * :code:`TrimmingMethod`: Specify the read trimming method for this SangerRead. The value must be “M1” (the default) or ‘M2’.
      * :code:`M1TrimmingCutoff`: The trimming cutoff for the Method 1. If TrimmingMethod is “M1”, then the default value is “0.0001”. Otherwise, the value must be “NULL”
      * :code:`M2CutoffQualityScore`: The trimming cutoff quality score for the Method 2. If TrimmingMethod is “M2”, then the default value is “20”. Otherwise, the value must be “NULL”. It works with M2SlidingWindowSize.
      * :code:`M2SlidingWindowSize`: The trimming sliding window size for the Method 2. If TrimmingMethod is “M2”, then the default value is “5”. Otherwise, the value must be “NULL”. It works with M2CutoffQualityScore.

   * Chromatogram Parameters:

      * :code:`baseNumPerRow`: It defines maximum base pairs in each row. The default value is “100”.
      * :code:`heightPerRow`: It defines the height of each row in chromatogram. The default value is “200”.
      * :code:`signalRatioCutoff`: The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is “0.33”.
      * :code:`showTrimmed`: The logical value storing whether to show trimmed base pairs in chromatogram. The default value is “TRUE”.

   Slots Definition
   ++++++++++++++++
