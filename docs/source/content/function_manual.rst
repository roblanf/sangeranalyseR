User Manual (functions)
=======================

This is the link to `sangeranalyseR user manual  <https://bioconductor.org/packages/devel/bioc/manuals/sangeranalyseR/man/sangeranalyseR.pdf>`_. Following are input parameters for **SangerRead**, **SangerContig**, and **SangerAlignment** constructors.

|

SangerRead Constructor Parameters
---------------------------------

.. code-block:: R

   SangerRead(inputSource = "ABIF",
              readFeature = "",
              readFileName = "",
              fastaReadName = "",
              geneticCode = GENETIC_CODE,
              TrimmingMethod = "M1",
              M1TrimmingCutoff = 0.0001,
              M2CutoffQualityScore = NULL,
              M2SlidingWindowSize = NULL,
              baseNumPerRow = 100,
              heightPerRow = 200,
              signalRatioCutoff = 0.33,
              showTrimmed = TRUE)

* **inputSource**: The input source of the raw file. It must be *"ABIF"* or *"FASTA"*. The default value is *"ABIF"*.
* **readFeature**: The direction of the Sanger read. The value must be *"Forward Read"* or *"Reverse Read"*.
* **readFileName**: The absolute filename of the target ABIF or FASTA file.
* **fastaReadName**:  If *"inputSource"* is *"FASTA"*, then this value has to be the name of the read inside the FASTA file; if *"inputSource"* is *"ABIF"*, then this value is *"NULL"* by default.
* **geneticCode**: Named character vector in the same format as *"GENETIC_CODE"* (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the *"getGeneticCode()"* function. The default is the standard code.
* **TrimmingMethod**: The read trimming method for the *SangerRead*. The value must be *"M1"* (the default) or *"M2"*, which represents *"method 1"* or *"method 2"* respectively. M1 is the modified Mott's trimming algorithm that can also be found in Phred/Phrap and Biopython. M2 is like trimmomatic's sliding window method.
* **M1TrimmingCutoff**: The cutoff for the trimming method 1. If *TrimmingMethod* is *"M1"*, then the default value is *"0.0001"*. Otherwise, the value must be *"NULL"*.
* **M2CutoffQualityScore**: The trimming cutoff quality score for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"20"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2SlidingWindowSize*.
* **M2SlidingWindowSize**: The trimming sliding window size for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"10"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2CutoffQualityScore*.
* **baseNumPerRow**: This parameter is related to chromatogram and defines maximum base pairs in each row. The default value is *"100"*.
* **heightPerRow**: This parameter is related to chromatogram and defines the height of each row in chromatogram. The default value is *"200"*.
* **signalRatioCutoff**: The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is *"0.33"*. This parameter is related to chromatogram.
* **showTrimmed**: The logical value storing whether to show trimmed base pairs in chromatogram. The default value is *"TRUE"*.

|

SangerContig Constructor Parameters
-----------------------------------

.. code-block:: R

   SangerContig(inputSource            = "ABIF",
                fastaFileName          = "",
                namesConversionCSV     = NULL,
                parentDirectory        = "",
                contigName             = "",
                suffixForwardRegExp    = "_F.ab1",
                suffixReverseRegExp    = "_R.ab1",
                TrimmingMethod         = "M1",
                M1TrimmingCutoff       = 0.0001,
                M2CutoffQualityScore   = NULL,
                M2SlidingWindowSize    = NULL,
                baseNumPerRow          = 100,
                heightPerRow           = 200,
                signalRatioCutoff      = 0.33,
                showTrimmed            = TRUE,
                refAminoAcidSeq        = "",
                minReadsNum            = 2,
                minReadLength          = 20,
                minFractionCall        = 0.5,
                maxFractionLost        = 0.5,
                geneticCode            = GENETIC_CODE,
                acceptStopCodons       = TRUE,
                readingFrame           = 1,
                processorsNum          = NULL)

* **inputSource**: The input source of the raw file. It must be *"ABIF"* or *"FASTA"*. The default value is *"ABIF"*.
* **fastaFileName**: If *"inputSource"* is *"FASTA"*, then this value has to be the name of the FASTA file; if *"inputSource"* is *"ABIF"*, then this value is *"NULL"* by default.
* **namesConversionCSV**: The absolute filename of CSV file that provides read names  following the naming regulation. If *"inputSource"* is *"FASTA"*, then users need to prepare the csv file or make sure the original names inside FASTA file are valid; if *"inputSource"* is *"ABIF"*, then this value is *"NULL"* by default.
* **parentDirectory**: The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerContig, all reads must be in the first layer in this directory.
* **contigName**: The contig name of all the reads in *"parentDirectory"*.
* **suffixForwardRegExp**: The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be "_F.ab1".
* **suffixReverseRegExp**: The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be "_R.ab1".
* **TrimmingMethod**: The read trimming method for the *SangerRead*. The value must be *"M1"* (the default) or *"M2"*, which represents *"method 1"* or *"method 2"* respectively. M1 is the modified Mott's trimming algorithm that can also be found in Phred/Phrap and Biopython. M2 is like trimmomatic's sliding window method.
* **M1TrimmingCutoff**: The cutoff for the trimming method 1. If *TrimmingMethod* is *"M1"*, then the default value is *"0.0001"*. Otherwise, the value must be *"NULL"*.
* **M2CutoffQualityScore**: The trimming cutoff quality score for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"20"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2SlidingWindowSize*.
* **M2SlidingWindowSize**: The trimming sliding window size for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"10"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2CutoffQualityScore*.
* **baseNumPerRow**: This parameter is related to chromatogram and defines maximum base pairs in each row. The default value is *"100"*.
* **heightPerRow**: This parameter is related to chromatogram and defines the height of each row in chromatogram. The default value is *"200"*.
* **signalRatioCutoff**: The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is *"0.33"*. This parameter is related to chromatogram.
* **showTrimmed**: The logical value storing whether to show trimmed base pairs in chromatogram. The default value is *"TRUE"*.
* **refAminoAcidSeq**: An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is "".
* **minReadsNum**: The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is *"2"*.
* **minReadLength**: Reads shorter than this will not be included in the readset. The default *"20"* means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
* **minFractionCall**: Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
* **maxFractionLost**: Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
* **geneticCode**: Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
* **acceptStopCodons**: The logical value \code{TRUE} or \code{FALSE}. \code{TRUE} (the defualt): keep all reads, regardless of whether they have stop codons; \code{FALSE}: reject reads with stop codons. If \code{FALSE} is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
* **readingFrame**: \code{1}, \code{2}, or \code{3}. Only used if \code{accept.stop.codons == FALSE}. This specifies the reading frame that is used to determine stop codons. If you use a \code{refAminoAcidSeq}, then the frame should always be \code{1}, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame.
* **processorsNum**: The number of processors to use, or NULL (the default) for all available processors.

|

SangerAlignment Constructor Parameters
--------------------------------------

.. code-block:: R

   SangerAlignment(inputSource            = "ABIF",
                   fastaFileName          = "",
                   namesConversionCSV     = NULL,
                   parentDirectory        = "",
                   suffixForwardRegExp    = "_F.ab1",
                   suffixReverseRegExp    = "_R.ab1",
                   TrimmingMethod         = "M1",
                   M1TrimmingCutoff       = 0.0001,
                   M2CutoffQualityScore   = NULL,
                   M2SlidingWindowSize    = NULL,
                   baseNumPerRow          = 100,
                   heightPerRow           = 200,
                   signalRatioCutoff      = 0.33,
                   showTrimmed            = TRUE,
                   refAminoAcidSeq        = "",
                   minReadsNum            = 2,
                   minReadLength          = 20,
                   minFractionCall        = 0.5,
                   maxFractionLost        = 0.5,
                   geneticCode            = GENETIC_CODE,
                   acceptStopCodons       = TRUE,
                   readingFrame           = 1,
                   minFractionCallSA      = 0.5,
                   maxFractionLostSA      = 0.5,
                   processorsNum          = NULL)

* **inputSource**: The input source of the raw file. It must be *"ABIF"* or *"FASTA"*. The default value is *"ABIF"*.
* **fastaFileName**: If *"inputSource"* is *"FASTA"*, then this value has to be the name of the FASTA file; if *"inputSource"* is *"ABIF"*, then this value is *"NULL"* by default.
* **namesConversionCSV**: The file path to the CSV file that provides read names that follow the naming regulation. If *"inputSource"* is *"FASTA"*, then users need to prepare the csv file or make sure the original names inside FASTA file are valid; if *"inputSource"* is *"ABIF"*, then this value is *"NULL"* by default.
* **parentDirectory**: The parent directory of all of the reads contained in ABIF format you wish to analyse. In SangerContig, all reads must be in the first layer in this directory.
* **suffixForwardRegExp**: The suffix of the filenames for forward reads in regular expression, i.e. reads that do not need to be reverse-complemented. For forward reads, it should be "_F.ab1".
* **suffixReverseRegExp**: The suffix of the filenames for reverse reads in regular expression, i.e. reads that need to be reverse-complemented. For revcerse reads, it should be "_R.ab1".
* **TrimmingMethod**: The read trimming method for the *SangerRead*. The value must be *"M1"* (the default) or *"M2"*, which represents *"method 1"* or *"method 2"* respectively. M1 is the modified Mott's trimming algorithm that can also be found in Phred/Phrap and Biopython. M2 is like trimmomatic's sliding window method.
* **M1TrimmingCutoff**: The cutoff for the trimming method 1. If *TrimmingMethod* is *"M1"*, then the default value is *"0.0001"*. Otherwise, the value must be *"NULL"*.
* **M2CutoffQualityScore**: The trimming cutoff quality score for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"20"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2SlidingWindowSize*.
* **M2SlidingWindowSize**: The trimming sliding window size for the trimming method 2. If *TrimmingMethod* is *"M2"*, then the default value is *"10"*. Otherwise, the value must be *"NULL"*. This parameter works with *M2CutoffQualityScore*.
* **baseNumPerRow**: This parameter is related to chromatogram and defines maximum base pairs in each row. The default value is *"100"*.
* **heightPerRow**: This parameter is related to chromatogram and defines the height of each row in chromatogram. The default value is *"200"*.
* **signalRatioCutoff**: The ratio of the height of a secondary peak to a primary peak. Secondary peaks higher than this ratio are annotated. Those below the ratio are excluded. The default value is *"0.33"*. This parameter is related to chromatogram.
* **showTrimmed**: The logical value storing whether to show trimmed base pairs in chromatogram. The default value is *"TRUE"*.
* **refAminoAcidSeq**: An amino acid reference sequence supplied as a string or an AAString object. If your sequences are protein-coding DNA seuqences, and you want to have frameshifts automatically detected and corrected, supply a reference amino acid sequence via this argument. If this argument is supplied, the sequences are then kept in frame for the alignment step. Fwd sequences are assumed to come from the sense (i.e. coding, or "+") strand. The default value is "".
* **minReadsNum**: The minimum number of reads required to make a consensus sequence, must be 2 or more. The default value is *"2"*.
* **minReadLength**: Reads shorter than this will not be included in the readset. The default *"20"* means that all reads with length of 20 or more will be included. Note that this is the length of a read after it has been trimmed.
* **minFractionCall**: Minimum fraction of the sequences required to call a consensus sequence for SangerContig at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
* **maxFractionLost**: Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerContig (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
* **geneticCode**: Named character vector in the same format as \code{GENETIC_CODE} (the default), which represents the standard genetic code. This is the code with which the function will attempt to translate your DNA sequences. You can get an appropriate vector with the getGeneticCode() function. The default is the standard code.
* **acceptStopCodons**: The logical value \code{TRUE} or \code{FALSE}. \code{TRUE} (the defualt): keep all reads, regardless of whether they have stop codons; \code{FALSE}: reject reads with stop codons. If \code{FALSE} is selected, then the number of stop codons is calculated after attempting to correct frameshift mutations (if applicable).
* **readingFrame**: \code{1}, \code{2}, or \code{3}. Only used if \code{accept.stop.codons == FALSE}. This specifies the reading frame that is used to determine stop codons. If you use a \code{refAminoAcidSeq}, then the frame should always be \code{1}, since all reads will be shifted to frame 1 during frameshift correction. Otherwise, you should select the appropriate reading frame.
* **minFractionCallSA**: Minimum fraction of the sequences required to call a consensus sequence for SangerAlignment at any given position (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.75 implying that 3/4 of all reads must be present in order to call a consensus.
* **maxFractionLostSA**: Numeric giving the maximum fraction of sequence information that can be lost in the consensus sequence for SangerAlignment (see the ConsensusSequence() function from DECIPHER for more information). Defaults to 0.5, implying that each consensus base can ignore at most 50 percent of the information at a given position.
* **processorsNum**: The number of processors to use, or NULL (the default) for all available processors.
