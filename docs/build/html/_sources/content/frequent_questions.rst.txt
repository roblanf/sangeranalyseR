Frequently Asked Questions
==========================

|

Q: What is the difference between two different trimming methods?
-----------------------------------------------------------------
A: In sangeranalyseR, we provide two trimming methods, *"M1"* (the default) or *"M2"*, which represents *"method 1"* or *"method 2"* respectively. M1 is the modified Mott's trimming algorithm that can also be found in Phred/Phrap and Biopython. M2 is like trimmomatic's sliding window method. If you want to set M1 as your trimming method, you need to assign **"TrimmingMethod"** to **"M1"** and **"M1TrimmingCutoff"** as the threshold that you want. Its default value is *"0.0001"*. In contrast, you can assign **"TrimmingMethod"** to **"M2"** if you want to set M2 as your trimming method. **"M2CutoffQualityScore"** and **"M2SlidingWindowSize"** are two parameters that control M2 trimming and their default values are *"20"* and *"10"* respectively.
