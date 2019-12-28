Quick Commands
==============
This page is the concise and core analysis steps of sangeranalyseR. For more details, please refer to :ref:`Beginner Guide`.

* **Step 1** : Input files preparation

   Put all your **AB1** files in a directory :code:`./tmp/`.

* **Step 2** : *SangerAlignment* creation

.. code-block:: R

   sangerAlignment <- SangerAlignment(parentDirectory     = "./tmp/",
                                      suffixForwardRegExp = "_F.ab1",
                                      suffixReverseRegExp = "_R.ab1")

* **Step 3** : Launching Shiny App

.. code-block:: R

   launchAppSA(sangerAlignment)


* **Step 4** : Writing FASTA file

.. code-block:: R

   writeFastaSA(sangerAlignment)


* **Step 5** : Generating report

.. code-block:: R

   generateReportSA(sangerAlignment)

|

For more detailed analysis steps, please choose one the following topics :

* :ref:`Beginner Guide`

* :ref:`Advanced User Guide - *SangerRead* (**AB1**)`

* :ref:`Advanced User Guide - *SangerContig* (**AB1**)`

* :ref:`Advanced User Guide - *SangerAlignment* (**AB1**)`

* :ref:`Advanced User Guide - *SangerRead* (**fasta**)`

* :ref:`Advanced User Guide - *SangerContig* (**fasta**)`

* :ref:`Advanced User Guide - *SangerAlignment* (**fasta**)`
