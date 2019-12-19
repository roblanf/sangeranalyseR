Getting Start
=============

System prerequisite
-------------------
* R >= 3.6.0 (current)
* Rstudio (recommend)

|

Installation
------------

Install the development version
+++++++++++++++++++++++++++++++
If you haven't install :code:`devtools` package before, please install it first.

.. code-block:: R

   install.packages("devtools")


Run the following code in your R console to install the newest version from Github.


.. code-block:: R

   library(devtools)
   install_github("roblanf/sangeranalyseR")
   library(sangeranalyseR)



Install from Bioconductor
+++++++++++++++++++++++++
After uploading to bioconductor !!!!

.. code-block:: R

   if (!requireNamespace("BiocManager"))
      install.packages("BiocManager")
      BiocManager::install()
   BiocManager::install("sangeranalyseR")


After installing :code:`sangeranalyseR`, load it in R console.

.. code-block:: R

   library(sangeranalyseR)

Now, you are ready to go !

|

Where to go from here ?
-----------------------
Please continue to the :ref:`Beginner Guide` for further usages.
