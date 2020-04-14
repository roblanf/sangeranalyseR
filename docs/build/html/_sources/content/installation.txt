Installation
===============

System requirements
-------------------
* R >= 3.6.0 (current)
* `Rstudio (recommended) <https://rstudio.com>`_

|


Install the development version
-------------------------------

If you haven't installed the :code:`devtools` package before, please install it first:

.. code-block:: R

   install.packages("devtools")


Then run the following code in your R console to install the newest version from Github.


.. code-block:: R

   library(devtools)
   install_github("roblanf/sangeranalyseR", ref = "develop")
   library(sangeranalyseR)

|

Install from Bioconductor
-------------------------
NB: This is currently a placeholder - the package isn't on Bioconductor yet...
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
Please continue to the :ref:`Quick Start Guide` or the more detailed :ref:`Beginners Guide`.
