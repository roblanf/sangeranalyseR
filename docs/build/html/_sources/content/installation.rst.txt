Installation
===============

System requirements
-------------------
* R >= 4.0.0 (current)
* `Rstudio (recommended) <https://rstudio.com>`_

|

Install from Bioconductor
-------------------------

sangeranalyseR is on `Bioconductor 3.12 development <https://bioconductor.org/packages/devel/bioc/html/sangeranalyseR.html>`_ now.

.. _sangeranalyseR_bioconductor:
.. figure::  ../image/bioconductor.png
   :align:   center

   Figure 1. sangeranalyseR on Bioconductor 3.12 development.

To install this package, start R (version "4.0") and enter:

.. code-block:: R

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    # The following initializes usage of Bioc devel
    BiocManager::install(version='devel')

    BiocManager::install("sangeranalyseR")

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

After installing :code:`sangeranalyseR`, load it in R console.

.. code-block:: R

   library(sangeranalyseR)

Now, you are ready to go !

|


Where to go from here ?
-----------------------
Please continue to the :ref:`Quick Start Guide` or the more detailed :ref:`Beginners Guide`.
