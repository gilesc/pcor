====
pcor
====

Simple C++ CLI program to calculate Pearson correlations for large matrices in
parallel using OpenMP. 

Installation
============

Requires armadillo library (Ubuntu/Debian package: "armadillo") and a
relatively recent C++ compiler (C++-0x support).

For local user install:

.. code-block:: bash

    make install

For sudo install:

.. code-block:: bash

    sudo make PREFIX=/usr install

Usage
=====

Input: a tab-delimited, labeled matrix on stdin. Example:

::

    	A	B	C
    X	1	2	3
    Y	4	5	6

There are two possible output modes, both on stdout, depending on whether the
``-n`` option is specified. If not, it outputs a full, all-vs-all columns
correlation matrix. If so, it outputs a tab-delimited format with the first
column corresponding to a column, and the subsequent columns corresponding, in
order, to the N most highly correlated input columns.

Use the program's help (``-h`` option) for further usage information.

Gotchas
=======

The computation is in parallel, so output order will not necessarily correspond
to the order of columns in the input matrix.

Max memory usage is approximately 2 * sizeof(double) * N, where N is the number
of elements in the input matrix. For most of the computation, it will use half
this amount.
