.. turbosoap_dscribe documentation master file, created by
   sphinx-quickstart on Wed Apr 29 16:38:44 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TurboSOAP_DScribe
=================

TurboSOAP-DScribe is a python package which provides TurboSOAP machine learning 
descriptor functionality to DScribe python package 
`https://singroup.github.io/dscribe/ <https://singroup.github.io/dscribe/>`_. 
Requires DScribe, Fortran compiler and LAPACK libraries.


Installation
------------
Easiest way to install is to install via pip.

.. code-block:: sh

    pip install turbosoap_dscribe

To install the latest development version, clone the source code from
github and install with pip from local file:

.. code-block:: sh

    git clone https://github.com/mcaroba/turbosoap_dscribe
    cd turbosoap_dscribe
    pip install .

If installation fails, please make sure that you have the prerequisites. Especially Fortran compiler and LAPACK libraries.
The setup.py script uses numpy.distutils to find Fortran compiler and LAPACK lbiraries.

TurboSOAP - a fast implementation of SOAP
=========================================
TurboSOAP is an implementation of the SOAP scheme which enables more control 
over the various aspects of the descriptor

Setup
-----
Instantiating the TurboSOAP descriptor is done as follows:

.. literalinclude:: ../turbosoap.py
    :language: python
    :lines: 1-3

The constructor takes the following parameters:

.. automethod:: dscribe.descriptors.turbosoap.TurboSOAP.__init__

After creation the TurboSOAP descriptor can applied to atomic structures:

.. literalinclude:: ../turbosoap.py
    :language: python
    :lines: 7-14

The TurboSOAP descriptor uses a lot of default parameters which are decent defaults but might not be 
physically or computationally optimal paremeters to represent atomic environments in specific cases. 
The parameters enable good accuracy and quick computation times if they are chosen judiciously.

In order to tune the parameters we need to use TurboSOAPSpecie to define them.

TurboSOAPSpecie
---------------
TurboSOAPSpecie can be used to define all the per-species parameters used by TurboSOAP.
A more involved example using TurboSOAPSpecies:

.. literalinclude:: ../turbosoap.py
    :language: python
    :lines: 16-24

.. automethod:: dscribe.descriptors.turbosoap.TurboSOAPSpecie.__init__

Comparing environments
----------------------
A quick example on how to compare environments. The quantity of interest is the similarity kernel between
two atomic environments which is between zero and one. Zero indicates that environments are not similar and one indicates that 
the environments are the same. The kernel is computed as a dot product between two TurboSOAP vectors:

.. literalinclude:: ../turbosoap.py
    :language: python
    :lines: 26-42

References
==========
If you use TurboSOAP, read and cite:
"Optimizing many-body atomic descriptors for enhanced computational 
performance of machine learning based interatomic potentials", 
Miguel Caro, Phys. Rev. B 100, 024112 (2019),
`https://doi.org/10.1103/PhysRevB.100.024112 <https://doi.org/10.1103/PhysRevB.100.024112>`_

In addition, read and cite the original SOAP paper:
"On representing chemical environments", Albert P. Bartók, Risi Kondor, and
Gábor Csányi, Phys. Rev. B 87, 184115, (2013),
`https://doi.org/10.1103/PhysRevB.87.184115 <https://doi.org/10.1103/PhysRevB.87.184115>`_.
