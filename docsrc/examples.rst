TurboSOAP Example
=================

Setup
-----
Instantiating the TurboSOAP descriptor is done as follows:

.. literalinclude:: ../better_example.py
    :language: python
    :lines: 1-3

After creation the TurboSOAP descriptor can applied to atomic structures:

.. literalinclude:: ../better_example.py
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

.. literalinclude:: ../better_example.py
    :language: python
    :lines: 16-24

.. automethod:: turbosoap_dscribe.TurboSOAPSpecie.__init__


Comparing environments
----------------------
A quick example on how to compare environments. The quantity of interest is the similarity kernel between
two atomic environments which is between zero and one. Zero indicates that environments are not similar and one indicates that 
the environments are the same. The kernel is computed as a dot product between two TurboSOAP vectors:

.. literalinclude:: ../better_example.py
    :language: python
    :lines: 26-42