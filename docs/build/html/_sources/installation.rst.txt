.. _installation:

===================
Installation Guide
===================

As of now, ONIX can only be installed from source.

----------------------------------------
Installing from source for Linux and Mac
----------------------------------------

In order to use ONIX coupled with OpenMC, it is necessary to first install OpenMC on your computer. Please refer to `OpenMC documentation <https://docs.openmc.org/en/stable/index.html>`_ to install OpenMC 

There is only one python dependency that is probably not installed on your computer that ONIX needs, the module `Networkx <https://networkx.github.io/>`_

To install Networkx, type the following:

.. code-block:: sh

    pip install networkx

Now we need to clone ONIX repository. You should go directly to the GitHub repository to get the latest version of `ONIX <https://github.com/jlanversin/ONIX>`_

.. code-block:: sh

    git clone https://github.com/jlanversin/ONIX.git
    git checkout master

Finally, navigate into ONIX folder and install ONIX on your computer:

.. code-block:: sh

    cd ONIX
    pip install .