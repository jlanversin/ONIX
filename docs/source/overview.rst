.. _overview:

================
Overview of ONIX
================

ONIX (for OpeN IsotopiX) is a state-of-the-art nuclear depletion software that is open-source. It can be used to model nuclear reactors simulation, estimate the production of fissile materials in reactors, support work in nuclear archaeology and for other applications. Because ONIX is open-source, anyone can use it for any type of applications (no export control or proprietary restrictions). This makes ONIX an ideal tools for innovative, community-based research on nuclear reactors and nuclear arms control.

The reference scientific article for ONIX and its benchmark validation can be found `here <https://www.sciencedirect.com/science/article/pii/S0306454920306009>`_.

-----------------------
Development and purpose
-----------------------

ONIX has been developed by `Julien de Troullioud de Lanversin <https://cisac.fsi.stanford.edu/people/julien-de-troullioud-de-lanversin>`_, `Moritz Kutt <https://ifsh.de/en/staff/kuett>`_ and `Alexander Glaser <https://sgs.princeton.edu/team/alex-glaser>`_.

The main impetus behind the development of ONIX was the lack of open-source simulation tools for research in nuclear reactors performance and nuclear arms control, especially for modelling isotopic depletion in nuclear systems. The vast majority of nuclear reactor physics codes are either proprietary or export-controlled which considerably restricts their distribution and their use. Unfortunately, these restrctions prevent community-based research and severily complicates collaboration between groups from different countries.

An important breakthrough was achieved with the development of `OpenMC <https://docs.openmc.org/en/stable/>`_, the first open-source Monte Carlo code for neutron transport. The original goal of ONIX was to provide an open-source software for isotopic depletion calculation that could be coupled with OpenMC in order to offer the first fully open-source reactor physics package (neutron transport + depletion calculation).

ONIX is also the first software designed to support research in nuclear archaelogy. Nuclear archaeology encompasses scientific methods and technics aimed at estimating the production of weapon-usable material in nuclear reactors (see `Fetter  <https://www.tandfonline.com/doi/abs/10.1080/08929889308426386>`_). Nuclear archaeology brings crucial technical contributions to nuclear arms control and disarmement efforts.

ONIX is written in Python 3 with an object-oriented design. Its Python API (Application Programming Interface) allows great flexibility, reusability, and readability for users and developers. As of now, the code is still at an early stage of development and many more features will be added. However, the software is fully functional and has been validated against established numerical and experimentl benchmarks. Our hope is that users and developers from all around the world will be interested in using ONIX and contributes to its improvement.

-------------------
Key functionalities
-------------------

ONIX is equipped with advanced algorithms that can undertake traditional tasks for a depletion code such as building the burnup matrix, solving the depletion equation, or computing the power of specific regions. The software includes a coupling interface to allow coupled simulations with OpenMC but it can also be used in a stand-alone mode. Thanks to its Python API, ONIX enables of great variety of output data processing.

Depletion capability
--------------------

The core function of ONIX is to deplete a system of nuclides exposed to a neutron flux or to compute their decays over time (when there are no neutron flux). To that end, ONIX uses its depletion solver module: Salamèche.

.. _fig-salameche:

.. figure:: ../_image/salameche.png
   :align: center
   :figclass: align-center

   Salamèche is used by ONIX to "burn" the system of isotopes 

Salamèche first builds the burnup matrix from the data provided by the nuclear data libraries. It then uses a Chebyshev Rational Approximation Method (`CRAM <https://www.tandfonline.com/doi/abs/10.13182/NSE15-26>`_) to solve the depletion equation and produce new isotopic densities.

Data libraries
--------------

Isomeric branching
------------------

Coupling with OpenMC
--------------------

Stand alone mode
----------------

Output data processing
----------------------

Scheduled changes in operation
------------------------------

--------------------------------
Nuclear archaeology capabilities
--------------------------------

Nuclear archaeology
-------------------

The NAX module in ONIX
----------------------

Example on a LWR reactor
------------------------

------------
Applications
------------

Nuclear reactor simulation
--------------------------

Fissile material production
---------------------------

Nuclear archaeology
-------------------

Nuclear fallouts
----------------

----------
Benchmarks
----------






