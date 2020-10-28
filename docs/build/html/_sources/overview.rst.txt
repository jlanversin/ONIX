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

ONIX is also the first software designed to support research in nuclear archaelogy. Nuclear archaeology encompasses scientific methods and technics aimed at estimating the production of weapon-usable materials in nuclear reactors (see `Fetter  <https://www.tandfonline.com/doi/abs/10.1080/08929889308426386>`_). Nuclear archaeology brings crucial technical contributions to nuclear arms control and disarmement efforts.

ONIX is written in Python 3 with an object-oriented design. Its Python API (Application Programming Interface) allows great flexibility, reusability, and readability for users and developers. As of now, the code is still at an early stage of development and many more features will be added. However, the software is fully functional and has been validated against established numerical and experimentl benchmarks. Our hope is that users and developers from all around the world will be interested in using ONIX and contributes to its improvement.

.. fig-workshop:

.. figure:: _image/workshop.JPG
   :width: 600
   :align: center
   :figclass: align-center

   The very first workshop on ONIX, August 2019 in Aix-la-Chapelle, Germany.

-------------------
Key functionalities
-------------------

ONIX is equipped with advanced algorithms that can undertake traditional tasks for a depletion code such as building the burnup matrix, solving the depletion equation, or computing the power of specific regions. The software includes a coupling interface to allow coupled simulations with OpenMC but it can also be used in a stand-alone mode. Thanks to its Python API, ONIX enables of great variety of output data processing.

Depletion capability with Salamèche
-----------------------------------

The core function of ONIX is to deplete a system of nuclides exposed to a neutron flux or to compute their decays over time (when there are no neutron flux). To that end, ONIX uses its depletion solver module: Salamèche.

.. _fig-salameche:

.. figure:: _image/salameche.png
   :width: 200
   :align: center
   :figclass: align-center

   Salamèche is used by ONIX to "burn" the system of isotopes 

Salamèche first builds the burnup matrix from the data provided by the nuclear data libraries. It then uses a Chebyshev Rational Approximation Method (`CRAM <https://www.tandfonline.com/doi/abs/10.13182/NSE15-26>`_) to solve the depletion equation and produce new isotopic densities.

.. _fig-matrix:

.. figure:: _image/matrix.png
   :width: 600
   :align: center
   :figclass: align-center

   Burnup matrix build by Salamèche 

Data libraries
--------------

ONIX can read and process any nuclear data library as long as it is converted into the specific format processed by ONIX. The software repository currently provides two of the most recent nuclear data libraries for decay and fission yields data, the ENDF/B-VIII.0 and JEFF3.2 libraries, converted into ONIX's format and readily usable.

ONIX can model all the reactions important for nuclear reactor simulation: :math:`\beta^{+/-}`, :math:`\alpha`, and :math:`\gamma` decay and the :math:`(n,fission)`, :math:`(n,\gamma)`, :math:`(n,2n)`, :math:`(n,3n)`, :math:`(n,\alpha)`, and :math:`(n,p)` reactions. In addition, the :math:`(n,t)` reaction (:math:`t` for tritium) is also modelled to allow for estimation of tritium production in nuclear reactors.

While ONIX can build burnup matrix of arbitrary size from any libraries provided, large burnup matrices tend to make simulation significantly slower. ONIX uses a script that scans through the nuclides network, identifies elements that will never be produced (isolated nuclides), and remove them to produce reduced, optimized libraries. The ENDF/B-VIII.0 and JEFF3.2 libraries provided with ONIX have already been "reduced" via this method.

.. _fig-chart:

.. figure:: _image/chart.png
   :width: 600
   :align: center
   :figclass: align-center

   A chart of the nuclides with ENDF/B-VIII.0 elements. Isotopes in red are those included in the reduced version of the library, those in black are removed.

To model simple network with few isotopes and reactions, ONIX provides a set of convenient Python classes and methods to manually build custom networks. For more details on this capability see :ref:`utils`.


Isomeric branching
------------------

Isomeric branchings indicate the fraction of product isotopes that are in an excited states after a reaction. Most depletion software uses constant value for these branchings. However, isomeric branching values depend on the neutron spectrum and, therefore, change with the evolution of the neutronics in a nuclear system. In a coupled mode with OpenMC, ONIX update isomeric branching values each time a new neutron spectrum is calculated by OpenMC. This allows for a much more accurate simulation of production of isomeres (isotopes in an excited state). This is important for depletion calculations as isomeres have different nuclear properties than the ground states isotope.

So far, ONIX computes isomeric branchings for :math:`(n,\gamma)` reactions only as they account for the majority of isomeres production.

Coupling with OpenMC
--------------------

One of the main advantages of ONIX is that it can be readily coupled with OpenMC. OpenMC is first tasked with computing one-group reaction rates and neutron flux in different burnup regions chosen by the user. This data is then fed into ONIX which will deplete the system for a certain length of time (or burnup level) and produce new isotopic densities for each burnup region. These new isotopic densities are passed to OpenMC for a second neutron transport calculation and the process continues iteratively until the desired time (or burnup level) is reached. For more details about the coupling between ONIX and OpenMC, please see `Julien thesis <https://search.proquest.com/openview/7de190dd2bf7f8f7017fde115e462bfb/1?pq-origsite=gscholar&cbl=18750&diss=y>`_ and `ONIX paper <https://www.sciencedirect.com/science/article/pii/S0306454920306009>`_.

.. _fig-couple:

.. figure:: _image/couple.png
   :width: 600
   :align: center
   :figclass: align-center

   A diagram describing the coupling between ONIX and OpenMC.

Standalone mode
----------------

ONIX can also be used in a standalone mode, without coupling with OpenMC. In that mode, the user is responsible for providing one-group cross sections to ONIX (either through the form of a library or via dedicated Python classes). Currently, ONIX does not provide prepared one-group cross sections library. One-group cross sections libraries can be produced when running coupled simulations with OpenMC.

Output data processing
----------------------

ONIX's simulations produce many different type of output data that can be conveniently processed and visualized with dedicated Python functions. Output data produced by ONIX include:

    - the evolution of the multiplication factor of the whole system
    - the evolution of the isotopic densities for each burnup region
    - the evolution of the neutron flux spectrum for each burnup region
    - the evolution of the isomeric branchings for each burnup region
    - the evolution of the one-group cross sections for each burnup region
    - the evolution of the power for each burnup region
    - the evolution of the one-group neutron flux for each burnup region
    - snapshots of the burnup matrix for each burnup step
    - OpenMC input and output at each burnup step
    - (optional) production and destruction terms for all isotopes for each burnup region at each burnup step

ONIX Python API includes a set of functions to process and visualize these output results. For instance, the command :meth:`onix.utils.plot_bucell_nuclide_network` can plot a diagram of the network of production and destruction terms for one nuclide at a specific burnup step for a specific burnup region as seen on Figure 5.

.. fig-network:

.. figure:: _image/network.png
   :width: 600
   :align: center
   :figclass: align-center

   A diagram illustrating the production and destruction terms of plutonium-239 in a LWR reactor.

Scheduled changes in operation
------------------------------

During the operation of a nuclear reactor, many changes might occur in the composition, density and temperature of the materials that are not the product of depletion or decay. For instance, the operator of the reactor might change the boron density in the water during operation, the temperature of the cladding might rise which would lead to a change in its density etc. All these changes can be modelled in ONIX with the following command:

	- :meth:`onix.sequence.set_isotopic_change` for change in the isotopic composition
	- :meth:`onix.sequence.set_density_change` for change in the overall density of a material
	- :meth:`onix.sequence.set_temperature_change` for change in the temperature of a material

Being able to implement scheduled changes is important when comparing simulations to experimental benchmarks.

--------------------------------
Nuclear archaeology capabilities
--------------------------------

ONIX has been designed to support research in nuclear archaeology. The :strong:`NAX` module (:strong:`N`uclear :strong:`A`rchaeology (:strong:`X` ?) contains multiple Python modules and classes that can automate calculations and tasks useful to nuclear archaeology.

Nuclear archaeology
-------------------

.. fig-network:

.. figure:: _image/indiana.png
   :width: 600
   :align: center
   :figclass: align-center

   Nuclear archaeologists seek to discover the past mysteries of nuclear reactors

Nuclear archaeology (see `Fetter <https://www.tandfonline.com/doi/abs/10.1080/08929889308426386>`_) is a set of scientific methods and technics which goals is to verify the past production of weapon-usable materials in a country's nuclear program. Typically, nuclear archaeology is interested in estimating past production of plutonium and tritium in nuclear reactors or Highly Enriched Uranium (HEU) in enrichment plants. The way nuclear archaeology calculates past plutonium and tritium production is by estimating the fluence of a reactor (i.e., the integral of the neutron flux over time) since the production of both elements is proportional to fluence. In order to estimate fluence, nuclear archaeology relies on specific isotopic ratios measured directly from the reactor. These isotopic ratios, however, need to be good *fluence indicators*, i.e., the change of the ratio when exposed to neutrons must be big enough so that the uncertainties affecting the deduced fluence are acceptably small.

It must be noted that these isotopic ratios will more likely be measured from samples in structural materials and not in the irradiated fuel itself. This is why ONIX is capable of modelling the depletion of non-nuclear materials in the reactor.

The NAX module in ONIX
----------------------

The NAX module automates and simplifies the task of finding the best fluence indicators for a given reactor design with a given operation history. The tasks the module automates are the following:

	- [A] Ientifying isotope chains from the same element that have at least two stable or long-lived members and where at least one member has a non-negligible neutron cross section
	- [B] Depleting these chains according to an operation history defined by the user
	- [C] Selecting the isotope ratios from these chains that would allow to deduce fluence with the lowest uncertainties

Nuclear archaeology is often concerned with operation histories made up of multiple successive fuel batches spanning decades. For such a case, a single, continuous coupled simulation would be prohibitively expensive in terms of computing time. To circumvent this problem, ONIX only requires one coupled simulation per type of fuel batch that constitute the operation history. Using the one-group neutronics parameters obtained from each of these "model" simulations, the NAX module will then deplete structural materials with a fast analytical integrator according to the complete operation history (constituted of repetition of the different batch types).

.. fig-piece-wise:

.. figure:: _image/piece-wise.png
   :width: 600
   :align: center
   :figclass: align-center

   ONIX depletes strucutral materials for long operation history using prepared one-group parameters calculated for each model batch types


Example on a LWR reactor
------------------------

To briefly illustrate its functionality, the NAX module is used here for a fuel pin cell of a typical light-water reactor (LWR). The reactor considered operated with two different types of fuel batches: a fuel load of MOX fuel manufactured from reprocessed thermal reactor UO\ :sub:`2` spent fuel (type 1 batch) and a fuel load of MOX fuel manufactured from weapon plutonium (type 2 batch). The operation history for this example spans over 20 years. More details can be found on the design and operation of this reactor in `this report <https://www.oecd-nea.org/jcms/pl_17872>`_.

To be able to select the best fluence indicators, ONIX can plot graphs showing the evolution of the relative error on fluence estimation associated with different isotopic ratios. Figure 8 presents such a graph for the fuel pin cell of the LWR.

.. fig-fluence-rel-error:

.. figure:: _image/fluence-rel-error.png
   :width: 600
   :align: center
   :figclass: align-center

   Relative error on fluence that would be achieved by measuring various isotopic ratios. Type 2 batches are colored in grey.

------------
Applications
------------

ONIX can be used for many different applications, from nuclear reactor simulations, plutonium production verification to nuclear fallouts isotopics.

Nuclear reactor simulation
--------------------------

ONIX various functionalities and its coupling with the neutron transport code OpenMC provides a reliable and complete reactor physics package to undertake detailed reactor core simulations.

For instance, ONIX was used to simulate the assembly of a VVER type nuclear reactor which design can be found in `this report <https://www.oecd-nea.org/jcms/pl_17750>`_. In this assembly, several UO\ :sub:`2` fuel rods contain gadolinium. Gadolinium is a strong neutron absorbant and lead to important spatial self-shielding effects. It is crucial for a depletion code to accurately predict the spatial distribution of gadolinium in fuel rods in order to correctly account for spatial self-shielding. Figure 9 shows that ONIX is able to accurately model the spatial distribution of gadolinium in fuel rods containing the neutron absorber.

.. fig-gd-ring:

.. figure:: _image/gd-ring.png
   :width: 600
   :align: center
   :figclass: align-center

   Spatial distribution of gadolinium-155 (plain line) and gadolinium-157 (dashed line) in fuel rods. ONIX results are compared with other established reactor physics software.

Weapon-usable materials production
----------------------------------

ONIX is an ideal tool to model the production of fissile material like plutonium and other elements used in nuclear weapons. The ability to model :math:`(n,t)` reactions allows ONIX to model the production of tritium in a reactor. Tritium is an essential component for modern nuclear warheads. Figure 10 presents possible tritium production in the Dimona reactor where Israel likely produces tritium for its nuclear arsenal.

.. fig-tritium:

.. figure:: _image/tritium.png
   :width: 600
   :align: center
   :figclass: align-center

   Possible production of tritium in the Dimona reactor. Different enrichment in lithium-6 in rods' sleeves lead to different production of tritium.


Nuclear archaeology
-------------------

ONIX is the only software that provides modules to support research in nuclear archaeology. The NAX module can be used to identify the best fluence indicators that could, if measured, provide accurate estimation on production of plutonium or tritium. The NAX module of ONIX has been used for North Korea's 5-MWe reactor where the country likely produces plutonium and tritium for its nuclear arsenal. Figure 11 shows the relative error on fluence associated with different isotopic ratios along the operational history of the reactor. The best ratios for fluence estimations are those with the lowest relative error on fluence.

.. fig-dprk-pu:

.. figure:: _image/dprk-pu.png
   :width: 600
   :align: center
   :figclass: align-center

   Relative error on fluence that would be achieved by measuring various isotopic ratios for the 5-MWe reactor.

Nuclear fallouts
----------------

*Content coming soon*

----------
Benchmarks
----------

A more interactive presentation of ONIX benchmark validation will be uploaded soon.

Right now, you can find benchmark validation of ONIX on `Julien thesis <https://search.proquest.com/openview/7de190dd2bf7f8f7017fde115e462bfb/1?pq-origsite=gscholar&cbl=18750&diss=y>`_ and `ONIX paper <https://www.sciencedirect.com/science/article/pii/S0306454920306009>`_.





