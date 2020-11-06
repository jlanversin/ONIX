-------------------------------------------
:mod:`onix` -- Salamèche: the burnup solver
-------------------------------------------

Burn
----

These functions operate the depletion of the System and the BUCells for each microsteps and macrosteps.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   onix.salameche.burn_step
   onix.salameche.burn_cell
   onix.salameche.burn_microstep


Matrix builder
--------------

These functions build the depletion matrix.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   onix.salameche.get_xs_mat
   onix.salameche.get_decay_mat
   onix.salameche.get_initial_vect

CRAM solver
-----------

These functions implement the CRAM solver and verify the consistency of calculated density values.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   onix.salameche.CRAM16
   onix.salameche.CRAM_reality_check
   onix.salameche.CRAM_density_check

