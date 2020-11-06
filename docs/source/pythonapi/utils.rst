--------------------------------------
:mod:`onix` -- Utilities and functions
--------------------------------------

Reaction classes
-----------------

These classes enable users to build their own custom nuclear libraries. Instantiations of these classes can be used as nuclear libraries for simulations.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myclass.rst

   onix.utils.decay_lib
   onix.utils.xs_lib
   onix.utils.fy_lib

Functions
---------

Convenient functions to be used directly by users and used in the source code.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   onix.utils.decay_to_halflife
   onix.utils.halflife_to_decay
   onix.utils.halflife_to_second
   onix.utils.is_int
   onix.utils.is_zamid
   onix.utils.is_name
   onix.utils.is_list_redundant
   onix.utils.get_list_redundant_elt
   onix.utils.zamid_list_to_name_list
   onix.utils.name_list_to_zamid_list
   onix.utils.get_zamid_z
   onix.utils.get_name_z
   onix.utils.get_zamid_a
   onix.utils.get_zamid_n
   onix.utils.get_zamid_s
   onix.utils.zamid_to_name
   onix.utils.name_to_zamid
   onix.utils.get_hm
   onix.utils.get_nucl_atomic_mass
   onix.utils.convert_mass_to_atom
   onix.utils.convert_atom_to_mass
   onix.utils.get_bu_sec_conv_factor
   onix.utils.get_decay_nucl
   onix.utils.get_xs_nucl
   onix.utils.get_fy_nucl
   onix.utils.is_lista_in_listb
   onix.utils.get_fy_parent_nucl
   onix.utils.is_number
   onix.utils.openmc_name_to_onix_name
   onix.utils.onix_name_to_openmc_name
   onix.utils.bu_namelist_to_mc_namelist
   onix.utils.mc_namelist_to_bu_namelist
   onix.utils.order_nuclide_per_z
   onix.utils.order_nuclide_name_per_z
   onix.utils.order_nuclide_per_a
   onix.utils.get_openmc_xs_nucl_list
   onix.utils.convert_spectrum_to_janis_weighting_format
   onix.utils.get_zamid_natural_abundance
   onix.utils.get_name_natural_abundance
   onix.utils.find_zamid_precursor
   onix.utils.interpolation_between_two_points



Data processors
---------------

Various functions to extract, read and visualize output results.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: myfunction.rst

   onix.utils.plot_bucell_nuclide_network
   onix.utils.plot_nuclide_dens_from_path
   onix.utils.plot_nuclide_group_dens_from_path
   onix.utils.plot_xs_time_evolution_from_path
   onix.utils.plot_xs_bu_evolution_from_path
   onix.utils.compare_xs_bu_evolution_from_path
   onix.utils.plot_kinf_from_path
   onix.utils.plot_flux_from_path
   onix.utils.plot_flux_spectrum_bu_evolution_from_path
   onix.utils.plot_lethargy_spectrum_bu_evolution_from_path
   onix.utils.plot_xs_dens_flux
   onix.utils.read_time_seq
   onix.utils.get_step_time_length_seq
   onix.utils.read_bu_seq
   onix.utils.read_kinf_seq
   onix.utils.read_flux
   onix.utils.read_flux_subseq
   onix.utils.get_fluence_seq
   onix.utils.get_fluence_subseq
   onix.utils.read_flux_spectrum
   onix.utils.read_energy_mid_points
   onix.utils.read_energy_bin_length
   onix.utils.read_dens
   onix.utils.get_total_density
   onix.utils.get_total_mass_density
   onix.utils.read_xs_seq
   onix.utils.get_time_averaged_xs
   onix.utils.get_time_averaged_flux
   onix.utils.get_tot_xs
   onix.utils.read_xs_nucl
   onix.utils.read_dens_nucl
   onix.utils.rank_nuclide_per_dens
   onix.utils.plot_matrix_bysign_from_compressed_matrix
   onix.utils.plot_nuclide_chart_color_per_nuclear_data
   onix.utils.plot_compare_two_nuclear_data_on_nuclide_chart
