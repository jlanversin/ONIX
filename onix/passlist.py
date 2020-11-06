"""Create list of passport, set mass, decay and xs"""
from .passport import Passport
from onix import data as d
import onix.utils as utils
import time

import math as m
import numpy as np
import os
import operator

class Passlist(object):

    """Passlist manages the list of nuclides' passports.

    It can for example order the list of nuclides in a certain order and is tasked with distributing nuclear data to individual nuclides.
    
    Parameters
    ----------
    nucl_list: list
        (Optiontal) List of nuclides for which Passlist should build the Passports list. A passlist object can be instantiated without providing a nuclides list

    """

    def __init__(self, *nucl_list):


        self._nucl_list = None
        self._passport_list = None

        if nucl_list:

            nucl_list = nucl_list[0]

            sample_elt = nucl_list[0]
            if utils.is_name(sample_elt):
                nucl_list = utils.name_list_to_zamid_list(nucl_list)

            passport_list = self.create_passport_list(nucl_list)

            self._nucl_list = nucl_list
            self._passport_list = passport_list
            self._set_mass(passport_list)
            self.zam_order_passport_list()
            self._set_zero_dens(passport_list)

    @property
    def nucl_list(self):
        """Returns the list of nuclides' names."""

        return self._nucl_list


    @property
    def passport_list(self):
        """Returns the list of Passport objects."""
        return self._passport_list

        
    def _get_zamid_passport_dict(self):
        """Convert the list of Passport objects into a dictionnary of Passport objects where keys are the zamid of the nuclides and entries are Passport objects."""

        passport_list = self.passport_list

        passport_dict = {}
        for i in passport_list:
            zamid = i.zamid
            passport_dict[zamid] = i

        return passport_dict

    def _get_name_passport_dict(self):
        """Convert the list of passport into a dictionnary of passports where entries are the zamid of the nuclides."""

        passport_list = self.passport_list

        passport_dict = {}
        for i in passport_list:
            name = i.name
            passport_dict[name] = i

        return passport_dict

    def azm_order_passport_list(self):
        """Orders the list of Passport objects after mass number first and then atomic number.
        """
        passport_list = self.passport_list
        print('AZM ORDER CALLED')
        changed = True
        while changed:
            changed = False
            for i in range(len(passport_list) - 1):
                a0 = passport_list[i].get_a
                a1 = passport_list[i+1].get_a
                z0 = passport_list[i].get_z
                z1 = passport_list[i+1].get_z
                state0 = passport_list[i].get_state()
                state1 = passport_list[i+1].get_state()
                azm0 = 1000*a0 + 10*z0 + state0
                azm1 = 1000*a1 + 10*z1 + state1
                if azm0 > azm1:
                    passport_list[i], passport_list[i+1] = passport_list[i+1], passport_list[i]
                    changed = True
        return None


    def zam_order_passport_list(self):
        """Orders the list of Passport objects after atomic number first and then mass number.
        """
        passport_list = self.passport_list
        changed = True
        while changed:
            changed = False
            for i in range(len(passport_list) - 1):
                passp0 = passport_list[i]
                passp1 = passport_list[i+1]
                zamid0 = int(passp0.zamid)
                zamid1 = int(passp1.zamid)

                if zamid0 > zamid1:
                    passport_list[i], passport_list[i+1] = passport_list[i+1], passport_list[i]
                    changed = True
        return None


    # def qsort(inlist):
    #     passport_list = self.passport_list
    #     if inlist == []: 
    #         return []
    #     else:
    #         pivot = inlist[0]
    #         lesser = qsort([x for x in inlist[1:] if x < pivot])
    #         greater = qsort([x for x in inlist[1:] if x >= pivot])
    #         return lesser + [pivot] + greater


    def zam_order_passport_list_2(self):
        """Orders the list of Passport objects after atomic number first and then mass number with Python list.sort() method.
        """
        passport_list = self.passport_list
        #passport_list = sorted(passport_list, key = lambda pp: int(pp.zamid))
        passport_list.sort(key = lambda pp: int(pp.zamid))

        #return ordered_passport_list


    def get_index_dict(self):
        """Returns a dictionnary where keys are z-a-m id of nuclides and entries are the index of their Passport in the Passport list.
        """
        passport_list = self.passport_list

        index_dict = {}

        index = 0
        for i in passport_list:
            index_dict[i.zamid] = index
            index = index + 1

        return index_dict

    def _add_nucl_list(self, nucl_list):

        sample_elt = nucl_list[0]
        if utils.is_name(sample_elt):
            nucl_list = utils.name_list_to_zamid_list(nucl_list)

        current_nucl_list = self.nucl_list
        new_extra_nuclides = list(set(nucl_list) - set(current_nucl_list))

        new_extra_passport_list = self.create_passport_list(new_extra_nuclides)
        self._set_mass(new_extra_passport_list)
        self._set_zero_dens(new_extra_passport_list)
        for nucl_pass in new_extra_passport_list:
            self.passport_list.append(nucl_pass)
            self.nucl_list.append(nucl_pass.zamid)

    def create_passport_list(self, nucl_list):
        """Create a Passport objects list after a provided list of nuclides.

        Parameters
        ----------
        nucl_list: list
            List of nuclides' names
        """
        passport_list = []
        for nucl in nucl_list:
            passport_list.append(Passport(nucl))

        return passport_list


    def _set_mass(self, passport_list):
        """Read and set the atomic mass for each nuclide in the passports list."""

        for i in range(len(passport_list)):

            nuc_pass = passport_list[i]
            zaid = nuc_pass.zamid[:-1]

            # Some (rare) nuclides (for example Ag133 from ENDFVIII) are not in the mass library
            # For these nuclides, the mass number is used for the mass
            if zaid not in d.default_atm_mass_lib:
                mass = int(zaid[-3:])

            else:
                mass = d.default_atm_mass_lib[zaid]
            
            nuc_pass.mass = mass


    def _set_decay(self, decay_lib_b, decay_lib_a):
        """Read and set the decay constants for each nuclide in the passports list."""

        passport_list = self.passport_list

        decay_b_dict = decay_lib_b
        decay_a_dict = decay_lib_a

        for i in range(len(passport_list)):

            nuc_pass = passport_list[i]
            zamid = nuc_pass.zamid
            if zamid in decay_b_dict:
                decay_b = decay_b_dict[zamid]
                decay_a = decay_a_dict[zamid]
                nuc_pass.set_decay(decay_a, decay_b)

    def _set_xs(self, xs_dict):
        """Read and set the cross sections for each nuclide in the passports list"""

        passport_list = self.passport_list

        for i in range(len(passport_list)):

            nuc_pass = passport_list[i]
            zamid = nuc_pass.zamid
            if zamid in xs_dict:
                xs = xs_dict[zamid]
                nuc_pass._set_xs(xs)

    def _overwrite_xs(self, xs_dict):
        """Read and set the cross sections for each nuclide in the passports list"""

        passport_list = self.passport_list

        for i in range(len(passport_list)):

            nuc_pass = passport_list[i]
            zamid = nuc_pass.zamid
            if zamid in xs_dict:
                xs = xs_dict[zamid]
                # If there is already a XS set, overwrite with None
                # _set_xs(xs) will automatically reinitialize xs_seq
                if nuc_pass.current_xs != None:
                    nuc_pass.current_xs = None
                nuc_pass._set_xs(xs)


    def _set_fy(self, fy_dict):
        """Read and set the fission yields for fission products in the passports list"""

        passport_list = self.passport_list

        for i in range(len(passport_list)):

            nuc_pass = passport_list[i]
            zamid = nuc_pass.zamid
            if zamid in fy_dict:
                fy = fy_dict[zamid]
                nuc_pass.fy = fy

    def _set_initial_dens(self, dens_dict):

        sample_key = list(dens_dict.keys())[0]
        passport_list = self.passport_list

        for i in range(len(passport_list)):
            nuc_pass = passport_list[i]

            if utils.is_zamid(sample_key):
                zamid = nuc_pass.zamid
                if zamid in dens_dict:
                    dens = dens_dict[zamid]
                    nuc_pass._set_initial_dens(dens)
                else:
                    nuc_pass._set_initial_dens(0.0)

            elif utils.is_name(sample_key):
                name = nuc_pass.name
                if name in dens_dict:
                    dens = dens_dict[name]
                    nuc_pass._set_initial_dens(dens)
                else:
                    nuc_pass._set_initial_dens(0.0)


    def _set_zero_dens(self, passport_list):

        for nucl in passport_list:
            nucl._set_initial_dens(0.0)
        


    # if the fission yield are updated, this should be regenerated
    # This function only add a child in the fission child list of parents if the fy is non-zero
    # Another way of doing would be to just add all fission_child from fy_lib to the child list of parents
    def _set_fission_child(self, fy_nucl, fy_parent):

        passport_dict = self._get_zamid_passport_dict()

        fission_child = fy_nucl
        fission_parent = fy_parent
        fission_child_per_parent = {}

        for parent in fission_parent:
            fission_child_per_parent[parent] = []

        for zamid in fission_child:
            nuc_pass = passport_dict[zamid]
            nuc_fy = nuc_pass.fy

            for parent in nuc_fy:
                if nuc_fy[parent][0] != 0:
                    fission_child_per_parent[parent].append(zamid)

        # Set the fission child list for each parent passport
        for parent in fission_child_per_parent:
            nuc_pass = passport_dict[parent]
            child_list = fission_child_per_parent[parent]
            nuc_pass.fission_child = child_list


    def neg_reac_warning(passport_list):

        print('NEGATIVE REACTION WARNING CALLED')

        for i in passport_list:
            nuc_name = i.name
            decay_a = i.decay_a
            xs_dic = i.xs

            if decay_a != None:

                for decay in decay_a:
                    decay_val = decay_a[decay]
                    if decay_val < 0:
                        raise Neg_decay('{} has its {} reaction negative with value {}'.format(nuc_name, decay, decay_val))

            if xs_dic != None:
        
                for xs in xs_dic:
                    xs_val = xs_dic[xs]
                    if xs_val < 0:
                        raise Neg_xs('{} has its {} reaction negative with value {}'.format(nuc_name, xs, xs_val))


    def order_name(self):

        passport_list = self.passport_list

        azm_sort(passport_list)

        N = len(passport_list)
        L = [None]*N

        for i in range(N):
            nuc_pass = passport_list[i]
            L[i] = nuc_pass.name

        return L















class Nuc_xs_not_found(Exception):
    """Raise when the user requests a cross-sections of a nuclide that is not in the nuclide set """
    pass

class Neg_decay(Exception):
    """Raise when a negative decay constant is found"""
    pass

class Neg_xs(Exception):
    """Raise when a negative cross-section is found"""
    pass
