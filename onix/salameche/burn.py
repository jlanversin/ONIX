
import numpy as np

from onix import data
from onix import passlist as pl
from . import mat_builder as mb
from . import cram
from . import py_pade   

# This function is in fact not used anymore.
# Might need to be removed
def burn(system):
    """Depletes a System Object according to its associated Sequence object.
    """
    sequence = system.sequence
    steps_number = sequence.steps_number
    # The norm in onix is that Step 0 is reserved for initial status
    # Step 1 is the first true step

    # This normalization resembles the step_normalization operated in
    # coupled mode. However, here, the flux density is the same in each 
    # cell (flux density set for system as a whole)
    # This normalization will just set the step flux and step pow dens
    # for each cell
    system.step_normalization()

    for s in range(1, steps_number+1):

        print ('\n\n\n\n STEP {}\n\n\n\n'.format(s))

        sequence._gen_step_folder(s)
        burn_step(system, s, mode='stand alone')

    system._gen_output_summary_folder()
    system._print_summary_allreacs_rank()
    system._print_summary_subdens()
    system._print_summary_dens()

def burn_step(system, s, mode):

    """Depletes the system for macrostep s.

    Parameters
    ----------
    system: onix.System
        System to be depleted
    s: int
        Macrostep number
    mode: str
        'stand alone' or 'couple'
    """

    bucell_list = system.get_bucell_list()
    reac_rank = system.reac_rank
    for bucell in bucell_list:
        print ('\n\n\n\n CELL {}\n\n\n\n'.format(bucell.name))
        # Create and set the folder corresponding to that cell
        bucell._set_folder()
        bucell._print_xs_lib()
        burn_cell(bucell, s, mode, reac_rank)
        bucell._change_isotope_density(s)
        bucell._change_total_density(s)
        bucell._print_substep_dens(s)

    if reac_rank == 'on':
        system._print_current_allreacs_rank()
    system._copy_cell_folders_to_step_folder(s)

def burn_cell(bucell, s, mode, reac_rank):

    """Depletes a BUCell for macrostep s.

    Parameters
    ----------
    bucell: onix.Cell
        BUCell to be depleted
    s: int
        Macrostep number
    mode: str
        'stand alone' or 'couple'
    reac_rank: str
        'on' or 'off'. If set to 'on', each BUCell will compute production and destruction terms ranking for each nuclide at every macrostep.
    """
    # Check if different nuclide lists are coherent with each other
    # This should not be called in burn_cell because bun_cell is in the sequence loop
    # This function should be only called once before burn_cell
    #bucell._check_nucl_list_consistency()

    # Initial print with initial density and initial time and bu
   # cell._print_initial_dens()

    passlist = bucell.passlist
    sequence = bucell.sequence
    flux_approximation = sequence.flux_approximation
    # microsteps_number is of length s-1
    microsteps_number = sequence.microsteps_number(s-1)
    
    B = mb.get_xs_mat(passlist)
    C = mb.get_decay_mat(passlist)
    N = mb.get_initial_vect(passlist)

    # Store B and C on compressed txt
    # there needs to be a new matrix print for every step (even substep)
    mb._print_all_mat_to_text(B, C, bucell, s)

    if flux_approximation == 'iv':
        for i in range(microsteps_number):
            N = burn_microstep(bucell, B, C, N, s, i, microsteps_number, mode, reac_rank)
    elif flux_approximation == 'pc':
        for i in range(microsteps_number):
            burn_substep_pc(bucell, B, C, N, s, i, microsteps_number, mode)
    elif flux_approximation == 'me':
        for i in range(microsteps_number):
            burn_substep_pcME4(bucell, B, C, N, s, i, microsteps_number, mode)

    bucell._set_step_dens()
    sequence._set_macrostep_bucell_bu()

    # At the end of this burn sequence, the flux and power

def burn_microstep(bucell, B, C, N, s, ss, ssn, mode, reac_rank):

    """Depletes a BUCell for microstep ss within macrostep s.

    Parameters
    ----------
    bucell: onix.Cell
        BUCell to be depleted
    B: numpy.array
        Neutron-induced reaction transmutation matrix
    C: numpy.array
        Decay matrix
    s: int
        Macrostep number
    ss: int
        Microstep number
    ssn: int
        Total number of microstep within macrostep s
    mode: str
        'stand alone' or 'couple'
    reac_rank: str
        'on' or 'off'. If set to 'on', each BUCell will compute production and destruction terms ranking for each nuclide at every macrostep.
    """

    print ('\n\n++++ Microstep {} ++++\n\n'.format(ss))

    bucell_id = bucell.id
    sequence = bucell.sequence
    norma = sequence.norma_unit
    #norma_value = sequence.norma_vector[s]

    # I decided that since pow dens sequence will only be updated dynamically
    # and not created at the beginning, then the conversion bu-time/time-bu
    # will also be done dynamically
    # In case of norma = bu, time need to be calculated now, before CRAM
    sequence._bucell_time_bu_substep_conversion(bucell, s, ss)

    # print ('s', s, 'ss', ss)
    # print ('time subseq mat', sequence._time_subseq_mat)
    # print ('bu subseq mat', sequence._bucell_bu_subseq_mat)
    time_point = sequence.time_subpoint(s, ss)
    bucell_bu_point = sequence.bucell_bu_subpoint(s, ss)
    time_substep = sequence.get_time_subintvl(s, ss)

    pow_dens = sequence.current_pow_dens
    flux = sequence.current_flux

    # If this is the first substep, no need to update the flux or pow dens
    # both have been already calculated 
    if ss != 0:
    # Check whether actinides are present in the cell
    # If not, then no flux/power update should be done
        act = bucell.check_act_presence()
        if act == 'yes':
        # Now that the density of nuclides is updated, calculate the new substep flux or power density
            if norma == 'power':
                flux = bucell._update_flux(pow_dens)

            elif norma == 'flux':
                pow_dens = bucell._update_pow_dens(flux)

    sequence._set_substep_flux(flux, s, ss)
    sequence._set_substep_pow_dens(pow_dens, s, ss)

    # print(('current_time_point', time_point))
    # print(('current bucell_bu', bucell_bu_point))
    A = (B*1e-24*flux + C)
    #A = B*1e-24*flux
    # print(('current flux', flux))
    # print(('current time_subintvl', time_substep))

    At = A*time_substep

    N = cram.CRAM16(At, N)
    #N = py_pade.pade(At, N)

    cram.CRAM_density_check(bucell, N)

    bucell._update_dens(N, ss, ssn)

    # Generate the allreacsdic for each nuclide
    if reac_rank == 'on':
        # print ('pika')
        # quit()
        bucell._set_allreacs_dic(s, ss, ssn)

    return N



























def burn_substep_pc(cell, B, C, N, s, i): # with predictor corrector

    cell_id = cell.id
    passlist = cell.passlist
    index_dic = cell.index_dic
    sequence = cell.sequence
    norma = sequence.norma_unit
    norma_value = sequence.norma_vector[s]
    power = norma_value

    time_point = sequence.time_subpoint(s,i)
    bu_point = sequence.bu_subpoint(s,i)
    time_substep = sequence.time_substep(s,i)

    Ni = N.copy()

    flux_pc = []

    for pc_i in range(2):

        flux_pc.append(sequence._update_flux(passlist, power))

        if pc_i ==1: # Corrector step

            flux_pc = (flux_pc[0] + flux_pc[1])/2

            # Pass the flux to sequence. If power norma, seq and subseq are updated. If flux norma, only subseq needs to be updated
            sequence._set_flux(flux_pc, s, i)

            # Generate the first set of reacs
            if s == 0 and i == 0:
                pl._set_allreacs_dic(passlist, sequence, s, i)


        A = (B*1e-24*flux_pc + C)
        #A = B*1e-24*flux
        print(('flux_pc', flux_pc))
        print(('time_substep', time_substep))

        At = A*time_substep

        N = cram.CRAM16(At, Ni)

        cram.CRAM_reality_check(cell, index_dic, N)

        pl._update_dens(N, passlist)

        if pc_i == 1:

          #  pl.print_dens_2(cell, s, i)

            # Store B and C on compressed txt
            # there needs to be a new matrix print for every step (even substep)
            mb._print_all_mat_to_text(B, C, flux_pc, cell, s, i)

            # Generate the allreacsdic for each nuclide
            pl._set_allreacs_dic(passlist, sequence, s, i)


def burn_substep_pcME4(cell, B, C, N, s, i): # with predictor corrector

    cell_id = cell.id
    passlist = cell.passlist
    index_dic = cell.index_dic
    sequence = cell.sequence
    norma = sequence.norma_unit
    norma_value = sequence.norma_vector[s]
    power = norma_value

    time_point = sequence.time_subpoint(s,i)
    bu_point = sequence.bu_subpoint(s,i)
    time_substep = sequence.time_substep(s,i)

    Ni = N.copy()

    flux_pc = []

    # Predictor
    print('Predictor')

    flux_pc.append(sequence._update_flux(passlist, power))

    A = (B*1e-24*flux_pc + C)
    #A = B*1e-24*flux
    print(('flux_pc', flux_pc))
    print(('time_substep', time_substep))

    At = A*time_substep

    N = cram.CRAM16(At, Ni)

    cram.CRAM_reality_check(cell, index_dic, N)

    pl._update_dens(N, passlist)



    # Corrector
    print('Corrector')

    flux_pc.append(sequence._update_flux(passlist, power))
    flux_av = (flux_pc[0] + flux_pc[1])/2

    # Pass the flux to sequence. If power norma, seq and subseq are updated. If flux norma, only subseq needs to be updated
    sequence._set_flux(flux_av, s, i)

    # Generate the first set of reacs
    if s == 0 and i == 0:
        pl._set_allreacs_dic(passlist, sequence, s, i)


    # Apct
    Apc = (B*1e-24*flux_av + C)
    #A = B*1e-24*flux
    print(('flux_pc', flux_pc))
    print(('time_substep', time_substep))
    Apct = A*time_substep

    # Correction from ME4
    BC_CB = np.matmul(B,C) - np.matmul(C,B)
    correction = BC_CB*((flux_pc[1] - flux_pc[0])*time_substep**2)/12

    Apct_corrected = Apct + correction

    np.set_printoptions(threshold='nan')
    print(("Apct", Apct))
    print(("determinant Apct", np.linalg.det(Apct)))
    print(("norm At", np.linalg.norm(Apct, 2)))

    print(("BC_CB", BC_CB))
    print(("determinant BC_CB", np.linalg.det(BC_CB)))
    print(("norm BC_CB", np.linalg.norm(BC_CB, 2)))




    print(("Correction", correction))
    print(("Apct_corrected", Apct_corrected))
    print(('determinant corrected',   np.linalg.det(Apct)))
    norm_corrected = np.linalg.norm(Apct_corrected, 2)
    print(('NORM corrected', norm_corrected))

    N = cram.CRAM16(Apct_corrected, Ni)

    cram.CRAM_reality_check(cell, index_dic, N)

    pl._update_dens(N, passlist)

  #  pl.print_dens_2(cell, s, i)

    # Store B and C on compressed txt
    # there needs to be a new matrix print for every step (even substep)
    mb._print_all_mat_to_text(B, C, flux_pc, cell, s, i)

    # Generate the allreacsdic for each nuclide
    pl._set_allreacs_dic(passlist, sequence, s, i)
























