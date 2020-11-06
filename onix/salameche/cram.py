"""Compute the solution of the matricial depletion equation using the CRAM method"""
import numpy as np
import warnings
import time

def CRAM16(At,N_0):
    """CRAM uses a Chebishev Rational Approximation Method of order 16 to compute the solution of the matricial depletion equation.

    Parameters
    ----------
    At: numpy.array
        Depletion matrix multiplied by the time interval over which nuclides are depleted
    N_0: numpy.array
        Initial nuclides' densities vector
    """

    print ('CRAM CALLED')
    t0 = time.time()

    lN = len(N_0)

    theta = np.array([
    -1.0843917078696988026e1 +1.9277446167181652284e1j,
    -5.2649713434426468895 +1.6220221473167927305e1j,
    +5.9481522689511774808 +3.5874573620183222829j,
    +3.5091036084149180974 +8.4361989858843750826j,
    +6.4161776990994341923 +1.1941223933701386874j,
    +1.4193758971856659786 +1.0925363484496722585e1j,
    +4.9931747377179963991 +5.9968817136039422260j,
    -1.4139284624888862114 +1.3497725698892745389e1j], dtype = np.complex256)

    alpha_0 = np.complex256(2.1248537104952237488e-16 + 0.0j) 

    alpha = np.array([
    -5.0901521865224915650e-7 -2.4220017652852287970e-5j,
    +2.1151742182466030907e-4 +4.3892969647380673918e-3j,
    +1.1339775178483930527e2 +1.0194721704215856450e2j,
    +1.5059585270023467528e1 -5.7514052776421819979j,
    -6.4500878025539646595e1 -2.2459440762652096056e2j,
    -1.4793007113557999718 +1.7686588323782937906j,
    -6.2518392463207918892e1 -1.1190391094283228480e1j,
    +4.1023136835410021273e-2 -1.5743466173455468191e-1j], dtype = np.complex256)

    l = len(theta)
    N = N_0*0
    _N = np.zeros((lN),dtype=np.complex128)

    for i in range(l):
        term1 = At - theta[i]*np.identity(np.shape(At)[0])
        term2 = alpha[i]*N_0
        _N += np.linalg.solve(term1,term2)
        
    N = 2*_N.real
    N = N + alpha_0*N_0
    # For some reason here N is still complex and not only real

    print('CRAM took:{} s'.format(time.time() - t0))

    return N.real

# CRAM is yielding non zero values for nuclides that should be at zero because no one is producing them
# This algorithm check which nuclide are in this situation and set their density to zero
def CRAM_reality_check(bucell, index_dic, N):
    """This functions checks against negative and extremely small densities (same as onix.salameche.CRAM_density_check). In addition, it compares the calculated new densities with the BUCell.leave attribute (this attribute enables to know which isotopes should be produced or not during depletion). The comparison enables the function to detect nuclides that have a non-zero density but should have a zero density. Likewise, it can detect nuclides that should be produced but have zero density.

    Parameters
    ----------
    bucell: onix.Cell
        BUCell being depleted
    index_dict: dict
        Dictionnary where keys are nuclides z-a-m id and entries are their indexes in the density vector
    N: numpy.array
        New density vector solution to the depletion equation
    """

    print('reality check called')
    passlist = bucell.passlist
    leaves = bucell.leaves
    fission_leaves = bucell.fission_leaves
    total_leaves = leaves + fission_leaves

    negative_count = 0
    small_count = 0
    intruder_count = 0
    missing_count = 0

    for nuc_pass in passlist:
        nuc_zamid = nuc_pass.zamid
        nuc_name = nuc_pass.name
        nuc_dens = nuc_pass.dens
        index = index_dic[nuc_pass.zamid]
        N_val = N[index]

        # if nuc_name == 'Au-200':
        #     print nuc_name, nuc_zamid
        #     print nuc_dens

        if N_val < 0:
          #  warnings.warn('NEGATIVE: Nuclide {}/{} has a negative density of {}'.format(nuc_name, nuc_zamid, N_val))
            N[index] = 0.0
            negative_count += 1

        N_val = N[index]

        if N_val < 1e-24:
          #  warnings.warn('TOO SMALL: Nuclide {}/{} has a density of {} below 1e-24'.format(nuc_name, nuc_zamid, N_val))
            N[index] = 0.0
            small_count += 1

        N_val = N[index]

        if nuc_zamid in total_leaves and N_val == 0:

          #  warnings.warn('MISSING: Nuclide {} has a density of 0 while it belongs to the creation tree'.format(nuc_name))
            missing_count += 1

        elif nuc_zamid not in total_leaves and N_val != 0:

          #  warnings.warn('INTRUDER: Nuclide {} has a density of {} while it is not in the creation tree'.format(nuc_name, N_val ))
            intruder_count += 1

    print(('There are {} negative'.format(negative_count)))
    print(('There are {} too small'.format(small_count)))
    print(('There are {} intruders'.format(intruder_count)))
    print(('There are {} missings'.format(missing_count)))

def CRAM_density_check(bucell, N):
    """This function checks for extremely low densities and negative densities. Densities below one atom per cubic centimeter are set to zero. Negative densities are produced by mathematical approximations inherent to the CRAM method and therefore do not bear any physical meaning. They are also set to zero.

    Parameters
    ----------
    bucell: onix.Cell
        BUCell being depleted
    N: numpy.array
        New density vector solution to the depletion equation
    """
    passlist = bucell.passlist
    index_dict = passlist.get_index_dict()
    passport_list = passlist.passport_list

    negative_count = 0
    small_count = 0

    for nuc_pass in passport_list:
        nuc_zamid = nuc_pass.zamid
        nuc_name = nuc_pass.name
        index = index_dict[nuc_pass.zamid]
        N_val = N[index]

        if N_val < 0:
            #warnings.warn('NEGATIVE: Nuclide {}/{} has a negative density of {}'.format(nuc_name, nuc_zamid, N_val))
            N[index] = 0.0
            negative_count += 1

        N_val = N[index]

        if N_val < 1e-24:
            #warnings.warn('TOO SMALL: Nuclide {}/{} has a density of {} below 1e-24'.format(nuc_name, nuc_zamid, N_val))
            N[index] = 0.0
            small_count += 1

    print(('There are {} negative'.format(negative_count)))
    print(('There are {} too small'.format(small_count)))





