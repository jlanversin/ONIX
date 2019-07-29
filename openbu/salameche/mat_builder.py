""" Uses the passport list to build the transmutation matrix"""
import numpy as np
import os

def _get_xs_mat(passlist):
    """ Build the cross section matrix"""

    # order the list by mass number
    #azm_sort(passport_list)

    # create a dic that stores the index of each nuclide in the list
    passport_list = passlist.passport_list
    index_dict = passlist.get_index_dict()

    N = len(passport_list)
    xs_mat = np.zeros((N,N))

    for row in range(N):
        nucli = passport_list[row]
        zaidi = nucli.zamid
        FAM = nucli.get_FAM()

        # Diagonal term
        nucli_xs = nucli.current_xs
        if nucli_xs != None:
            xs_mat[row][row] = -nucli_xs['removal'][0]
        # Non diagonal terms (non fission)
        xs_parent = nucli.xs_parent
        for j in xs_parent:
            # If the parent nuclide is not in passport_list, skip
            if xs_parent[j] not in index_dict:
                continue
            index = index_dict[xs_parent[j]]
            parent_pass = passport_list[index]
            parent_xs = parent_pass.current_xs
            # If there are no xs for the parent nuclide, skip
            if parent_xs == None:
                continue

            # In the xs_parent dic, reaction that start from an excited state have an 'X' at the beginning of the reaction name. This 
            # is not the case in the xs dic. Therefore, we need to remove the starting X to match the reaction name in parent_xs
            if parent_pass.state == 1:
                j = j[1:]

            # If the specific xs that produce nucli does not exist for the parent, skip
            if j not in parent_xs:
                continue

            xs_val = parent_xs[j][0]
            # If the xs that produces nucli is zero, skip
            if xs_val == 0.0:
                continue
            xs_mat[row][index] = xs_val



        ####################### WARNING, TO BE REMOVED ##########################
        ######FOR ARGONNE BENCHMARK ONLY

        # if zaidi =='611490':
        #     indexpm148X = index_dict['611481']
        #     pm148X_passport = passport_list[indexpm148X]
        #     indexpm148X_xs_val = pm148X_passport.current_xs['ngamma'][0]
        #     xs_mat[row][indexpm148X] = indexpm148X_xs_val

        # if zaidi =='952430':
        #     indexam242X = index_dict['952421']
        #     am242X_passport = passport_list[indexam242X]
        #     indexam242X_xs_val = am242X_passport.current_xs['ngamma'][0]
        #     xs_mat[row][indexam242X] = indexam242X_xs_val
        ####################### WARNING, TO BE REMOVED ##########################


        # fission products yields  
        if FAM == 'FP':
            fy = nucli.fy
            for j in fy:
                # If the parent nuclide not in passport_list, skip
                if j not in index_dict:
                    continue
                # If the fission yield is zero, skip
                if fy[j][0] == 0.0:
                    continue
                index = index_dict[j]
                parent_pass = passport_list[index]
                parent_xs = parent_pass.current_xs
                # If there are no xs for the parent nuclide, skip
                if parent_xs == None:
                    continue
                # if fission is not in parent cross section
                if 'fission' not in parent_xs:
                    continue
                fission_val = parent_xs['fission'][0]
                # If the fission that produces nucli is zero, skip
                if fission_val == 0.0:
                    continue
                val = fission_val*fy[j][0]*1e-2

                # if j == '922350' and zaidi == '571490':

                #     print ('PIIIIIIIIIKKKKKAAAAAAAAAAAAAAAAA')
                #     print (fission_val)
                #     print (fy[j][0])
                #     print (val)

                xs_mat[row][index] = val

    return xs_mat






def _get_decay_mat(passlist):
    """ Build the cross section matrix"""

    # order the list by mass number
    #azm_sort(passport_list)

    # create a dic that stores the index of each nuclide in the list
    passport_list = passlist.passport_list
    index_dict = passlist.get_index_dict()

    N = len(passport_list)
    
    decay_mat = np.zeros((N,N))

    for row in range(N):
        nucli = passport_list[row]
        zaidi = nucli.zamid


        # Diagonal term
        nucli_decay = nucli.decay_a
        if nucli_decay != None and nucli_decay != 'stable':
            decay_mat[row][row] = -nucli_decay['total decay']

        # Non diagonal terms
        decay_parent = nucli.decay_parent
        for j in decay_parent:
            # If the parent nuclide is not in passport_list, skip
            if decay_parent[j] not in index_dict:
                continue

            index = index_dict[decay_parent[j]]
            parent_pass = passport_list[index]
            parent_decay = parent_pass.decay_a
            # If there are no decay for the parent nuclide, skip
            if parent_decay == None or parent_decay == 'stable':
                continue

            # In the decay_parent dic, reaction that start from an excited state have an 'X' at the beginning of the reaction name. This 
            # is not the case in the parent_decay dic. Therefore, we need to remove the starting X to match the reaction name in parent_decay
            if parent_pass.state == 1:
                j = j[1:]

            # If i is not in parent_decay
            if j not in parent_decay:
                continue

            decay_val = parent_decay[j]
            # If the decay that produces nucli is zero, skip
            if decay_val == 0.0:
                continue
            decay_mat[row][index] = decay_val

    return decay_mat


def _print_all_mat_to_text(xs_mat, decay_mat, cell, s):

    mat_folder_path = _get_mat_folder_path(cell)

    _gen_mat_folder(mat_folder_path)

    print ('get_xs_mat_text_1')
    xs_mat_txt = _get_xs_mat_text_1(xs_mat, cell)
   # xsphi_mat_txt = _get_xs_mat_text_2(xs_mat, cell, flux)
    print ('get_decay_mat_text')
    decay_mat_txt = _get_decay_mat_text(decay_mat, cell)

    xs_mat_name = mat_folder_path + '/xs_mat'
 #   xsphi_mat_name = mat_folder_path + '/{}_xsphi_mat'.format(step_point)
    decay_mat_name = mat_folder_path + '/decay_mat'

    xs_mat = open(xs_mat_name, 'w')
  #  xsphi_mat = open(xsphi_mat_name, 'w')
    decay_mat = open(decay_mat_name, 'w')

    xs_mat.write(xs_mat_txt)
  #  xsphi_mat.write(xsphi_mat_txt)
    decay_mat.write(decay_mat_txt)

    xs_mat.close()
  #  xsphi_mat.close()
    decay_mat.close()


def _gen_mat_folder(path):

    mat_folder_path = path

    if os.path.exists(mat_folder_path):
        shutil.rmtree(mat_folder_path)

    os.makedirs(mat_folder_path)


def _get_mat_folder_path(cell):

    cell_folder_path = cell.folder_path
    mat_folder = 'matrix'

    mat_folder_path = cell_folder_path + '/' + mat_folder

    return mat_folder_path


def _get_xs_mat_text_1(xs_mat, cell):

    passlist = cell.passlist
    passport_list = passlist.passport_list
    N = len(passport_list)
    Btxt = ''
    for row in range(N):
        nuc_pass = passport_list[row]
        nuc_zamid = nuc_pass.zamid
        Btxt += '{}|{}:'.format(nuc_zamid, row)

        # Diag terms
        diag_val = xs_mat[row][row]
        Btxt += ' {} {},'.format(row, diag_val)

        # Non diag terms
        for col in range(N):
            if col == row:
                continue
            if xs_mat[row][col] != 0.0:
                Btxt += ' {} {},'.format(col, xs_mat[row][col])

        Btxt += '\n'

    return Btxt

    # B = open('xs_mat', 'w')
    # B.write(Btxt)

def _get_xs_mat_text_2(xs_mat, cell, flux):

    passlist = cell.passlist
    passport_list = passlist.passport_list
    N = len(passport_list)
    Btxt = ''
    for row in range(N):
        nuc_pass = passport_list[row]
        nuc_zamid = nuc_pass.zamid
        Btxt += '{}|{}:'.format(nuc_zamid, row)

        # Diag terms
        diag_val = xs_mat[row][row]*flux*1e-24
        Btxt += ' {} {},'.format(row, diag_val)

        # Non diag terms
        for col in range(N):
            if col == row:
                continue
            if xs_mat[row][col] != 0.0:
                Btxt += ' {} {},'.format(col, xs_mat[row][col]*flux*1e-24)

        Btxt += '\n'

    return Btxt
    # B = open('xsphi_mat', 'w')
    # B.write(Btxt)


def _get_decay_mat_text(decay_mat, cell):


    passlist = cell.passlist
    passport_list = passlist.passport_list
    N = len(passport_list)
    Ctxt = ''
    for row in range(N):
        nuc_pass = passport_list[row]
        nuc_zamid = nuc_pass.zamid
        Ctxt += '{}|{}:'.format(nuc_zamid, row)

        # Diag terms
        diag_val = decay_mat[row][row]
        Ctxt += ' {} {},'.format(row, diag_val)

        # Non diag terms
        for col in range(N):
            if col == row:
                continue
            if decay_mat[row][col] != 0.0:
                Ctxt += ' {} {},'.format(col, decay_mat[row][col])

        Ctxt += '\n'

    return Ctxt
    # C = open(txt_mat_name, 'w')
    # C.write(Ctxt)


# Old methods that create matrixes from txt versions
# Need to decide if keep or not
# There are not used anymore within the code so far

def xs_mat_from_txt(Btxt_name):

    Btxt_rel_path = '/data/{}'.format(Btxt_name)
    dir_path = os.path.dirname(__file__)
    Btxt_path = dir_path + Btxt_rel_path
    Btxt = open(Btxt_path, 'r')
    B_line = Btxt.readlines()

    N = len(B_line)
    xs_mat = np.zeros((N,N))
    for line in B_line:
        data = line.split('|')[1]
        row = int(data.split(':')[0])
        col_data = data.split(':')[1].split(',')
        col_data.remove('\n') #The format of Ctxt and Btxt contains a last coma that needs to be removed
        for data in col_data:
            col = int(data.split()[0])
            val = float(data.split()[1])
            xs_mat[row][col] =  val

    return xs_mat

def decay_mat_from_txt(Ctxt_name):

    Ctxt_rel_path = '/data/{}'.format(Ctxt_name)
    dir_path = os.path.dirname(__file__)
    Ctxt_path = dir_path + Ctxt_rel_path
    Ctxt = open(Ctxt_path, 'r')
    C_line = Ctxt.readlines()

    N = len(C_line)
    decay_mat = np.zeros((N,N))
    for line in C_line:
        data = line.split('|')[1]
        row = int(data.split(':')[0])
        col_data = data.split(':')[1].split(',')
        col_data.remove('\n') #The format of Ctxt and Btxt contains a last coma that needs to be removed
        for data in col_data:
            col = int(data.split()[0])
            val = float(data.split()[1])
            decay_mat[row][col] =  val

    return decay_mat


def initial_vect_from_txt(mattxt_name, passdic):

    mattxt_rel_path = '/data/{}'.format(mattxt_name)
    dir_path = os.path.dirname(__file__)
    mattxt_path = dir_path + mattxt_rel_path
    mattxt = open(mattxt_path, 'r')
    mat_line = mattxt.readlines()

    N = len(mat_line)
    vect = np.zeros((N))
    i = 0
    for line in mat_line:
        zamid = line.split('|')[0]
        if zamid in passdic:
            vect[i] = passdic[zamid].current_dens
        i += 1

    return vect

def nucl_list_from_txt(mattxt_name):

    mattxt_rel_path = '/data/{}'.format(mattxt_name)
    dir_path = os.path.dirname(__file__)
    mattxt_path = dir_path + mattxt_rel_path
    mattxt = open(mattxt_path, 'r')
    mat_line = mattxt.readlines()

    N = len(mat_line)
    nucl_list = []
    for line in mat_line:
        zamid = line.split('|')[0]
        nucl_list.append(zamid)

    return nucl_list



def _get_initial_vect(passlist):

    passport_list = passlist.passport_list

    N = len(passport_list)
    
    vect = np.zeros((N))
    for i in range(N):
        nuc_pass = passport_list[i]
        vect[i] = nuc_pass.current_dens

    return vect
