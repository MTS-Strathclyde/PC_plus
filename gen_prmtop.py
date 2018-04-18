#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 17:48:40 2015

@author: max
"""

from __future__ import print_function, division

import os
import argparse
import sys
import glob
import shlex
try:
    import pybel
except ImportError:
    print("Pybel not found, some features may not work correctly!")
import subprocess
try:
    import parmed.tools as pact
    from parmed.amber import AmberParm
except ImportError as e:
    print('Make sure you are using AmberTools17!')
    raise e

MIN_SCRIPT = """Normal minimization
   &cntrl
        imin=1,      ! perform minimization
        maxcyc=10000,  ! The maximum number of cycles of minimization
        drms=1e-3,   ! RMS force
        ntmin=3,     ! xmin algorithm
        ntb=0,       ! no periodic boundary
        cut=999.,    ! non-bonded cutoff
        ntpr=5       ! printing frequency
   /
"""



RUNLEAP = """source {leaprc}
mol = loadmol2 "{name}.mol2"
check mol
loadamberparams "{name}.frcmod"
SaveAmberParm mol "{name}.prmtop" "{name}.incrd"
SavePdb mol "{name}.pdb"
quit
"""
#loadamberparams {name}.frcmod
def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description=""" Prepare AMBER prmtop file.""")
    #Positional args
    parser.add_argument('file', metavar='molec.pdb',
                        help="""Input file. Either structure file [pdb or mol2]
                        or existing prmtop file.""")
    #Optional args
    parser.add_argument('-f', '--moltype',
                        help="""Molecule type [pdb]""", default='pdb')
    parser.add_argument('-c', '--molcharge',
                        help="""Charge of the solute [0]""", default=0,
                        type=int)
    parser.add_argument('-cm', '--charge_model',
                        help="""Charge model to use. All antechamber options as
                        well as opls will work [bcc].""", default='bcc')                     
    parser.add_argument('-lj', '--lennard_jones',
                        help="""Lennard jones parameters to use (opls, gaff, gaff2, sybyl)
                        [gaff]""", default='gaff')                     
    parser.add_argument('--scale_r',
                        help="""Scale all lennard jones radiuses . [1.0]""",
                        default=1.0, type=float)
    parser.add_argument('--scale_eps',
                        help="""Scale all epsilon values. [1.0]""",
                        default=1.0, type=float)
    parser.add_argument('--multiplicity',
                        help="""Multiplicity of the solute [1]""", default=1,
                        type=int)
    parser.add_argument('--clean_up',
                        help=""" How should auxiliary files be treated:
                        0 - delete nothing;
                        1 - delete all [default];
                        """, default=1, type=int)
    parser.add_argument('--charge_f',
                        help="""Supply a charge file to the antechamber. A
                        file should contain a list of atomic partial charges
                        appearing in the same ordear as are atoms in tge pdb
                        file. One row shouldn't contain more than 8 charges.
                        """)
    parser.add_argument('--scale_chg',
                        help="""Scale all charges predicted by the
                        antechamber by a certain value . [1.0]""",
                        default=1.0, type=float)
    parser.add_argument('--input_chg',
                        help="""Manually input charge for every atom type""",
                        action='store_true')
    parser.add_argument('--input_lj',
                        help="""Manually input LJ parameters for every atom type""",
                        action='store_true')
    parser.add_argument('--lj_radius_type',
                        help="""Specify type of LJ radius to input (rmin/2 or 
                        sigma) [rmin/2]""",
                        default='rmin/2')
    parser.add_argument('-n', '--new_name',
                        help="""Name of the new prmtop. By default will overwrite
                        existing one.""")
    parser.add_argument('--nomod',
                        help=""" Use this to simply run antechamber (and 
                        minimization) and exit. Doesn't change default
                        atomtypes.""",
                        action='store_true')
    parser.add_argument('--minimize',
                        help=""" Minimize solute using gradient descent.""",
                        action='store_true')
    parser.add_argument('--ffld_out',
                        help=""" File containing ffld server output that script
                        will use to assign opls radii and charges. lj option 
                        must be set to opls for this to work!""")
    return parser.parse_args(argv)



def minimize_solute(name, prmtop_name, args):
    """ Minimize the solute structure using Sander. 
    The pdb file in the calculation directory **will** be overwritten.
    
    Parameters
    ----------

    name : string
        Calculation name
    
    prmtop_name : string, default None
        Topology file for calculation. Path should be given
        relative to the directory in which pdb file is located.
        
    args : Namespace
        Command line arguments

    """
    min_script_name = 'min_script.input'
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    print('Minimizing solute structure')
    # use or create restart (incrd) file
    rst_name = os.path.join(name + '.incrd')
    if not os.path.isfile(rst_name):
        subprocess.check_output(['antechamber',
                         '-i', '{}.pdb'.format(no_p_name),
                         '-fi', 'pdb',
                         '-o', '{}.incrd'.format(no_p_name), #output file
                         '-fo', 'rst',   #output format
                         ],
                         cwd=p)
    with open(os.path.join(p, min_script_name), 'w') as f:
        f.write(MIN_SCRIPT)
    # minimize solute and overwrite restart file
    print(subprocess.check_output(['sander',
                                       '-O', #overwrite files
                                       '-i', min_script_name,
                                       '-p', prmtop_name,
                                       '-c', '{}.incrd'.format(no_p_name),
                                       '-r', '{}.incrd'.format(no_p_name)],
                                       cwd=p))
    with open(rst_name) as f:
        rst_text = f.read()
    # converst restart file to pdb and write
    p = subprocess.Popen(['ambpdb',
                          '-p', prmtop_name], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    pdb_min = p.communicate(input=rst_text)[0]
    #print(pdb_min)
    with open(name + '.pdb', 'w') as f:
        f.write(pdb_min)



def generate_prmtop(name, args):
    """Generate topology file using GAFF charges scaled by the supplied
    factor.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    logfile : A writable file object
        A file to which calculation std. output will be written

    molcharge : int
        Charge of the solute

    multiplicity : int
        Multiplicty of the solute

    charge_f : string
        A name of the file containing atomic partial charges readable by
        antechamber. The file will be supplied to the antechamber through
        the option rc. It should contain partial charges appearing in the
        same order as in the pdb file with not more than 8 charges per row.

    Returns
    -------
    out: string
        Returns the name of prepared topology file
    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    #Firstly we use antechamber to recognize atom and bonding types, and
    #generate topology. Even if we want opls lj paramters, we still need
    #to assign atom types first to generate .prmtop file and only then
    #change them.
    if args.input_chg or args.charge_model == 'opls':
        subprocess.check_output(["antechamber",
                         '-i', args.file,
                         '-fi', args.moltype,
                         '-at', args.lennard_jones, 
                         '-o', '{}.mol2'.format(no_p_name), #output file
                         '-fo', 'mol2',   #output format describing each residue
                         '-s', '2',    #status info ; 2 means verbose
                         '-dr', 'n'   #do not check molecule for "correctness"
                         ],
                         cwd=p)
    elif args.charge_f:
        raise ValueError('This option is bugged - needs fixing')
        args.charge_f = os.path.relpath(args.charge_f, p) #path relative to calc dir
        subprocess.check_output(['antechamber',
                         '-i', args.file,
                         '-fi', args.moltype,
                         '-at', args.lennard_jones,
                         '-o', '{}.mol2'.format(no_p_name), #output file
                         '-fo', 'mol2',   #output format describing each residue
                         '-c', 'rc',      #charge method: read in charge
                         '-cf', args.charge_f, #file with charges
                         '-s', '2',    #status info ; 2 means verbose
                         ],
                         cwd=p)
    else:
        subprocess.check_output(['antechamber',
                      '-i',  args.file,
                      '-fi', args.moltype,
                      '-at', args.lennard_jones,
                      '-o', '{}.mol2'.format(no_p_name), #output file
                      '-fo', 'mol2',  #output format
                      '-c', args.charge_model,      #charge method 
                      '-s', '2',    #status info ; 2 means verbose
                      '-nc', str(args.molcharge),   #Net molecule charge
                      '-m', str(args.multiplicity)   #Multiplicity
                      ],
                      cwd=p)
#    #Run parmchk to generate missing gaff force field parameters
    try:
        subprocess.check_output(['parmchk2',
                         '-i', '{}.mol2'.format(no_p_name),
                         '-f', 'mol2',
                         '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                         cwd=p)
    except  subprocess.CalledProcessError:
        # try falling back on parmchk1
        subprocess.check_output(['parmchk',
                         '-i', '{}.mol2'.format(no_p_name),
                         '-f', 'mol2',
                         '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                         cwd=p)
    #Run tleap to generate topology and coordinates for the molecule
    if args.lennard_jones == 'gaff2':
        leaprc = 'leaprc.gaff2'
    else:
        leaprc = 'leaprc.gaff'
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'w') as f:
        f.write(RUNLEAP.format(name=no_p_name, leaprc=leaprc))
    subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    prmtop_name = '{}.prmtop'.format(no_p_name)
    return prmtop_name


def run_ffld_server(args, name):
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    if args.ffld_out:
        with open(args.ffld_out) as f:
            out = f.readlines()
    else:
        schrod_path = os.environ['SCHRODINGER']
        command = """{}/utilities/ffld_server -ipdb {}.pdb -print_parameters"""\
                    .format(schrod_path,name)
        out = subprocess.check_output(shlex.split(command)).splitlines()
    # check if ff parms were generated 
    parm_data_start = None
    for i, l in enumerate(out):
        if l.startswith('OPLSAA FORCE FIELD TYPE ASSIGNED'):
            parm_data_start = i + 4
    if not parm_data_start:
        # try to use mol2 for generation
        print('Failed to assigne parameters for pdb, trying mol2')
        # gen mol2 using babel - it recognizes atoms
        mol = pybel.readfile('pdb', name + '.pdb').next()
        mol_mol2 = mol.write('mol2')
        # in babel mol2 there is a problem atoms called CLx are regarded as
        # carbons, not as clorines
        fixed_mol2 = []
        for l in mol_mol2.splitlines():
            ls = l.split()
            if len(ls) == 9: # atom row
                if ls[1].startswith('CL') or ls[1].startswith('Cl'):
                    ls[5] = 'Cl'
                    l = '{:>7}  {:<7}{:>10}{:>10}{:>10} {:<8}{:<3}{:<8}{:>10}'.format(*ls)
                if ls[1].startswith('BR') or ls[1].startswith('Br'):
                    ls[5] = 'Br'
                    l = '{:>7}  {:<7}{:>10}{:>10}{:>10} {:<8}{:<3}{:<8}{:>10}'.format(*ls)                
            fixed_mol2.append(l)
        fixed_mol2 = '\n'.join(fixed_mol2)
        with open(name + '.mol2', 'w') as f:
            f.write(fixed_mol2)
        command = """{}/utilities/ffld_server -imol2 {}.mol2 -print_parameters"""\
                    .format(schrod_path,name)
        out = subprocess.check_output(shlex.split(command)).splitlines()
        # check again
        for i, l in enumerate(out):
            if l.startswith('OPLSAA FORCE FIELD TYPE ASSIGNED'):
                parm_data_start = i + 4
        if not parm_data_start:
            raise ValueError('Failed to assign oplsaa parameters to {}'.format(name))
    return out, parm_data_start


def get_opls_parameters(args, name):
    out, parm_data_start = run_ffld_server(args, name)
    radii = []
    epss = []
    chgs = []    
    for l in out[parm_data_start:]:
        if not l.startswith('-----'):
            l = l.split()
            radii.append(float(l[5])/2*2**(1./6))   # rmin/2, A
            epss.append(float(l[6]))                # kcal/mol 
            chgs.append(float(l[4]))                # e
        else:
            break
    return radii, epss, chgs
    

def get_usr_input(parm_name, atname, old_value):
    usr_value = raw_input('Provide new {} for {} (blank=keep old) [{:f}]: '.\
                        format(parm_name, atname, old_value))
    if usr_value:
        return float(usr_value)
    else:
        return old_value    
    
    
def get_chargef_charges(charge_f):
    with open(charge_f) as f:
        charges = f.read().split()
    return charges
    

def prepare_prmtop(args, name):
    """ Places appropriate prmtop file into the calculation folder and scales
    it's charges if necessary.

    Parameters
    ----------

    args : Namespace
        Command line arguments

    Returns
    -------

    out : string
        Path to prmtop file.

    """
    # The atom specification in ParmedActions is terrible and requires the use
    # of masks. A good description of what it is can be found in Amber14 manual
    # section 28.2.3
    if args.file.endswith('.prmtop'):
        prmtop_name = args.file
    else:
        prmtop_name = generate_prmtop(name, args)
    if args.new_name:
        new_name = args.new_name
    else:
        new_name = prmtop_name
    if args.minimize:
        minimize_solute(name, prmtop_name, args)
    if args.nomod:
        # we are done
        return 0
    parm = AmberParm(prmtop_name)
    if args.charge_model == 'opls' or args.lennard_jones == 'opls':
        opls_radii, opls_epss, opls_chgs = get_opls_parameters(args, name)
    # account for scenario when charge_f is submitted along with existing prmtop
    if args.charge_f and args.file.endswith('.prmtop'):
        chargef_charges = get_chargef_charges(args.charge_f)
    #iterate over atoms
    for i, atom in enumerate(parm.atoms):
        attyp, atname, attchg = atom.type, atom.name, float(atom.charge)
        #print(attchg)
        nbidx = parm.LJ_types[attyp]
        lj_r = float(parm.LJ_radius[nbidx - 1])
        lj_eps = float(parm.LJ_depth[nbidx - 1])
        # change attyp to atnmae
        act = pact.change(parm, '@{} AMBER_ATOM_TYPE {}'.format(atname, atname))
        act.execute()
        # deal with chgs
        if args.input_chg:
            print()
            attchg = get_usr_input('charge', atname, attchg)
        elif args.charge_model == 'opls':
            attchg = opls_chgs[i]
        elif args.charge_f and args.file.endswith('.prmtop'):
            attchg = float(chargef_charges[i])
        attchg = attchg * args.scale_chg
        act = pact.change(parm, '@{} charge {:f}'.format(atname, float(attchg)))
        act.execute()
        # deal with lj
        if args.input_lj:
            if args.lj_radius_type == 'sigma':
                lj_r = lj_r*2./(2**(1./6))  # convert to sigma
            lj_r = get_usr_input('lj_r', atname, lj_r)
            if args.lj_radius_type == 'sigma':
                lj_r = lj_r/2.*(2**(1./6))
            lj_eps = get_usr_input('lj_eps', atname, lj_eps)
        elif args.lennard_jones== 'opls':
            lj_r = opls_radii[i]
            lj_eps = opls_epss[i]
        lj_r = lj_r * args.scale_r
        lj_eps = lj_eps * args.scale_eps
        #print(lj_r, lj_eps)
        act = pact.changeLJSingleType(parm, '@{} {:f} {:f}'.format(atname, lj_r, lj_eps))
        act.execute()
    #parm.overwrite = True
    parm.write_parm(new_name)


def main(argv):
    args = process_command_line(argv)
    name = os.path.splitext(args.file)[0]
    prepare_prmtop(args, name)
    #clean
    level = args.clean_up
    p, no_p_name = os.path.split(name)
    to_del1_glob = ['ANTECHAMBER*', 'sqm*']
    to_del1_files = [ no_p_name + '.frcmod',
                     'ATOMTYPE.INF',
                    'leap.log', 'NEWPDB.PDB', 'PREP.INF', 'runleap.in']
    will_be_deleted_list = []
    if level == 1:
        for wildcard in to_del1_glob:
            will_be_deleted_list.extend(glob.glob(os.path.join(p, wildcard)))
        will_be_deleted_list.extend([os.path.join(p, f) for f in \
                                                            to_del1_files])
    for f in will_be_deleted_list:
        try:
            os.unlink(f)
        except OSError as e:
            if e.errno == 2:
                pass
            else:
                raise e


if __name__ == '__main__':
    main(sys.argv[1:])
    

