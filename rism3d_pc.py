#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 16:04:49 2014

@author: Maksim Misin (mishin1991@gmail.com)

Compute solvation free energy with rism3d.singlpnt and apply PC corrections. 
The script also can prepare topology, minimize solute and run generate 
susceptibility files using 1D-RISM. The output is written to separate .log and 
results.txt files.

As an input takes pdb file compatible with antechamber. For example:
ATOM      1  C1  MOL     1       3.537   1.423   0.000  1.00  0.00
ATOM      2  H1  MOL     1       4.089   2.224   0.496  1.00  0.00
ATOM      3  H2  MOL     1       4.222   0.611  -0.254  1.00  0.00
ATOM      4  H3  MOL     1       2.759   1.049   0.669  1.00  0.00
ATOM      5  H4  MOL     1       3.077   1.810  -0.912  1.00  0.00
TER
END

To run the simmulation simply type:
python rism3d_pc.py molecule.pdb

The script requires working installations of python2.7 and AmberTools 12+

For more information run:
python rism3d_isc.py -h

If you find the script useful, please cite:

Misin, M.; Fedorov, M.; Palmer, D. Hydration Free Energies of Ionic Species 
by Molecular Theory and Simulation, J. Phys. Chem. B. 2016. 
http://dx.doi.org/10.1021/acs.jpcb.5b10809

Misin, M.; Fedorov, M. V.; Palmer, D. S. Accurate Hydration Free Energies at 
a Wide Range of Temperatures from 3D-RISM. J. Chem. Phys. 2015, 142, 091105. 
http://dx.doi.org/10.1063/1.4914315

an article where PC+ correction was originally proposed:

Sergiievskyi, V.; Jeanmairet, G.; Levesque, M.; Borgis, D. Solvation 
Free-Energy Pressure Corrections in the Three Dimensional Reference Interaction
Site Model. J. Chem. Phys. 2015, 143, 184116. 
http://dx.doi.org/10.1063/1.4935065

as well all related amber RISM programs.



    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import sys
import argparse
import subprocess
import distutils.spawn
import shutil
import os
import glob
import datetime
import time
import threading
import signal
import numpy as np


## Non RISM globals ##
__version__ = '2016.1'

REQUIRED_EXECUTABLES = ['antechamber', 'parmchk', 'tleap', 'rism3d.snglpnt',
                        'rism1d', 'sander', 'ambpdb']

IS_UNIX = os.name == 'posix'

## Constants ##
K_B = 1.9872041E-3 # boltzmann const in kcal/mol/K
N_A = 6.022141e23 # avogadro's constant


## RISM-related globals ##

SOLV_SUCEPT_SCRPT = """#!/bin/csh -f

cat > {name1d}.inp <<EOF
&PARAMETERS
    THEORY='{rism1d}', CLOSUR='{closure}',           !Theory
    NR=16384, DR=0.025,                    !Grid
    OUTLST='xCGT', rout=0,                 !Output
    NIS=20, DELVV=0.3, TOLVV=1.e-12,       !MDIIS
    KSAVE=-1, KSHOW=1, maxstep=10000,      !Check pointing and iterations
    SMEAR=1, ADBCOR=0.5,                   !Electrostatics
    TEMPER={temp}, DIEps={diel},           !bulk solvent properties
    NSP=1
/
    &SPECIES                               !SPC water
    DENSITY={conc}d0,
    MODEL="$AMBERHOME/dat/rism1d/model/{smodel}.mdl"
/
EOF

rism1d {name1d} > {name1d}.out || goto error

"""

RUNLEAP = """source leaprc.gaff
mol = loadmol2 "{name}.mol2"
check mol
loadamberparams "{name}.frcmod"
SaveAmberParm mol "{name}.prmtop" "{name}.incrd"
SavePdb MOL "{name}.pdb"
quit
"""

MIN_SCRIPT = """Normal minimization
   &cntrl
        imin=1,      ! perform minimization
        maxcyc=200,  ! The maximum number of cycles of minimization
        drms=1e-3,   ! RMS force
        ntmin=3,     ! xmin algorithm
        ntb=0,       ! no periodic boundary
        cut=999.,    ! non-bonded cutoff
        ntpr=5       ! printing frequency
   /
"""
MIN_SCRIPT_RISM = """Minimization with rism
   &cntrl
	imin=1, maxcyc=200, drms=1e-3, ntmin=3,ntb=0,cut=999.,
	ntpr=5, irism=1
   /
   &rism
	closure='{closure}', buffer=25, tolerance=1e-4,solvcut=9999
   /
"""


MIN_SCRIPT_NAME = 'min.input'
RISM1D_NAME = '{smodel}_{temp}'
RESULTS_NAME = 'results.txt'

RESULTS = """dGsolv(closure)= {exchem} kcal/mol
dGsolv(GF)= {gf} kcal/mol
PMV= {pmv} A^3

dGsolv(PC+)= {PC_plus} kcal/mol
dGsolv(PC)= {PC} kcal/mol

P_minus_ideal_gas_pressure= {pressure_plus} kcal/mol/A^3
P= {pressure} kcal/mol/A^3

"""


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
    parser = argparse.ArgumentParser(description="""Run 3D-RISM single point
            calculation. This script is a wrapper around Amber rism3d.snglpnt program,
            designed to simplify calculations. It also computes pressure corrections
            to RISM solvation free energy.
            """)
    #molecule options
    molecule_options = parser.add_argument_group('Solute options',
                            """Options related to a solute molecule. 
                            Calculation requires only pdb and prmtop files.
                            If prmtop file is not present, script will try to 
                            create this file using antechamber and
                            tleap. By default AM1-BCC charges and GAFF
                            will be used.""")
    molecule_options.add_argument('file', metavar='molec.pdb',
                        help="""Input solute file. Must be in pdb format
                        acceptable by Antechamber. Must have a .pdb
                        extension.""")
    molecule_options.add_argument('-p', '--prmtop',
                        help="""Path to parameters and topology (prmtop) file
                        of solute.""")
    molecule_options.add_argument('--scale_chg',
                        help="""Scale all solute by this value prmtop file [1.0].""",
                        default=1.0, type=float)
    molecule_options.add_argument('-c', '--molcharge',
                        help="""Charge of the solute [0].""", default=0,
                        type=int)
    molecule_options.add_argument('--multiplicity',
                        help="""Multiplicity of the solute [1].""", default=1,
                        type=int)
    molecule_options.add_argument('--minimize',
                        help=""" Minimize solute before performing
                        3D-RISM calculation using either gradient descent (min) or
                        RISM minimization using sander (rism, not recommended).
                        If no keywords are provided minimization is not performed.
                        """)
    #1drism
    rism1d_options = parser.add_argument_group('1D-RISM options',
                            """3D-RISM calculation requires xvv file 
                            (site-site susceptibilities) to run.
                            If such file is not provided script can run a 
                            1D-RISM to try to genrate it, but will asume that
                            the solvent is pure water. Some of the 1D-RISM related
                            options are in 3D-RISM group. """)
    rism1d_options.add_argument('-x', '--xvv',
                        help="""Path to an existing xvv file. This will skip 1D-RISM
                        calculation and all related parameters
                        will be ignored. Solvent density as well as calculation
                        temeprature will be read from this file.""")
    rism1d_options.add_argument('--smodel',
                        help="""Solvent model for 1D-RISM calculation available in
                        "$AMBERHOME/dat/rism1d/model/{smodel}.mdl" [SPC].""",
                        default="SPC")
    rism1d_options.add_argument('--rism1d',
                        help="""Type of 1D-RISM theory. Only DRISM has been
                        extensively tested [DRISM].""",
                        default="DRISM")
    #3drism
    rism3d_options = parser.add_argument_group('3D-RISM options',
                            """Options related to the main calculation.""")
    rism3d_options.add_argument('--closure',
                        help="""Brdige closure for 3D-RISM and 1D-RISM calculation
                        if it is necessary. Either HNC, KH, PSEn
                        (n in PSEn should be an integer) or a list of them
                        for sequential convergence. For 1D-RISM calculation
                        only last closure will be used [PSE3].""",
                        default=["PSE3"], nargs='*')
    rism3d_options.add_argument('-t', '--temperature',
                        help="""Temperature in K at which calculation will be
                        preformed. If xvv file was provided, this option will be
                        used only for naming directory [298.15].""",
                        default=298.15, type=float)
    rism3d_options.add_argument('--clean_up',
                        help=""" How should auxiliary files be treated:
                        0 - delete nothing;
                        1 - delete some [default];
                        2 - delete all but input, results, and log.
                        """, default=1, type=int)
    rism3d_options.add_argument('--dir_name',
                        help="""Custom name for produced calculation
                        directory. The default one is:
                        {mol_name}_{temperature}.""")
    rism3d_options.add_argument('--timeout',
                        help=""" Minutes after which 3D-RISM calculation
                        will be killed. Use 0 for no timeout. [0]. Only works
                        on Unix-like system.
                        """, default=0, type=float)
    rism3d_options.add_argument('--tolerance',
                        help=""" Maximum residual values for 3D-RISM solution
                        convergence. If many closures a list of closures
                        can be supplied [1E-5].""",
                        default=['1e-5'], nargs='*')
    rism3d_options.add_argument('--write_g',
                        help="""Write radial distribution functions produced
                        in 3D-RISM calculation.""",
                        action='store_true')
    rism3d_options.add_argument('--write_c',
                        help="""Write direct correlation function produced
                        in 3D-RISM calculation.""",
                        action='store_true')
    rism3d_options.add_argument('--write_u',
                        help="""Write solute solvent potential energy grid.""",
                        action='store_true')
    rism3d_options.add_argument('--write_h',
                        help="""Write total correlation function in k space.""",
                        action='store_true')
    rism3d_options.add_argument('--write_asymp',
                        help="""Write asymptotics of total and direct 
                        correlation fuctions in real space.""",
                        action='store_true')
    rism3d_options.add_argument('--noasympcorr',
                        help=""" Thermodynamics of 3D-RISM is calculated without
                        long-range asymptotics.""",
                        action='store_true')
    rism3d_options.add_argument('--buffer',
                        help="""Minimum distance between the solute and the
                        edge of the solvent box in A for 3D-RISM
                        calculation [25].""",
                        default=25, type=float)
    rism3d_options.add_argument('--solvbox',
                        help="""Size of the x, y, and z dimensions of the box in
                        Angstroms. Specifying this parameter overrides buffer.""",
                        nargs=3)
    rism3d_options.add_argument('--grdsp',
                        help="""Linear grid spacings for x, y and z
                        dimensions. Should be separated with spaces. Units: A
                        [0.5 0.5 0.5].""",
                        default=(.5, .5, .5), nargs=3)
    rism3d_options.add_argument('--polar_decomp',
                        help="""Decomposes solvation free energy into polar
                        and non-polar components.""",
                        action='store_true')
    rism3d_options.add_argument('--verbose3d',
                        help="""Verbosity of 3D-RISM calculation. 0 - print
                        nothing; 1 - print iterations; 2 - print all [2].""",
                        default=2, type=int)
    rism3d_options.add_argument('--maxstep3d',
                        help="""Maximum number of iterations in 3D-RISM
                        calculation [500].""",
                        default=500, type=int)
    rism3d_options.add_argument('--rism3d_path',
                        help="""Specify absolute path or exact name of rism3d.sngpnt
                        [rism3d.snglpnt].""", default='rism3d.snglpnt')
    return parser.parse_args(argv)


def water_dielectric_const(T):
    """Return water dielectric constant for temperature 253.15K < T < 383.15K.
    Uses interpolation equation (eq. 9) for static dielectri constant found in
    the doucment by The International Association for the Properties of
    Water and Steam from 2011
    <http://www.iapws.org/relguide/LiquidWater.pdf>`__.
    Pressure = 0.1 MPa

    Parameters
    ----------
    T : float
        Temperature in K

    Returns
    -------
    e : float
        Water dielectric constant at T

    Examples
    --------
    >>> round(water_dielectric_const(273.15), 3)
    87.927
    >>> round(water_dielectric_const(298.15), 3)
    78.375
    >>> round(water_dielectric_const(375), 3)
    55.266
    """
    if not 253.15 <= T <= 383.15:
        raise ValueError("Temperature is outside of allowed range.")
    T_star = T/300.0
    coefs = [-43.7527, 299.504, -399.364, 221.327]
    exp_f = [-0.05, -1.47, -2.11, -2.31]
    e = 0
    for i in range(4):
        e += coefs[i]*T_star**(exp_f[i])
    return e


def water_concentration(T):
    """Return water concentration for temperature range 253.15K < T < 383.15K.
    Uses interpolation equation (eq. 2) for specific volume found in
    the doucment by The International Association for the Properties of
    Water and Steam from 2011
    <http://www.iapws.org/relguide/LiquidWater.pdf>`__.
    Pressure = 0.1 MPa

    Parameters
    ----------
    T : float
        Temperature in K

    Returns
    -------
    conc : float
        Water conentration at T in mol/l

    Examples
    --------
    >>> round(water_concentration(273.15), 3)
    55.498
    >>> round(water_concentration(298.15), 3)
    55.343
    """
    if not 253.15 <= T <= 383.15:
        raise ValueError("Temperature is outside of allowed range.")
    p0 = 10.0**5    # Pa
    R = 8.31464     # J/mol/K
    Tr = 10.0
    Ta = 593.0
    Tb = 232.0
    a = [1.93763157E-2,
         6.74458446E+3,
        -2.22521604E+5,
         1.00231247E+8,
        -1.63552118E+9,
         8.32299658E+9]
    b = [5.78545292E-3,
        -1.53195665E-2,
         3.11337859E-2,
        -4.23546241E-2,
         3.38713507E-2,
        -1.19946761E-2]
    n = [None, 4., 5., 7., 8., 9.]
    m = [1., 2., 3., 4., 5., 6.]
    def alpha(T):
        return Tr/(Ta - T)
    def beta(T):
        return Tr/(T - Tb)
    coef = a[0] + b[0]*beta(T)**m[0]
    for i in range(1, 6):
        coef += a[i]*alpha(T)**n[i] + b[i]*beta(T)**m[i]
    v0 = R*Tr/p0*coef  # m3/mol
    return 1/(v0*1000)    # mol/L


class Xvv(object):
    """ Wrapper around xvvfile used to compute 3d-rism pressure """
    def __init__(self, fname):
        """ Read xvvfile and set instance attributes 
        
        Parameters
        ----------
        
        fname : string
            Path to a valid xvv file
        """
        self.fname = fname
        self.ngrid = None
        self.nsites = None
        self.nspecies = None
        self.temperature = None
        self.dr = None
        self.atom_names = None
        self.densities = None
        self.xvv_data = None
        self.multiplicities = None
        self.unique_sites_per_species = None
        self.total_sites_per_species = None
        self.species_densities = None
        self.normalized_densities = None
        self._read_xvvfile()
        self._compute_species_properties()

    def _read_xvvfile(self):
        with open(self.fname) as f:
            lines = f.readlines()
        tot_lines = len(lines)
        for i, line in enumerate(lines):
            line = line.split()
            if len(line) <= 1:
                continue
            if line[1] == 'POINTERS':
                data = map(int, lines[i+2].split())
                self.ngrid, self.nsites, self.nspecies = data
            if line[1] == 'MTV':
                self.multiplicities = map(int, lines[i+2].split())
            if line[1] == 'NVSP':
                self.unique_sites_per_species = map(int, lines[i+2].split())
            if line[1] == 'THERMO':
                data = lines[i+2].split()
                self.temperature = float(data[0]) # K
                self.dr = float(data[4]) # Angstrom
            if line[1] == 'ATOM_NAME':
                data = lines[i+2].strip()
                #split into groups of 4
                self.atom_names = [data[i:i+4].strip() for i in range(0, len(data), 4)]
            if line[1] == 'RHOV' and len(line) == 2:
                self.densities = map(float, lines[i+2].split())
                #are there more lines with density?
                counter = 3
                while lines[i+counter].startswith(' '):
                    self.densities.extend(map(float, lines[i+counter].split()))
                    counter += 1
                try:
                    assert len(self.densities) == len(self.atom_names)
                except AssertionError:
                    print('Inconsistent number of densities and atom names')
                    print(self.densities)
                    print(self.atom_names)
                    raise ValueError
            if line[1] == 'XVV' and len(line) == 2:
                self.xvv_data = []
                xvv_ind = i + 2
                while xvv_ind < tot_lines and not lines[xvv_ind].startswith('%'):
                    self.xvv_data.extend(lines[xvv_ind].split())
                    xvv_ind += 1
                break
        assert len(self.xvv_data) == self.ngrid*self.nsites*self.nsites
        self.xvv_data = np.array(self.xvv_data, dtype=float)
        self.xvv_data = np.reshape(self.xvv_data,
                                   (self.ngrid, self.nsites, self.nsites),
                                   order='F')

    def _compute_species_properties(self):
        self.normalized_densities = []
        for density, multiplicity in zip(self.densities, self.multiplicities):
            self.normalized_densities.append(density/multiplicity)
        self.species_densities = []
        self.total_sites_per_species = []
        pointer = 0 
        for sp_sites in self.unique_sites_per_species:
            pointer += sp_sites
            total_sites = sum(self.multiplicities[pointer - sp_sites:pointer])
            self.total_sites_per_species.append(total_sites)
            self.species_densities.append(self.normalized_densities[pointer - 1])
        assert len(self.species_densities) == self.nspecies
    
    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. The pressure can be pretty
            sensitive to it. It is recommended to experiment with a couple of
            k values or better, plot dependency of pressure on it to see
            which value works best.
            
        Return
        ------
        pressures : tuple of floats
            Tuple containeing two pressures.
            First element is 3D-RISM pressure (used in PC), second element is
            3D-RISM pressure minus ideal gas pressure (used in PC+).
            Both have units of kcal/mol/A^3.
        """
        xvv_k = self.xvv_data[k,:,:]
        density_vec = np.array(self.normalized_densities)
        mult_vec = np.array(self.multiplicities)
        # Z_k from sergievskyi's article
        z_k = mult_vec/density_vec*(np.identity(self.nsites) - np.linalg.inv(xvv_k))
        z_k_sum_densities2 = np.sum(density_vec*z_k*density_vec.T)
        densities_times_sites = [sites*dens for sites, dens in zip(self.total_sites_per_species,
                                                                   self.species_densities)]
        pressure = sum(densities_times_sites) - .5*z_k_sum_densities2
        pressure = pressure*self.temperature*K_B
        ideal_pressure  = sum(self.species_densities)*K_B*self.temperature
        return pressure, pressure - ideal_pressure


class RunCmd(threading.Thread):
    """ Will only work on Unix-like systems. And sometimes it will not work
    even there. """
    def __init__(self, cmd, timeout, logfile, cwd='.'):
        """ Run subprocess for a fixed ammount of time and then kill it. Writes
        both stdout and stderr to logfile as the process continues.

        Parameters
        ----------
        cmd : list
            Command to execute.

        timeout : float
            Time in minutes after which process will be killed. Pass 0 to
            let process run indefinitely.

        logfile : A writable file object
            A file to which calculation stdout and stderr will be written

        cwd : string, default .
            Directory in which process will be run
        """
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.logfile = logfile
        self.cwd = cwd
        if timeout == 0:
            self.timeout = None
        else:
            self.timeout = timeout*60 #convert to seconds

    def run(self):
        self.p = subprocess.Popen(self.cmd, preexec_fn=os.setsid,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT,
                           cwd=self.cwd)
        lines_iterator = iter(self.p.stdout.readline, "") #terminates when readline returns ""
        for line in lines_iterator:
            self.logfile.write(line)
            self.logfile.flush()

    def run_and_timeout(self):
        """Run subprocess

        Returns
        -------
        out : int
            Exit code of subprocess
        """
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            os.killpg(self.p.pid, signal.SIGTERM)
            self.join()
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print(self.cmd)
            print('Has run out of time')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            raise RuntimeError("3D-RISM calc didn't finish in time.")


def prepare_calc_directory(mol_path, T, dir_name=None):
    """Copy pdb file into the directory with the same name. If such directory
    doesn't exist it will try to create it.

    Parameters
    ----------
    mol_path : string
        Path to solute

    T : float
        A calculation temperature

    dir_name : string or None
        A name of directory in which calculation is produced. If None is
        supplied a default value will be used.

    Returns
    -------
    name: string
        Full path to pdb file without extension
    dir_name: string
        Name of the calculation directory
    """
    pdb_path, name_without_path = os.path.split(mol_path)
    if not dir_name:
        dir_name = os.path.join(pdb_path, name_without_path[:-4] + '_' + str(T))
    try:
        os.mkdir(dir_name)
    except OSError, e:
        if e.errno == 17:
            pass
        else:
            raise e
    name = os.path.join(dir_name, name_without_path)
    shutil.copy(mol_path, name)
    return name[:-4], dir_name


def prepare_logfile(name, argv):
    """Create a logfile which will be used throught the calculation.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension
    
    argv : list
        Command used to start script

    Returns
    -------
    out: file object
        Returns logfile.
    """
    p, _ = os.path.split(name)
    if p == '':
        p = '.'
    log_name = '{}.log'.format(name)
    logfile = open(log_name, 'w')
    logfile.write(str(datetime.datetime.now())+'\n')     # timestamp
    logfile.write(' '.join(argv) + '\n')
    logfile.flush()
    return logfile


def generate_prmtop(name, logfile, molcharge=0, multiplicity=1):
    """Generate topology file using GAFF and AM1-BCC charges, scaled by the supplied
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

    Returns
    -------
    out: string
        Returns the name of prepared topology file
    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    #Firstly we use antechamber to recognize atom and bonding types, and
    #generate topology
    ante_out = subprocess.check_output(['antechamber',
                     '-i', '{}.pdb'.format(no_p_name),
                     '-fi', 'pdb',
                     '-o', '{}.mol2'.format(no_p_name), #output file
                     '-fo', 'mol2',   #output format describing each residue
                     '-c', 'bcc',      #charge method  (AM1-BCC)
                     '-s', '2',    #status info ; 2 means verbose
                     '-nc', str(molcharge),   #Net molecule charge
                     '-m', str(multiplicity)   #Multiplicity
                     ],
                     cwd=p)
    logfile.write(ante_out)
    #Run parmchk to generate missing gaff force field parameters
    parm_out = subprocess.check_output(['parmchk2',
                     '-i', '{}.mol2'.format(no_p_name),
                     '-f', 'mol2',
                     '-o', '{}.frcmod'.format(no_p_name)], #file with missing FF params
                     cwd=p)
    logfile.write(parm_out)
    logfile.flush()
    #Run tleap to generate topology and coordinates for the molecule
    leap_input_name = os.path.join(p, 'runleap.in')
    with open(leap_input_name, 'w') as f:
        f.write(RUNLEAP.format(name=no_p_name))
    leap_out = subprocess.check_output(['tleap', '-f', 'runleap.in'], cwd=p)
    logfile.write(leap_out)
    logfile.flush()
    prmtop_name = '{}.prmtop'.format(no_p_name)
    return prmtop_name


def check_consistency(prmtop_name, name):
    """ Check if the ordering of atoms is the same both in pdb file and in
    prmtop file.
    
    Parameters
    ----------
    
    prmtop_name : string
        Path to prmtop file
    
    name : string
        Calculation name
    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    pdb_atom_list = []
    with open(name + '.pdb') as f:
        for line in f:
            if line.startswith('ATOM'):
                pdb_atom_list.append(line[12:16])  # atom name field
    with open(os.path.join(p, prmtop_name)) as f:
        prmtop_atom_string = ''
        n=4 # number of columns for each atname in prmtop
        lastcol=80 # last column in each prmtop line
        for line in f:
            if line.startswith('%FLAG ATOM_NAME'):
                f.next() # skip one line
                atom_name_row = next(f).strip('\n')
                while not atom_name_row.startswith('%'):
                    prmtop_atom_string += atom_name_row[:lastcol+1]
                    atom_name_row = next(f).strip('\n')
        # split string into characters of length n
        prmtop_atom_list = [prmtop_atom_string[i:i+n] \
                                for i in range(0, len(prmtop_atom_string), n)]
    for i, (pdb_aname, prmtop_aname) in \
                             enumerate(zip(pdb_atom_list, prmtop_atom_list)):
        if pdb_aname.strip() != prmtop_aname.strip():
            raise ValueError("The name of atom number {} in pdb file is: {} and \
in prmtop file: {}. Check the consistency between two files.".format(i,
                        pdb_aname.strip(), prmtop_aname.strip()))


def prepare_prmtop(args, name, dir_name, logfile):
    """ Places appropriate prmtop file into the calculation folder and scales
    it's charges if necessary.

    Parameters
    ----------

    args : Namespace
        Command line arguments

    name : string
        Calculation name

    dir_name : string
        Name of calculation directory

    logfile : File_object
        Calculation log

    Returns
    -------

    out : string
        Path to prmtop file.

    """
    # Copy prmtop file, because we might want to change it (change charges)
    if not args.prmtop:
        print('Running AM1-BCC calculation...')
        prmtop_name = generate_prmtop(name, logfile, args.molcharge, args.multiplicity)
    else:
        print('Reading user provided prmtop file')
        try:
            shutil.copy(args.prmtop, dir_name)
            prmtop_name = os.path.split(args.prmtop)[1]  # names are relative to calc directory
        except shutil.Error:
            # most likely error is due to src = destination
            prmtop_name = os.path.split(args.prmtop)[1]
        check_consistency(prmtop_name, name)
    #Open file and scale all the charges
    prm_lines = []
    with open(os.path.join(dir_name, prmtop_name), 'r') as f:
        for line in f:
            if line.startswith('%FLAG CHARGE'):
                prm_lines.append(line)
                prm_lines.append(next(f)) # skip format line
                next_line = next(f)
                while next_line.startswith(' '):
                    chrgs = next_line.split()
                    new_chrgs = ['{: .8E}'.format(float(chg)*args.scale_chg) for chg in chrgs]
                    prm_lines.append( ' ' + ' '.join(new_chrgs) + '\n')
                    next_line = next(f)
                prm_lines.append(next_line)
            else:
                prm_lines.append(line)
    with open(os.path.join(dir_name, prmtop_name), 'w') as f:
        f.writelines(prm_lines)
    return prmtop_name
    
    
def minimize_solute(name, logfile, prmtop_name, args, xvv):
    """ Minimize the solute structure using Sander. 
    The pdb file in the calculation directory **will** be overwritten.
    
    Parameters
    ----------

    name : string
        Calculation name

    logfile : File_object
        Calculation log
    
    prmtop_name : string, default None
        Topology file for calculation. Path should be given
        relative to the directory in which pdb file is located.
        
    args : Namespace
        Command line arguments

    xvv: string
        Name of an existing xvv file.

    """
    p, no_p_name = os.path.split(name)
    if p == '':
        p = '.'
    if args.minimize == 'min':
        min_script = MIN_SCRIPT
    elif args.minimize == 'rism':
        min_script = MIN_SCRIPT_RISM
    else:
        raise ValueError('Unknown minimization type. Use either rism or min.')
    print('Minimizing solute structure')
    # use or create restart (incrd) file
    rst_name = os.path.join(name + '.incrd')
    if not os.path.isfile(rst_name):
        conv_out = subprocess.check_output(['antechamber',
                         '-i', '{}.pdb'.format(no_p_name),
                         '-fi', 'pdb',
                         '-o', '{}.incrd'.format(no_p_name), #output file
                         '-fo', 'rst',   #output format
                         ],
                         cwd=p)
        logfile.write(conv_out)
    with open(os.path.join(p, MIN_SCRIPT_NAME), 'w') as f:
        f.write(min_script.format(closure=args.closure[-1]))
    # minimize solute and overwrite restart file
    min_out = subprocess.check_output(['sander',
                                       '-O', #overwrite files
                                       '-i', MIN_SCRIPT_NAME,
                                       '-p', prmtop_name,
                                       '-c', '{}.incrd'.format(no_p_name),
                                       '-r', '{}.incrd'.format(no_p_name),
                                       '-xvv', os.path.relpath(xvv, p)
                                       ],
                                       cwd=p)
    logfile.write(min_out)
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


def run_rism1d(name, logfile, T=298.15, smodel="SPC", rism1d="DRISM",
                        closure="PSE3", xvv=None):
    """Generate xvv file at a given temperature. 

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    T : float, default 298.15
        A calculation temperature

    smodel : string, default SPC
        Water model available in $AMBERHOME/dat/rism1d/model/{smodel}.mdl

    rism1d : string, default DRISM
        Type of 1D-RISM theory. Only DRISM has been extensively tested

    closure : string, default HNC
        Brdige closure which will be used in both 1D-RISM simmulation

    xvv : string
        Path to an existing xvv file. If supplied, 1D-RISM calulation will
        be skipped. Otherwise, returned path will be changed relative
        to the 3D-RISM calculation directory.

    Returns
    -------
    xvv : string
        Path to a created (existing) xvv file relative to calculation dir.
    """
    p, _ = os.path.split(name)
    rism1d_name = None
    if p == '':
        p = '.'
    if not xvv:
        #Generate solvent susceptibility file
        print('Running 1D-RISM calculation...')
        rism1d_name = RISM1D_NAME.format(smodel=smodel, temp=T)
        xvv_script_name_no_p = '{}.sh'.format(rism1d_name)
        xvv_script_name = os.path.join(p, xvv_script_name_no_p)
        diel = round(water_dielectric_const(T), 3)
        conc = round(water_concentration(T), 3)
        succ_srcirpt = SOLV_SUCEPT_SCRPT.format(temp=T, diel=diel, conc=conc,
                                                smodel=smodel, rism1d=rism1d,
                                                closure=closure,
                                                name1d=rism1d_name)
        with open(xvv_script_name, 'w') as f:
            f.write(succ_srcirpt)
        xvv_out = subprocess.check_output(['bash', xvv_script_name_no_p],
                                          cwd=p)
        logfile.write(xvv_out)
        logfile.flush()
        xvv = '{}.xvv'.format(rism1d_name)
    else:
        if abs(T - 298.15) > 1.0e-4:
            print('Warning: xvv file submitted')
            print('Temperature passed from command line will be ignored!')
        print('Using provided xvv file...')
        xvv = os.path.relpath(xvv, p) # path relative to the calc. directory
    return xvv


class RISM3D_Singlpnt(object):
    """ A class used to assist setting up 3D-RISM calculation.

    Init is used to specify temperature as well as non-standard names for
    topology (prmtop) or water susceptibility (xvv) files.

    The calculation details like closure or tolerance are defined in
    setup_calculation method.
    """
    def __init__(self, name, T, logfile, prmtop_name=None, xvv=None):
        """ Create a class for running rism3d.snglpnt.

        Parameters
        ----------
        name : string
            Full path to pdb file without extension

        T : float
            A calculation temperature

        logfile : A writable file object
            A file to which calculation std. output will be written

        prmtop_name : string, default None
            A name of topology file for calculation. If
            it is not specified defaults to name.prmtop. Path should be given
            relative to the directory in which pdb file is located.

        xvv : string, default None
            A name of susceptibility file for this calculation. Defaults
            to water_T.xvv, where T is calculation temperature rounded to
            two digits. Path should be given relative to the directory in
            which pdb file is located.

        """
        self.name = name
        self.T = T
        self.p, self.no_p_name = os.path.split(name)
        self.logfile = logfile
        if prmtop_name:
            self.prmtop_name = prmtop_name
        else:
            self.prmtop_name = '{}.prmtop'.format(self.no_p_name)
        if xvv:
            self.xvv_name = xvv
        else:
            self.xvv_name = 'water_{temp}.xvv'.format(temp=self.T)
        self.run_flags_list = None

    def setup_calculation(self, closure='hnc', write_g=False, write_h=False,
                          write_c=False,
                          write_u=False, write_asymp=False,
                          noasympcorr=False,
                          buffer_distance=25.0,
                          solvbox=False,
                          grdspc=(0.5, 0.5, 0.5),
                          tolerance=1e-5, polar_decomp=False,
                          verbose=0, maxstep=500, 
                          rism3d_path='rism3d.snglpnt'):
        """ Setup calculation rism3d.snglpnt. calculation.

        More details on each of the parameter can be found in AmberTools
        manual RISM section.

        Parameters
        ----------
        closure : string, default hnc
            Allowed closure values are kh, hnc, pseN. Here N is an
            integer.

        write_g : boolean, default False
            Specifies whether program will write radial distribution
            functions.

        write_h : boolean, default False
            Specifies whether program will write total correlation
            functions in k space.

        write_c : boolean, default False
            Specifies wheter program will write direct correlation
            functions.
            
        write_u : boolean, default False
            Specifies wheter program will write potential energy
            grid.
            
        write_asymp : boolean, default False
            Write asymptotics of total and direct correlation fuctions in 
            real space.         
            
        noasympcorr : boolean, default False
            Don't use long range corrections to compute thermodynamics.

        buffer_distance : float, default 25.0
            Minimum distance between the solute and the edge of
            the solvent box in A.

        solvbox : array-like (should contain 3 floats)
            Size of the box in x y and z directions. Overrides buffer_distance.

        grdsp: array-like (should contain 3 floats), default (0.5, 0.5, 0.5)
            Comma separated linear grid spacings for x, y and z dimensions.

        tolerance: float, default 1e-10
            Maximum residual values for solution convergence.

        polar_decomp: boolean, default False
            Decomposes solvation free energy into polar and non-polar
            components

        verbose: int, default 0
            Either 0, 1 or 2. Determines verbosity of caluclation.

        maxstep: int, default 1000
            Number of iterations in 3D-RISM calculation.
            
        rism3d_path : str, default rism3d.snglpnt
            Absolute path or exact name of rism3d.snglpnt program
        """
        grdspc = ','.join(map(str, grdspc))
        if solvbox:
            solvbox = ','.join(map(str, solvbox))
        self.run_flags_list = [rism3d_path,
             '--pdb', '{}.pdb'.format(self.no_p_name),
             '--prmtop', self.prmtop_name,
             '--xvv', self.xvv_name,
             '--grdspc', grdspc,]
        self.run_flags_list.extend(['--tolerance'] + tolerance)
        self.run_flags_list.extend(['--closure'] + closure)
        if solvbox:
            self.run_flags_list.extend(['--solvbox', solvbox])
        else:
            self.run_flags_list.extend(['--buffer', str(buffer_distance)])
        if write_g:
            self.run_flags_list.extend(['--guv',
                                         'g_{}'.format(self.no_p_name)])
        if write_h:
            self.run_flags_list.extend(['--huv',
                                         'h_{}'.format(self.no_p_name)])
        if write_c:
            self.run_flags_list.extend(['--cuv',
                                         'c_{}'.format(self.no_p_name)])
        if write_u:
            self.run_flags_list.extend(['--uuv',
                                         'u_{}'.format(self.no_p_name)])
        if write_asymp:
            self.run_flags_list.extend(['--asymp',
                                         'a_{}'.format(self.no_p_name)])
        if noasympcorr:
            self.run_flags_list.extend(['--noasympcorr'])
        if polar_decomp:
            self.run_flags_list.extend(['--polarDecomp'])
        if verbose:
            self.run_flags_list.extend(['--verbose'])
            self.run_flags_list.extend(['{}'.format(verbose)])
        if maxstep:
            self.run_flags_list.extend(['--maxstep'])
            self.run_flags_list.extend(['{}'.format(maxstep)])

    def run_calculation_and_log(self, timeout=30):
        """Run 3D-RISM single point calculation and log.

        Parameters
        ----------
        timeout : float, defult 30
            When HNC calculations get stuck they tend
            to run for a long amount of time. If calculation is
            running more than 30min most likely it is stuck.
            This option records 3D-RISM caclulation PID and
            kills it after supplied number of minutes. Doesn't work
            on windows. And in some other cases as well.
        """
        start_time = time.time()
        #print(self.run_flags_list)
        if IS_UNIX:  # use timeout
            run3drism = RunCmd(self.run_flags_list, timeout,
                               self.logfile, cwd=self.p)
            run3drism.run_and_timeout()
        else:    # windows
            rism_out = subprocess.check_output(self.run_flags_list, cwd=self.p)
            self.logfile.write(rism_out)
        self.logfile.flush()
        #write timestamp and close
        end_time = time.time()
        self.logfile.write(str(datetime.datetime.now()) + '\n')
        runtime = end_time - start_time
        self.logfile.write('3D-RISM runtime: {:.0f}'.format(runtime))
        self.logfile.flush()
        self.logfile.close()


def clean_up(name, T, level):
    """Delete junk.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    T: float
        A calculation temperature

    level : {0, 1, 2}
        0 - delete nothing;
        1 - delete ANTECHAMBER*, all water but .sh and .therm,
        .frcmod, .mol2, NEWPDB.PDB, PREP.INF, ATOMTYPE.INF, runleap.in
        sqm*, leap.log;
        2 - delete ALL but RESULTS_NAME and logfile - not recommended.
    """
    p, no_p_name = os.path.split(name)
    water_name = 'water_{}'.format(T)
    to_del1_glob = ['ANTECHAMBER*', 'sqm*', 'water*vv*']
    to_del1_files = [no_p_name + '.mol2', no_p_name + '.frcmod',
                    water_name + '.inp', water_name + '.out',
                    water_name + '.sav', 'ATOMTYPE.INF',
                    'leap.log', 'NEWPDB.PDB', 'PREP.INF', 'runleap.in',
                    MIN_SCRIPT_NAME, 'mdout', 'mdinfo']
    will_be_deleted_list = []
    if level == 1:
        for wildcard in to_del1_glob:
            will_be_deleted_list.extend(glob.glob(os.path.join(p, wildcard)))
        will_be_deleted_list.extend([os.path.join(p, f) for f in \
                                                            to_del1_files])
    if level == 2:
        all_files = os.listdir(p)
        all_files.remove(RESULTS_NAME)
        log_name = '{}.log'.format(no_p_name)
        all_files.remove(log_name)
        will_be_deleted_list.extend([os.path.join(p, f) for f in all_files])
    for f in will_be_deleted_list:
        try:
            os.unlink(f)
        except OSError, e:
            if e.errno == 2:
                pass
            else:
                raise e


def write_results(name, xvv_obj):
    """ Parses log file and writes free energies and corrections to
    results.txt.

    Parameters
    ----------
    name : string
        Full path to pdb file without extension

    xvv_obj : Xvv class instance
        Wrapper around xvv file used for calculation
    """
    p, _ = os.path.split(name)
    log_name = '{}.log'.format(name)
    exchem = None
    pmv = None
    with open(log_name, 'r') as f:
        for line in f:
            if line[0:11] == "rism_exchem":
                exchem = float(line.split()[1])
            if line[0:11] == "rism_exchGF":
                gf = float(line.split()[1])
            if line[0:11] == "rism_volume":
                pmv = float(line.split()[1])
    if not pmv:
        raise ValueError("Cannot find pmv value in log file. Most likely calculation didn't converge.")
    # compute PC
    pres, pres_plus = xvv_obj.compute_3drism_pressures() # [kcal/mol/A^3]
    PC = exchem - pres*pmv # pressure correction [kcal/mol] 
    PC_plus = exchem -pres_plus*pmv # pressure correction plus [kcal/mol]
    #Write and print results
    results = RESULTS.format(exchem=exchem, gf=gf, pmv=pmv, PC=PC, PC_plus=PC_plus,
                             pressure=pres, pressure_plus=pres_plus)
    with open(os.path.join(p, RESULTS_NAME), 'w') as f:
        f.write(results)
    print('Calculation has finished')
    print('RISM exchem={} kcal/mol'.format(exchem))
    print('PC+ dG*(solv)={} kcal/mol'.format(PC_plus))
    print('Detailed output can be found in {}.log'.format(name))
    return PC_plus


def main(argv):
    args = process_command_line(argv)
    for executable in REQUIRED_EXECUTABLES:
        if not distutils.spawn.find_executable(executable):
            raise NameError("{} is not found!".format(executable))
    print('Starting SFE calculation for {}'.format(args.file))
    name, dir_name = prepare_calc_directory(args.file, args.temperature, args.dir_name)
    logfile = prepare_logfile(name, argv)
    prmtop_name = prepare_prmtop(args, name, dir_name, logfile)
    #xvv is the path to xvv file relative to calc directory
    xvv = run_rism1d(name, logfile, args.temperature, args.smodel,
                     args.rism1d, args.closure[-1], xvv=args.xvv)
    xvv_obj = Xvv(os.path.join(dir_name, xvv))
    if args.minimize:
        minimize_solute(name, logfile, prmtop_name, args, xvv)
    rism_calc = RISM3D_Singlpnt(name, xvv_obj.temperature,
                                logfile, prmtop_name=prmtop_name,
                                xvv=xvv)
    rism_calc.setup_calculation(args.closure, 
                                write_g=args.write_g, 
                                write_h=args.write_h,
                                write_c=args.write_c,
                                write_u=args.write_u,
                                write_asymp=args.write_asymp,
                                noasympcorr=args.noasympcorr,
                                buffer_distance=args.buffer,
                                solvbox=args.solvbox,
                                grdspc=args.grdsp,
                                tolerance=args.tolerance,
                                polar_decomp=args.polar_decomp,
                                verbose=args.verbose3d,
                                maxstep=args.maxstep3d,
                                rism3d_path=args.rism3d_path)
    print('Running 3D-RISM calculation...')
    rism_calc.run_calculation_and_log(args.timeout)
    pc_plus = write_results(name, xvv_obj)
    clean_up(name, xvv_obj.temperature, args.clean_up)
    return pc_plus
    

if __name__ == '__main__':
    main(sys.argv[1:])
    


