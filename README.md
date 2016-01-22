PC+
===

Python script to simplify computation of 3D-RISM pressure corrected (PC+/PC) hydration free energy. The 3D-RISM pressure is computed using the equation 20. from the following article: 

_Sergiievskyi, V.; Jeanmairet, G.; Levesque, M.; Borgis, D. Solvation Free-Energy Pressure Corrections in the Three Dimensional Reference Interaction Site Model. J. Chem. Phys. 2015, 143, 184116. http://dx.doi.org/10.1063/1.4935065_

If you find the script useful please cite the following works:

_Misin, M.; Fedorov, M.; Palmer, D. Hydration Free Energies of Ionic Species by Molecular Theory and Simulation, J. Phys. Chem. B. 2016. http://dx.doi.org/10.1021/acs.jpcb.5b10809_

and

_Misin, M.; Fedorov, M. V.; Palmer, D. S. Accurate Hydration Free Energies at a Wide Range of Temperatures from 3D-RISM. J. Chem. Phys. 2015, 142, 091105. http://dx.doi.org/10.1063/1.4914315_


Authors
-------
Maksim Mišin <mishin1991@gmail.com>

Usage examples
--------------

As an input script takes a pdb file with single solute. For example (methane.pdb):
```text
ATOM      1  C1  MOL     1       3.537   1.423   0.000  1.00  0.00
ATOM      2  H1  MOL     1       4.089   2.224   0.496  1.00  0.00
ATOM      3  H2  MOL     1       4.222   0.611  -0.254  1.00  0.00
ATOM      4  H3  MOL     1       2.759   1.049   0.669  1.00  0.00
ATOM      5  H4  MOL     1       3.077   1.810  -0.912  1.00  0.00
TER
END
```

1) 298.15 K methane hydration free energy calculation:

```
$ python rism3d_pc.py methane.pdb
Starting SFE calculation for methane.pdb at T=298.15 K
Running AM1-BCC calculation...
Running 1D-RISM calculation...
Running 3D-RISM calculation...
Calculation has finished
RISM exchem=8.57636028 kcal/mol
PC+ dG*(hyd)=1.07928069898 kcal/mol
Detailed output can be found in methane_298.15/methane.log
```

2) 350 K calculation with tip3p water
	
```
$ python rism3d_pc.py methane.pdb -t 350 --smodel TP3
Starting SFE calculation for methane.pdb at T=350.0 K
Running AM1-BCC calculation...
Running 1D-RISM calculation...
Running 3D-RISM calculation...
Calculation has finished
RISM exchem=9.01341914 kcal/mol
PC+ dG*(hyd)=1.35551076507 kcal/mol
Detailed output can be found in methane_350.0/methane.log
```
	
3) Using existing topology and susceptibility (xvv) files

```
$ python rism3d_pc.py methane.pdb -p methane.prmtop -x water_nacl.xvv
Starting SFE calculation for methane.pdb at T=298.15 K
Reading user provided prmtop file
Using provided xvv file...
Running 3D-RISM calculation...
Calculation has finished
RISM exchem=8.67168499 kcal/mol
PC+ dG*(hyd)=1.1103088454 kcal/mol
Detailed output can be found in methane_298.15/methane.log
```

Prerequisites
-------------

The script requires:

* Python 2.7 or later: http://www.python.org/
* AmberTools13 or later: http://ambermd.org/


Get some help
-------------

    $ python rism3d_pc.py -h
    usage: rism3d_pc.py [-h] [-p PRMTOP] [--scale_chg SCALE_CHG] [-c MOLCHARGE]
                        [--multiplicity MULTIPLICITY] [--minimize MINIMIZE]
                        [-x XVV] [--smodel SMODEL] [--rism1d RISM1D]
                        [--closure [CLOSURE [CLOSURE ...]]] [-t TEMPERATURE]
                        [--clean_up CLEAN_UP] [--dir_name DIR_NAME]
                        [--timeout TIMEOUT]
                        [--tolerance [TOLERANCE [TOLERANCE ...]]] [--write_g]
                        [--write_c] [--write_u] [--write_asymp] [--noasympcorr]
                        [--buffer BUFFER] [--solvbox SOLVBOX SOLVBOX SOLVBOX]
                        [--grdsp GRDSP GRDSP GRDSP] [--polar_decomp]
                        [--verbose3d VERBOSE3D] [--maxstep3d MAXSTEP3D]
                        molec.pdb

    Run 3D-RISM single point calculation. This script is a wrapper around Amber
    rism3d.snglpnt program, designed to simplify calculations. It also computes
    pressure corrections to RISM solvation free energy.

    optional arguments:
      -h, --help            show this help message and exit

    Solute options:
      Options related to a solute molecule. Calculation requires only pdb and
      prmtop files. If prmtop file is not present, script will try to create
      this file using antechamber and tleap. By default AM1-BCC charges and GAFF
      will be used.

      molec.pdb             Input solute file. Must be in pdb format acceptable by
                            Antechamber. Must have a .pdb extension.
      -p PRMTOP, --prmtop PRMTOP
                            Path to parameters and topology (prmtop) file of
                            solute.
      --scale_chg SCALE_CHG
                            Scale all solute by this value prmtop file [1.0].
      -c MOLCHARGE, --molcharge MOLCHARGE
                            Charge of the solute [0].
      --multiplicity MULTIPLICITY
                            Multiplicity of the solute [1].
      --minimize MINIMIZE   Minimize solute before performing 3D-RISM calculation
                            using either gradient descent (min) or RISM
                            minimization using sander (rism, not recommended). If
                            no keywords are provided minimization is not
                            performed.

    1D-RISM options:
      3D-RISM calculation requires xvv file (site-site susceptibilities) to run.
      If such file is not provided script can run a 1D-RISM to try to genrate
      it, but will asume that the solvent is pure water. Some of the 1D-RISM
      related options are in 3D-RISM group.

      -x XVV, --xvv XVV     Path to an existing xvv file. This will skip 1D-RISM
                            calculation and all related parameters will be
                            ignored. Solvent density as well as calculation
                            temeprature will be read from this file.
      --smodel SMODEL       Solvent model for 1D-RISM calculation available in
                            "$AMBERHOME/dat/rism1d/model/{smodel}.mdl" [SPC].
      --rism1d RISM1D       Type of 1D-RISM theory. Only DRISM has been
                            extensively tested [DRISM].

    3D-RISM options:
      Options related to the main calculation.

      --closure [CLOSURE [CLOSURE ...]]
                            Brdige closure for 3D-RISM and 1D-RISM calculation if
                            it is necessary. Either HNC, KH, PSEn (n in PSEn
                            should be an integer) or a list of them for sequential
                            convergence. For 1D-RISM calculation only last closure
                            will be used [PSE3].
      -t TEMPERATURE, --temperature TEMPERATURE
                            Temperature in K at which calculation will be
                            preformed. If xvv file was provided, this option will
                            be used only for naming directory [298.15].
      --clean_up CLEAN_UP   How should auxiliary files be treated: 0 - delete
                            nothing; 1 - delete some [default]; 2 - delete all but
                            input, results, and log.
      --dir_name DIR_NAME   Custom name for produced calculation directory. The
                            default one is: {mol_name}_{temperature}.
      --timeout TIMEOUT     Minutes after which 3D-RISM calculation will be
                            killed. Use 0 for no timeout. [0]. Only works on Unix-
                            like system.
      --tolerance [TOLERANCE [TOLERANCE ...]]
                            Maximum residual values for 3D-RISM solution
                            convergence. If many closures a list of closures can
                            be supplied [1E-5].
      --write_g             Write radial distribution functions produced in 3D-
                            RISM calculation.
      --write_c             Write direct correlation function produced in 3D-RISM
                            calculation.
      --write_u             Write solute solvent potential energy grid.
      --write_asymp         Write asymptotics of total and direct correlation
                            fuctions in real space.
      --noasympcorr         Thermodynamics of 3D-RISM is calculated without long-
                            range asymptotics.
      --buffer BUFFER       Minimum distance between the solute and the edge of
                            the solvent box in A for 3D-RISM calculation [25].
      --solvbox SOLVBOX SOLVBOX SOLVBOX
                            Size of the x, y, and z dimensions of the box in
                            Angstroms. Specifying this parameter overrides buffer.
      --grdsp GRDSP GRDSP GRDSP
                            Linear grid spacings for x, y and z dimensions. Should
                            be separated with spaces. Units: A [0.5 0.5 0.5].
      --polar_decomp        Decomposes solvation free energy into polar and non-
                            polar components.
      --verbose3d VERBOSE3D
                            Verbosity of 3D-RISM calculation. 0 - print nothing; 1
                            - print iterations; 2 - print all [2].
      --maxstep3d MAXSTEP3D
                            Maximum number of iterations in 3D-RISM calculation
                            [500].


Notes
-----
The script has been tested only on Ubuntu, but it should work on most Linux distributions and on Mac OS X. To make it Windows-friendly one would probably need to change the names of executable programs and add them to the PATH.



