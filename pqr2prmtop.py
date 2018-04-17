#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon, 16 Apr, 2018

@author: David M. Rogers
"""

from __future__ import print_function, division
from math import floor

import os
import argparse
import sys
import glob

try:
    #import parmed.tools as pact
    from parmed.topologyobjects import Atom,ResidueList,AtomType
    from parmed.structure import Structure
    from parmed.charmm import CharmmPsfFile
except ImportError as e:
    print('Make sure you are using AmberTools17!')
    raise e

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
    parser = argparse.ArgumentParser(description=""" Prepare Simple AMBER prmtop file.""")
    #Positional args
    parser.add_argument('pqr', metavar='molec.pqr',
                        help="""Input file, a PDB with all atoms followed by
                        3 space-delimited columns containing charge,
                        radius (Rmin), and epsilon params in units of
                        electrons, Angstrom, and kcal/mol (respectively).""")
    parser.add_argument('out', metavar='molec.prmtop',
                        help="""Name of the prmtop file to output.
                                No bonded interactions will be added.
                        """)
    return parser.parse_args(argv)

def read_pqr(filename):
    """ Read the atom/residue/chain names, coordinates and extra fields
        from a pqr file.

    Parameters
    ----------

    filename : PDB filename
        Path to the input file

    Returns
    -------

    name 
        [4-char]

    res
        [4-char]

    resn
        [int]

    chain
        [1-char]

    x
        [3-float]

    qre
        [float*] -- all remaining space-delimited floats

    """
    name  = [] # 4-char
    res   = [] # 4-char
    resn  = [] # int
    chain = [] # 1-char
    x     = [] # 3-float
    qre   = [] # float*
    with open(filename) as f:
      for line in f.xreadlines():
        if line[0:6] not in ["ATOM  ", "HETATM"]:
            continue

        name.append(line[12:16])
        res.append(line[17:21])
        chain.append(line[21])
        resn.append(int(line[22:26]))
        x.append( [ float(line[30:38]), float(line[38:46]),
                    float(line[46:54]) ] )
        qre.append(map(float, line[54:].split()))
    return name,res,resn,chain,x,qre

def mk_structure(args, name):
    name, res, resn, chain, x, qre = read_pqr(name)
    assert all(len(z) == 3 for z in qre), "All ATOM lines must end with 3 parameters: chg / e0, Rmin / Ang, and eps / kcal/mol."
    Q = sum(z[0] for z in qre)
    if abs(Q - floor(Q+0.5)) > 1e-6:
        print("Warning! Total charge on molecule is not an integer: %e"%Q)

    s = Structure()
    for n,(q,r,e) in zip(name, qre):
        s.add_atom(atom(n, q, e, r), 'MOL', 1, 'A')

    return s

def main(argv):
    args = process_command_line(argv)
    s = mk_structure(args, args.pqr)
    s.save(args.out, format='amber', overwrite=True)

type_serial = 0

# Create an atom type especially for this chg/LJ parameter combo
# Uses global parameter, type_serial!
def atom(name, chg, eps, rmin):
    global type_serial
    type_serial += 1
    at = AtomType(name, type_serial, 1.0, 1)
    at.set_lj_params(eps, rmin)
    a = Atom(name=name, type=name, charge=chg, mass=1.0,# Below 2 choices are
             solvent_radius=0.5*(rmin + 3.5537)/1.7, # not clearly defined.
             screen=0.8)
    a.atom_type = at
    return a

if __name__ == '__main__':
    main(sys.argv[1:])

