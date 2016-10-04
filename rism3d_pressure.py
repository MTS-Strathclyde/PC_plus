#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 14:57:09 2016

Compute pressure using varous methods.
@author: max
"""

import numpy as np
import sys
import pubfft
from scipy import fftpack
import matplotlib.pyplot as plt
import argparse

K_B = 1.9872041E-3 # boltzmann const in kcal/mol/K
N_A = 6.022141e23 # avogadro's constant


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
    parser = argparse.ArgumentParser(description="""Calculate 3D-RISM 
            pressure and a number of related quantities.""")
    #Positional args
    parser.add_argument('xvv', metavar='solvent.xvv',
                        help="""Site-site susceptibility file in Amber
                        1D-RISM format.""")
    #Optional args
    parser.add_argument('-t', '--therm',
                        help=""" File with 1D-RISM thermodynamic output""")    
    parser.add_argument('-c', '--cvv',
                        help=""" File with 1D-RISM cvv functions""")    
    parser.add_argument('-k', '--k_number',
                        help=""" Compute pressure at k = dk*k_number [0]""",
                        default=0, type=int)
    parser.add_argument('--plot_xvvk',
                        help=""" Create a plot of X_{ij} (k) for each i and j""",
                        action='store_true')
    parser.add_argument('--save_xvvk',
                        help=""" Name of the file in which xvv_{ij} (k) will be
                        stored""")
    parser.add_argument('--plot_xvvr',
                        help=""" Create a plot of X_{ij} (r) for each i and j""",
                        action='store_true')
    parser.add_argument('--save_xvvr',
                        help=""" Name of the file in which xvv_{ij} (r) will be
                        stored""")
    parser.add_argument('--plot_zr',
                        help=""" Create a plot of z_{ij} (r) for each i and j""",
                        action='store_true')
    parser.add_argument('--plot_zk',
                        help=""" Create a plot of z_{ij} (k) for each i and j""",
                        action='store_true')
    parser.add_argument('--save_zr',
                        help=""" Name of the file in which z_{ij} (r) will be
                        stored""")
    parser.add_argument('--save_zk',
                        help=""" Name of the file in which z_{ij} (k) will be
                        stored""")
    return parser.parse_args(argv)



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
        # Unique sites per species
        self.unique_sites_per_species = None
        # Actual number of sites per species
        self.total_sites_per_species = None
        # Density of each species
        self.species_densities = None
        # Sites from the same species have densities equal to densities of species
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

    def get_k(self):
        """ Return array of k values"""
        dk = np.pi/(self.ngrid*self.dr)
        return np.array([i*dk for i in range(self.ngrid)])


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


    def _get_compressibility_from_therm(self, therm_p):
        """ Return solvent compressibility.
    
        Parameters
        ----------
        therm_p : string
            path to .therm file
    
        Returns
        -------
        compres : float
            Units: 1/MPa
        """
        with open(therm_p, 'r') as f:
            therm_lines = f.readlines()
        compres = float(therm_lines[2].split()[-1])
        units = therm_lines[2].split()[-2]
        if units == '[10e-4/MPa]':
            # Old version of Ambertools
            return compres*10e-4
        if units == '[1/kPa]':
            # # !! This is ambertools 14, where compressiblity is bugged.
            # Units are shown to be [1/kPa], while in reality compres thea are in [1/MPa]
            # http://archive.ambermd.org/201503/0651.html
            return compres
        if units == '[1/MPa]':
            # This is ambertools 15
            # All is good
            return compres
        else:
            raise ValueError('Unknown compressiblity format, check *.therm file')


    def compute_3drism_pressures(self, k=0):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 20 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        k : int
            Which k value to use to compute pressure. 1 is recommended.
            
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
        #print self.total_sites_per_species
        ideal_pressure  = sum(self.species_densities)*K_B*self.temperature
        #print 'under_pressure',pressure - 2*ideal_pressure
        return pressure - ideal_pressure, pressure
        
    def compute_compres_pressure(self, therm_p):
        """ Compute 3drism pressure using loaded xvv file.
        Uses equation 21 from the article by Sergiievskyi et al. 
        (http://dx.doi.org/10.1063/1.4935065). 

        Parameters
        ----------
        therm_p : string
            path to .therm file
            
        Return
        ------
        pressures : tuple of floats
            Tuple containeing two pressures.
            First element is 3D-RISM pressure minus ideal gas pressure (used in PC+),
            second is 3D-RISM pressure (used in PC).
            Both have units of kcal/mol/A^3.
        """
        compres = self._get_compressibility_from_therm(therm_p) # 1/MPa
        # 4.184e24 A^3 = kcal/MPa
#        C_k0_rho2 = sum(self.species_densities) - \
#                    (4.184e24*self.temperature*K_B/N_A*compres)**(-1) #density**2*c(k=0)
        C_k0_rho = 1 - \
                    (sum(self.species_densities)*1.0e24*self.temperature*K_B/N_A*compres*4184)**(-1) #density**2*c(k=0)
        C_k0_rho2 = C_k0_rho*sum(self.species_densities)
        # pressure using eq. 21
        # works for mixtures as well!
        pressure = K_B*self.temperature*sum(np.array(self.species_densities)\
                   *(np.array(self.total_sites_per_species)+1))/2.\
                 - K_B*self.temperature/2.*C_k0_rho2
        ideal_pressure = sum(self.species_densities)*K_B*self.temperature
        return pressure - ideal_pressure, pressure
        
    def compute_cvv_pressure(self, cvv_fname):
        cvv = np.loadtxt(cvv_fname).astype(float)
        dr = cvv[1,0] - cvv[0,0]
        C_k0_rho2 = 0
        for i in range(self.nsites):
            for j in range(i, self.nsites):
                interaction = i + j + 1
                integral = np.sum(cvv[:, interaction]*cvv[:, 0]**2*4*np.pi)*dr
                if i != j:
                    integral = integral*2
                rho2 = self.multiplicities[i]*self.multiplicities[j]*\
                       self.normalized_densities[i]*self.normalized_densities[j]
                C_k0_rho2 += rho2*integral
        pressure = K_B*self.temperature*sum(np.array(self.species_densities)\
                   *(np.array(self.total_sites_per_species)+1))/2.\
                 - K_B*self.temperature/2.*C_k0_rho2
        ideal_pressure = sum(self.species_densities)*K_B*self.temperature
        return pressure - ideal_pressure, pressure
        
    def compute_zk(self):
        """ Return z(k) matrix """
        zk = []
        dk = np.pi/(self.ngrid*self.dr)
        k = np.array([i*dk for i in range(self.ngrid)])
        mult_vec = np.array(self.multiplicities)
        density_vec = np.array(self.normalized_densities)
#        print density_vec
#        for k in range(self.ngrid):
#            if k%100 == 0:
#                print k
#            xvv_k = self.xvv_data[k,:,:]
#            # Z_k from sergievskyi's article
#            zk_i = mult_vec/density_vec*(np.identity(self.nsites) - np.linalg.inv(xvv_k))
#            #print zk_i.shape
#            zk.append(zk_i)
        zk = [mult_vec/density_vec*(np.identity(self.nsites) - \
                          np.linalg.inv(self.xvv_data[i,:,:])) \
                                                    for i in range(self.ngrid)]
        zk = np.array(zk)
        return k, np.array(zk)
        
    def compute_xvvr(self):
        """ Return xvv(r) matrix """
        r = np.array([i*self.dr for i in range(self.ngrid)])
        k = self.get_k()
        xvvr = [["" for i in range(self.nsites)] for j in range(self.nsites)]
        for i in range(self.nsites):
            for j in range(self.nsites):
                xvvk_ij = self.xvv_data[:,i,j]
                xvvr_ij = pubfft.sinfti(xvvk_ij*k, self.dr, -1)/r
#                n_pots_for_interp = 6
#                r_for_interp = r[1:n_pots_for_interp+1]
#                xvvr_for_interp = xvvr_ij[:n_pots_for_interp]
#                poly_coefs = np.polyfit(r_for_interp, xvvr_for_interp, 3)
#                poly_f = np.poly1d(poly_coefs)
#                xvvr[i][j] = [poly_f(0)]
                xvvr[i][j] = xvvr_ij
        return r, np.swapaxes(xvvr, 0, 2)

    def compute_zr(self):
        """ Return z(r) matrix """
        r = np.array([i*self.dr for i in range(self.ngrid)])
        k, zk = self.compute_zk()
        print 'computed zk',zk.shape
        zr = [["" for i in range(self.nsites)] for j in range(self.nsites)]
        for i in range(self.nsites):
            for j in range(self.nsites):
                zk_ij = zk[1:,i,j]
                zr_ij = pubfft.sinfti(zk_ij*k[1:], self.dr, -1)/r[1:]
                #zr_ij = np.abs(fftpack.fft(zk_ij))
                n_pots_for_interp = 6
                r_for_interp = r[1:n_pots_for_interp+1]
                zr_for_interp = zr_ij[:n_pots_for_interp]
                poly_coefs = np.polyfit(r_for_interp, zr_for_interp, 3)
                poly_f = np.poly1d(poly_coefs)
                zr[i][j] = [poly_f(0)]
                zr[i][j].extend(zr_ij)
        return r, np.swapaxes(zr, 0, 2)
        
        
def plot_and_save(x, func, xvv_inst, plot=False, fname=False):
    mat = x.T
    for m in range(xvv_inst.nsites):
        for n in range(m+1):
            if fname:
                mat = np.c_[mat, func[:,m,n].T]
            if plot:
                plt.plot(x, func[:,m,n], 
                         label='{}-{}'.format(xvv_inst.atom_names[m], 
                                              xvv_inst.atom_names[n]))
    if fname:
        np.savetxt(fname, mat)
    if plot:
        plt.legend()
        plt.savefig('graph.png', dpi=300)
        plt.show()    
        

def main(argv):
    args = process_command_line(argv)
    xvv_inst = Xvv(args.xvv)
    pres_plus, pres = xvv_inst.compute_3drism_pressures(args.k_number)
    print 'Pressure+: {:.6f} kcal/mol/A^3'.format(pres_plus)
    print 'Pressure:  {:.6f} kcal/mol/A^3'.format(pres)
    kcal_per_a_cubed_to_bar = 4.184e28/N_A
    print 'Pressure+: {:.6f} bar'.format(pres_plus*kcal_per_a_cubed_to_bar)
    print 'Pressure:  {:.6f} bar'.format(pres*kcal_per_a_cubed_to_bar)
    if args.therm:
        cpres_plus, cpres = xvv_inst.compute_compres_pressure(args.therm)
        print 'compres Pressure+: {:.6f} kcal/mol/A^3'.format(cpres_plus)
        print 'compres Pressure:  {:.6f} kcal/mol/A^3'.format(cpres)
    if args.cvv:
        cvvpres_plus, cvvpres = xvv_inst.compute_cvv_pressure(args.cvv)
        print 'cvv Pressure+: {:.6f} kcal/mol/A^3'.format(cvvpres_plus)
        print 'cvv Pressure:  {:.6f} kcal/mol/A^3'.format(cvvpres)
    analyze_xvvr = args.plot_xvvr or args.save_xvvr
    if analyze_xvvr:
        r, xvvr = xvv_inst.compute_xvvr()
        plot_and_save(r, xvvr, xvv_inst, args.plot_xvvr, args.save_xvvr)  
    analyze_xvvk = args.plot_xvvk or args.save_xvvk
    if analyze_xvvk:
        k = xvv_inst.get_k()
        plot_and_save(k, xvv_inst.xvv_data, xvv_inst, args.plot_xvvk, args.save_xvvk)
    analyze_zr = args.plot_zr or args.save_zr
    if analyze_zr:
        r, zr = xvv_inst.compute_zr()
        plot_and_save(r, zr, xvv_inst, args.plot_zr, args.save_zr)
    analyze_zk = args.plot_zk or args.save_zk
    if analyze_zk:
        k, zk = xvv_inst.compute_zk()
        plot_and_save(k, zk, xvv_inst, args.plot_zk, args.save_zk)
    
    #
    #pressures = []
    #for i in range(40):
    #    pressures.append(xvv_inst.compute_3drism_pressures(i))
    #plt.plot(pressures, '-o')
    #plt.show()
        
        
if __name__ == '__main__':
    main(sys.argv[1:])
    
    
    
