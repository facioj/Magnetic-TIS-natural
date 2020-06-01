#!/usr/bin/python

import os
import numpy as np
from scipy.interpolate import interp1d

#some functions to  integrate density of states, find chemical potentials and DOS at a certain energy for a bunch of calculations leaving in folder of naming "folder_name_..."
def uncAddition(de,A=1.,B=1.,sigmaA=0.,sigmaB=0.):
    a = de/2
    b = de/2
    f =  a * A + b * B
    sigmaf = (a * sigmaA)**2 + (b * sigmaB)**2 
    return f, sigmaf


def identify_spin_blocks(data):
    i = 0
    line_0 = [] 

    for fila in data:
        if 'spin' in fila:
            line_0.append(i)
        i += 1

    assert len(line_0) == 2

    block1 = data[line_0[0]+1:line_0[1]-2]
    block2 = data[line_0[1]+1:len(data)-2]

    return block1,block2

class dos_tools():

    def find_index_in_gap(self):
        i=1
        for fila in self.data[1:]:
            vals = map(eval,fila.split())
            if(vals[1]**2 < 1e-16):
                break
            i += 1
        return i

    def find_top_valence_bands(self):
        i=self.index_in_gap
        while(i>0):
            fila = self.data[i]
            vals = map(eval,fila.split())
            if(i==self.index_in_gap and vals[1]>1e-8):
               assert 0
            elif(vals[1]>1e-8):
               return vals[0]
            i-=1

    def find_bottom_conduction_bands(self):
        i=0
        for fila in self.data[self.index_in_gap:]:
           vals = map(eval,fila.split())
           if(i==0 and vals[1]>1e-14):
               assert 0
           elif(vals[1]>1e-14):
               self.bottom_conduction_bands = vals[0]
               return vals[0]
           i+=1


    def trap_method(self,energies,f):
        """
        Args:

         - f: is a list of values to integrate with error estimation
         - energies: energy mesh

        Returns:
         - list containing partial integration and error 

        """
        assert len(energies)==len(f)
        factor = 1./self.volume
        de = (energies[1]-energies[0])*factor
        int_=0.
        int_StD=0.
        self.n_list = []
        en = []
        n = []
        for i in range(1,len(energies)):
            trapzPartAve,trapzPartStd = uncAddition(de=de,A=f[i-1][0],sigmaA=f[i-1][1], B=f[i][0],sigmaB=f[i][1])
            int_ += trapzPartAve 
            int_StD += trapzPartStd
            error_e = np.sqrt(int_StD)
            dos_e = 0.5  *(f[i-1][0] + f[i][0])
            en.append((energies[i]+energies[i-1])/2)
            n.append(int_)
            self.n_list.append([(energies[i]+energies[i-1])/2,int_,error_e,dos_e])
                              #mu, electron density, error in electron density, DOS

        fit_function = interp1d(en,n)
        return self.n_list,fit_function


    def dos_data_and_error(self,errorfile=None,bottom_error=None):
        #identify spin blocks
        block1,block2 = identify_spin_blocks(self.data)

        dos = []
        energies = []
        f_data = []
        if(self.holes==False):
            for i in range(len(block1)):
                fila_s1 = block1[i]
                fila_s2 = block2[i]
                vals_1 = map(eval,fila_s1.split())
                vals_2 = map(eval,fila_s2.split())
                assert vals_1[0] == vals_2[0]

                energy = vals_1[0]
                total_dos = abs(vals_1[1])+abs(vals_2[1])

                if(energy > self.bottom_conduction_bands):
                    energies.append(energy-self.bottom_conduction_bands)
                    dos.append(total_dos) 

            dos_estimate2 = []
            energies2 = []
            if(errorfile != None):
                file = open(errorfile,'r')
                data_error = file.readlines()
                file.close()
                bl_1,bl_2=identify_spin_blocks(data_error)
                for i in range(len(bl_1)):
                    fila_s1 = bl_1[i]
                    fila_s2 = bl_2[i]
                    vals_1 = map(eval,fila_s1.split())
                    vals_2 = map(eval,fila_s2.split())
                    assert vals_1[0] == vals_2[0]

                    energy = vals_1[0]
                    total_dos_estimate_2 = abs(vals_1[1])+abs(vals_2[1])
                
                    if(energy > bottom_error):
                        energies2.append(energy-bottom_error)
                        dos_estimate2.append(total_dos_estimate_2) 

                f  = interp1d(energies, dos)
                f2  = interp1d(energies2, dos_estimate2)
                new_energies = []
                for i in range(len(dos)):
                    x = energies[i]
                    if(x < max(energies2) and x > min(energies2)):
                        er = abs(f(x)-f2(x))
                        new_energies.append(x)
                        f_data.append([dos[i],er])
                return new_energies,f_data
            else:
                for i in range(len(dos)):
                    f_data.append([dos[i],0])

        return energies,f_data

 
    def int_dos(self,directory="."):
        factor = 1./self.volume 


        #compute step energy
        fila0=  self.data[5]
        vals0 = map(eval,fila0.split())
        fila1=  self.data[6]
        vals1 = map(eval,fila1.split())
        de = factor * (vals1[0] - vals0[0]) 

        bottom_conducting_bands = self.find_bottom_conduction_bands() 
        top_valence_bands = self.find_top_valence_bands() 

        #identify spin blocks
        block1,block2 = identify_spin_blocks(self.data)
        #print "#Two spin blocks identified of lenghts:", len(block1),len(block2)

        self.n_list = [] #list containing energy,integral and density
        #NEED TO CHANGE THIS INTO: PREPARE DATA FOR CALLING TETRA METHOD
        INT=0
        if(self.holes==False):
            for i in range(len(block1)):
                fila_s1 = block1[i]
                fila_s2 = block2[i]
                vals_1 = map(eval,fila_s1.split())
                vals_2 = map(eval,fila_s2.split())
                assert vals_1[0] == vals_2[0]

                energy = vals_1[0]

                if(energy > bottom_conducting_bands):
                  INT += de * (abs(vals_1[1])+abs(vals_2[1]))
                  self.n_list.append([energy,INT,abs(vals_1[1])+abs(vals_2[1])])
        else:
            i = self.index_in_gap
            while(i>0):
                fila_s1 = block1[i]
                fila_s2 = block2[i]
                vals_1 = map(eval,fila_s1.split())
                vals_2 = map(eval,fila_s2.split())
                assert vals_1[0] == vals_2[0]

                energy = vals_1[0]

                if(energy < top_valence_bands):
                  INT += de * (abs(vals_1[1])+abs(vals_2[1]))
                  self.n_list.append([energy,INT,abs(vals_1[1])+abs(vals_2[1])])
                i-=1

        if(directory != None):
            file = open("""%(directory)s/int_vs_mu.dat"""%locals(),'w')
            if(self.holes==True):
                file = open("""%(directory)s/int_vs_mu_holes.dat"""%locals(),'w')

            print >> file, "#mu, n [10^24cm^{-3}], dos"
            for pair in self.n_list:
                if(self.holes==True):
                  print >> file, pair[0]-top_valence_bands,pair[1],pair[2]
                else:
                  print >> file, pair[0]-bottom_conducting_bands,pair[1],pair[2]
            file.close()

        return self.n_list



    def __init__(self,file_name,volume,holes=False,directory=None):
        self.file_name = file_name
        self.volume = volume
        self.holes = holes

        file = open(self.file_name,'r')
        self.data = file.readlines()
        file.close()

        self.index_in_gap = self.find_index_in_gap()
        self.find_bottom_conduction_bands()
      #  self.int_dos(directory)

    def find_mu(self,searched_for):
        closest_to_searched_for = min(self.n_list, key = lambda t: abs(t[1]-searched_for))

        mu = closest_to_searched_for[0]
        dos = closest_to_searched_for[3]
        error = (closest_to_searched_for[1]-searched_for)**2
           
        return mu,dos,error


    def find_mu_int(self,searched_for):
        closest_to_searched_for = min(self.n_list, key = lambda t: abs(t[1]-searched_for))
        #build n vs e interpolation function

        #build dos vs e interpolation function
        mu = closest_to_searched_for[0]
        dos = closest_to_searched_for[3]
        error = (closest_to_searched_for[1]-searched_for)**2
           
        return mu,dos,error

    def shifted_dos(self):
        xs = []
        ys = []
        for fila in self.data:
            fs = fila.split()
            if(len(fs)==2):
                vals = map(eval,fila.split())
                xs.append(vals[0]-self.bottom_conduction_bands)
                ys.append(vals[1])
            if(len(fs)==0):
                break
        return xs,ys

    def make_dos_vs_density(self,file_output=None,densities=[1e18+i*1e18 for i in range(200)]):
    
        if(self.holes):
            file_name_ending = "_holes.dat"
        else:
            file_name_ending = ""

        if(file_output==None):
            file_output= """dos_vs_density_%(file_name_ending)s.dat"""%locals()
        file = open(file_output,'w')

        ans = []
        for dens in densities:
            mu,dos,error = self.find_mu(dens/1e24)
            ans.append([dens,dos,mu,error])
            print >> file, dens,dos,mu-self.bottom_conduction_bands,error
        return ans

def estimate_error(d1,d2,energy_interval=[-0.3,0.3],label='$D(\\varepsilon)$',label_error='Error',legend=None):

    import matplotlib.pyplot as plt

    x1,y1 = d1.shifted_dos()
    x2,y2 = d2.shifted_dos()

    Nps = 0
    for x in x1:
        if(x < energy_interval[1] and x > energy_interval[0]):
            Nps+=1

    x = np.linspace(energy_interval[0], energy_interval[1], Nps)
    f1  = interp1d(x1, y1)
    f2  = interp1d(x2, y2)

    fig, ax = plt.subplots()
    ax.set_xlabel('$\\varepsilon$ (eV)')
    plt.xlim(energy_interval[0],energy_interval[1])
    ax.plot(x, 0*f1(x), color='gray', linestyle='--',linewidth=0.5)
    ax.plot(x, f1(x), color='black',label=label)
    ax.fill_between(x, f1(x)-0.5*abs(f2(x)-f1(x)), f1(x)+0.5*abs(f2(x)-f1(x)) , facecolor='green', alpha=0.5, label=label_error)
    ax.legend(loc='upper left',title=legend)
    plt.show()


def compare_dos(d001,d001ap,d100,d100ap,energy_interval=[-0.3,0.3],label1='$D(\\varepsilon)$',label_error='Error',legend=None):
    """ Make plot comparing d001 and d100, including for each as errror estimation the difference with d001ap and d100ap
    """
    
    import matplotlib.pyplot as plt

    x001,y001 = d001.shifted_dos()
    x100,y100 = d100.shifted_dos()

    x001_ap,y001_ap = d001ap.shifted_dos()
    x100_ap,y100_ap = d100ap.shifted_dos()

    Nps = 0
    for x in x001:
        if(x < energy_interval[1] and x > energy_interval[0]):
            Nps+=1

    x = np.linspace(energy_interval[0], energy_interval[1], Nps)
    f001  = interp1d(x001, y001)
    f100  = interp1d(x100, y100)
    f001_ap  = interp1d(x001_ap, y001_ap)
    f100_ap  = interp1d(x100_ap, y100_ap)

    fig, ax = plt.subplots()
    ax.set_xlabel('$\\varepsilon$ (eV)')
    plt.xlim(energy_interval[0],energy_interval[1])
    ax.plot(x, 0*f001(x), color='gray', linestyle='--',linewidth=0.5)
    ax.plot(x, f001(x), color='black',label="001")
    ax.plot(x, f100(x), color='red',label="100")
    ax.fill_between(x, f001(x)-0.5*abs(f001(x)-f001_ap(x)), f001(x)+0.5*abs(f001(x)-f001_ap(x)) , facecolor='gray', alpha=0.5, label=label_error)
    ax.fill_between(x, f100(x)-0.5*abs(f100(x)-f100_ap(x)), f100(x)+0.5*abs(f100(x)-f100_ap(x)) , facecolor='orange', alpha=0.5, label=label_error)
    ax.legend(loc='upper left',title=legend)
    plt.show()

def anisotropy(d001,d001ap,d100,d100ap,energy_interval=[-0.3,0.3],label1='$D(\\varepsilon)$',label_error=None,legend=None,lab_mesh_1=None,lab_mesh_2=None):
    """ Make plot comparing d001 and d100, including for each as errror estimation the difference with d001ap and d100ap
    """
    
    import matplotlib.pyplot as plt

    #first, compare DOS
    x001,y001 = d001.shifted_dos()
    x100,y100 = d100.shifted_dos()

    x001_ap,y001_ap = d001ap.shifted_dos()
    x100_ap,y100_ap = d100ap.shifted_dos()

    Nps = 0
    x=[]
    for e in x001:
        if(e < energy_interval[1] and e > energy_interval[0]):
            x.append(e)
            Nps+=1

    f001  = interp1d(x001, y001)
    f100  = interp1d(x100, y100)
    f001_ap  = interp1d(x001_ap, y001_ap)
    f100_ap  = interp1d(x100_ap, y100_ap)

    fig, axs = plt.subplots(3,figsize=(4.9, 6))
    fig.subplots_adjust(hspace=0.4)

    axs[0].set_title(legend)
    axs[0].set_xlabel('$\\varepsilon$ (eV)')
    axs[0].set_ylabel('$D(\\varepsilon)$ (arb. unit.)')
    axs[0].set_xlim([energy_interval[0], energy_interval[1]])
    axs[0].plot(x, 0*f001(x), color='gray', linestyle='--',linewidth=0.5)
    axs[0].plot(x, f001(x), color='black',label="$m || \hat{z}$")
    axs[0].plot(x, f100(x), color='red',label="$m || \hat{x}$")
    axs[0].fill_between(x, f001(x)-0.5*abs(f001(x)-f001_ap(x)), f001(x)+0.5*abs(f001(x)-f001_ap(x)) , facecolor='gray', alpha=0.5, label=label_error)
    axs[0].fill_between(x, f100(x)-0.5*abs(f100(x)-f100_ap(x)), f100(x)+0.5*abs(f100(x)-f100_ap(x)) , facecolor='orange', alpha=0.5, label=label_error)
    axs[0].legend(loc='lower right')

    #second, compare DOS integrals
    en_001,f_data_001 = d001.dos_data_and_error(errorfile=d001ap.file_name,bottom_error=d001ap.bottom_conduction_bands)
    print en_001[0],en_001[-1]
    n_list_001,n_vs_e_001 = d001.trap_method(en_001,f_data_001)
    en_100,f_data_100 = d100.dos_data_and_error(errorfile=d100ap.file_name,bottom_error=d100ap.bottom_conduction_bands)
    print en_100[0],en_100[-1]
    n_list_100,n_vs_e_100 = d100.trap_method(en_100,f_data_100)

    xp=x[2:] 
    axs[1].set_xlabel('$\\varepsilon$ (eV)')
    axs[1].set_ylabel('$n$ (carrier/cm$^3$)')
    axs[1].set_xlim([0, energy_interval[1]])
    axs[1].set_yscale(value="log")
    axs[1].plot(xp, n_vs_e_001(xp)*1e24, color='black',label="001")
    axs[1].plot(xp, n_vs_e_100(xp)*1e24, color='red',label="100")
    file = open("""n_vs_energy_%(legend)s.dat"""%locals(),'w')

    for i in range(len(xp)):
        print >> file, xp[i],n_vs_e_001(xp[i])*1e24, n_vs_e_100(xp[i])*1e24
    file.close()

    #third, anisotropy ratio

    densities = []
    dens = 1e17
    while(dens < max(n_vs_e_100(xp)*1e24)):
        densities.append(dens)
        dens *= 1.1
    print densities

    e,f=d001ap.dos_data_and_error()
    en_ap,n_vs_e_001_ap = d001ap.trap_method(e,f)
    e,f=d100ap.dos_data_and_error()
    en_ap,n_vs_e_100_ap = d100ap.trap_method(e,f)

    ratio = []
    ratio2 = []
    dens = []
    y001 = n_vs_e_001(xp)
    y100 = n_vs_e_100(xp)
    y001_ap = n_vs_e_001_ap(xp)
    y100_ap = n_vs_e_100_ap(xp)
    doses = []
    doses_ap = []
    for d in densities:
        y001_d = np.array([abs(y*1e24-d) for y in y001])
        index_min=np.argmin(y001_d)
        mu = xp[index_min]
        dos_001 = f001(mu)

        y100_d = np.array([abs(y*1e24-d) for y in y100])
        index_min=np.argmin(y100_d)
        mu = xp[index_min]
        dos_100 = f100(mu)

        y001_d = np.array([abs(y*1e24-d) for y in y001_ap])
        index_min=np.argmin(y001_d)
        mu = xp[index_min]
        dos_001_ap = f001_ap(mu)

        y100_d = np.array([abs(y*1e24-d) for y in y100_ap])
        index_min=np.argmin(y100_d)
        mu = xp[index_min]
        dos_100_ap = f100_ap(mu)
        doses.append([dos_001,dos_100])
        doses_ap.append([dos_001_ap,dos_100_ap])
        ratio.append((dos_100/dos_001)**2)
        ratio2.append((dos_100_ap/dos_001_ap)**2)
        dens.append(d)

    file = open("""dos_vs_density_%(legend)s.dat"""%locals(),'w')
    for i in range(len(ratio)):
        print >> file, densities[i],ratio[i], ratio2[i], doses[i][0],doses[i][1],doses_ap[i][0],doses_ap[i][1] 
    file.close()
    #axs[2].set_xlabel('$\\n$ (electron/cm$^3$)')
    #axs[1].set_xlim([energy_interval[0], energy_interval[1]])
    axs[2].set_xscale(value="log")
    axs[2].set_xlabel('$n$ (carrier/cm$^3$)')
    axs[2].set_ylabel('$(D_{x}/D_{z})^2$')
    axs[2].plot(densities, [1 for i in range(len(densities))], color='gray', linestyle='--',linewidth=0.5)
    axs[2].plot(dens, ratio, color='black',label=lab_mesh_1,linestyle='dashed',fillstyle=None,markersize=2)
    axs[2].plot(dens, ratio2, color='blue',label=lab_mesh_2,linestyle='dashed',fillstyle=None,markersize=2)
    axs[2].set_ylim([0,6])
    axs[2].set_xlim([1e17,1e22])
    axs[2].yaxis.set_ticks([0,1,3,5])
    axs[2].legend(loc='upper right')
#    axs[2].set_xlim([min(densities),max(densities)])
    axs[2].fill_between(dens, ratio, ratio2 , facecolor='green', alpha=0.5, label=label_error)
  #  axs[1].fill_between(x, n_vs_e_100(x)-0.5*n_vs_e_100_error(x), n_vs_e_100(x)+0.5*n_vs_e_100_error(x) , facecolor='orange', alpha=0.5, label=label_error)
    plt.show()



