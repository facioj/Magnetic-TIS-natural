#!/usr/bin/python
import numpy as np

def get_func(k_center,enk,I,gamma,gamma_k):

    def lorentzian_k(k):
       return  gamma_k**2 / ( (k-k_center)**2 + gamma_k**2)

    def lorentzian(k,omega):
       return I * gamma**2 / ( (omega-enk)**2 + gamma**2) * lorentzian_k(k)

    return lorentzian

def compute_functions(file_name,E_max,E_min,gamma,gamma_k):
    file = open(file_name,'r')
    data = file.readlines()
    file.close()

    print "Computing functions"
    all_lor = []
    for fila in data[2:]:
      vals = map(eval,fila.split())
      if(len(vals)>1):
        if(vals[1]>E_min and vals[1] < E_max):
            ek = vals[1]
            weight = 0
            for j in range(2,len(vals)):
                weight += vals[j]
            A = get_func(vals[0],vals[1],weight,gamma,gamma_k)
            all_lor.append(A)

    return all_lor

def simulate_spectrum(file_name =  "bands_full.dat", k0=0,kf= 0.576,Nk=200,E_max=0.5,E_min=-0.5,Ne=200,gamma_k=0.002,gamma=0.004,lambda_0=20,orbital="pz",suffix=""):
    dk = (kf-k0)/Nk
    #initialize A
    momenta = np.linspace(k0, kf+dk, Nk)
    energies = np.linspace(E_min, E_max, Ne)
    A_final = []
    for i_k in range(len(momenta)):
        I_e = []
        for j_e in range(len(energies)):
            I_e.append(0)
        A_final.append(I_e)

    all_lor = compute_functions(file_name,E_max,E_min,gamma,gamma_k)
    print "Evaluating functions"
    for func in all_lor:
        s = np.vectorize(func)
        A = s(momenta[:,None],energies[None,:])
        A_final += A

    file_output = """A_gammak_%(gamma_k)s_gammae_%(gamma)s_Nk_%(Nk)s_Ne_%(Ne)s_lambda_%(lambda_0)s_%(orbital)s%(suffix)s"""%locals()
    file = open(file_output,'w')
    for i in range(len(momenta)):
      for j in range(len(energies)):
          print >> file,momenta[i],energies[j],A_final[i][j]
      print >> file,""
    file.close()

    return file_output

def sum_spectrum(file1,file2,file_output):
    file = open(file1,'r')
    data1 = file.readlines()
    file.close()
    file = open(file2,'r')
    data2 = file.readlines()
    file.close()

    assert len(data1) == len(data2)

    file = open(file_output,'w')

    for i in range(len(data1)):
        f1 = data1[i].split()
        f2 = data2[i].split()
        if(len(f1)>0):
            v1 = map(eval,f1)
            v2 = map(eval,f2)
            print >> file, v1[0],v1[1],v1[2]+v2[2]
        if(len(f1)==0):
            print >> file,"\n",
    file.close()


def plot_spectrum(x=None, y=None, z=None,file_name=None,title=None,save_to="out.png",xgamma=0,factor=1):
    import matplotlib.pyplot as plt

    fig= plt.figure(figsize=(3,4.5))

    xp = []
    yp = []
    zp = []

    if(x==None):
        if(file_name==None):
            print "Nothing to do"
            return 0
        else:
            file = open(file_name,'r')
            data = file.readlines()
            J=1
            K=1
            for fila in data:
                sp = fila.split()
                if(len(sp)==3):
                    vals = map(eval,sp)
                    if(J==1):
                        xp.append((vals[0]-xgamma)*factor)
                        J=0
                    if(K==1):
                        yp.append(vals[1])
                    zp.append(vals[2])
                if(len(sp)==0):
                    K=0
                    J=1
            x = np.array(xp)
            y = np.array(yp)
            z = np.array(zp)
            z = np.reshape(z,(len(x),len(y)))
            z = np.transpose(z)

    xv, yv = np.meshgrid(x,y)
    plt.xlabel('$k_y$ ($\AA^{-1}$)')
    plt.ylabel('$\\varepsilon$ (eV)')
    plt.xlim([-0.2,0.2])
    plt.ylim([-0.3,0.3])
    plt.title(title)
    plt.pcolormesh(xv, yv, z,cmap=plt.cm.YlOrBr)
    plt.savefig(save_to,dpi=300)
 



