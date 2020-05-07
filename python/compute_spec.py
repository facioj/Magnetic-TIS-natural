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

def simulate_spectrum(file_name =  "bands_full.dat", k0=0,kf= 0.576,Nk=400,E_max=0.3,E_min=-0.3,Ne=400,gamma_k=0.002,gamma=0.004,lambda_0=20,orbital="pz"):
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

    file_output = """A_gammak_%(gamma_k)s_gammae_%(gamma)s_Nk_%(Nk)s_Ne_%(Ne)s_lambda_%(lambda_0)s_%(orbital)s"""%locals()
    file = open(file_output,'w')
    for i in range(len(momenta)):
      for j in range(len(energies)):
          print >> file,momenta[i],energies[j],A_final[i][j]
      print >> file,""
    file.close()

    return file_output

def plot_spectrum(x=None, y=None, z=None,file_name=None,title=None,save_to="out.png"):
    import matplotlib.pyplot as plt

    fig= plt.figure(figsize=(3,5))

    xp = []
    yp = []
    zp = []
    factor = 3.1415 / 4.3615 

    xgamma = 0.2 #np.sqrt(1./3+1./9)*factor

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
                        xp.append(vals[0]*factor-xgamma)
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
    plt.title(title)
    plt.pcolormesh(xv, yv, z,cmap=plt.cm.YlOrBr)
    plt.pcolormesh(-xv, yv, z,cmap=plt.cm.YlOrBr)
    plt.savefig(save_to,dpi=300)
 

#plot_spectrum(file_name="A_gammak_0.002_gammae_0.01_Nk_400_Ne_400_lambda_20_pz",title="$p_{3/2+1/2}, \lambda=10\AA$",save_to="U_3_32_12_lambda_10.png")
#assert 0

file_output = simulate_spectrum(file_name="+bwsum_orb_32_12_lambda_20",lambda_0=20,orbital="32_12",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{3/2+1/2}, \lambda=10\AA$",save_to="U_3_32_12_lambda_10.png")

file_output = simulate_spectrum(file_name="+bwsum_orb_32_-12_lambda_20",lambda_0=20,orbital="32_-12",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{3/2-1/2}, \lambda=10\AA$",save_to="U_3_32_-12_lambda_10.png")

file_output = simulate_spectrum(file_name="+bwsum_orb_12_12_lambda_20",lambda_0=20,orbital="12_12x",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{1/2+1/2}, \lambda=10\AA$",save_to="U_3_12_12_lambda_10.png")

file_output = simulate_spectrum(file_name="+bwsum_orb_12_-12_lambda_20",lambda_0=20,orbital="12_-12",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{1/2-1/2}, \lambda=10\AA$",save_to="U_3_12_-12_lambda_10.png")


file_output = simulate_spectrum(file_name="+bwsum_orb_32_32_lambda_20",lambda_0=20,orbital="32_32",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{3/2+3/2}, \lambda=10\AA$",save_to="U_3_32_32_lambda_10.png")

file_output = simulate_spectrum(file_name="+bwsum_orb_32_-32_lambda_20",lambda_0=20,orbital="32_-32",gamma=0.01)
plot_spectrum(file_name=file_output,title="$p_{3/2-3/2}, \lambda=10\AA$",save_to="U_3_32_-32_lambda_10.png")

