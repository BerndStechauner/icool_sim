import matplotlib
import numpy as np
import scipy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
from scipy.interpolate import spline
from scipy.stats import norm
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.optimize import curve_fit

fig, ax = plt.subplots()

#define_Lattice_fuction_for solenoid___________________________________
#______Solenoid Matrix_________________________________________________
#Bz in [T]
#Brho = Regidity: longitudinal momentum in MeV devided 299.792 [MeV]
def SOL(L,Brho,Bz):
    K = Bz/(2*Brho)
    C = np.cos(K*L)
    S = np.sin(K*L)
    return [[C, S/K], [-K*S, C]]
#_______________________________________________________________________
#define_Twiss_transforamtion_matrix_____________________________________
def TWISS(M):
    C = M[0][0]
    S = M[0][1]
    Cprime = M[1][0]
    Sprime = M[1][1]
    return [[C*C, -2*C*S, S*S], [-C*Cprime, S*Cprime + Sprime*C, -S*Sprime], [Cprime*Cprime, -2*Sprime*Cprime, Sprime*Sprime]]
#_______________________________________________________________________
#______Solenoid strength________________________________________________
def k(B,p):
    return 0.299792458 * B * 0.5 / p
#_______________________________________________________________________

B00 = 48                #high field in [T]
B_low0 = 4              #low field in [T] (offset B)
L0 = 2.2                #solenoid length in [m]
sym0 = 2.195
pz = 0.1765             #beam momentum in [GeV]
parnum = 10000          #particle number for for001.dat
#Absorber material
AbsMaterial = 'LH'
#Initial emittance in [m]
norm_emi = 0.0003
step0 = 0.01

#Optimization range of solenoid length
R = np.arange(0.01,0.1, 0.0001)

def f(R):
    B0 = B00
    B_low = B_low0       #low field in [T]
    L = L0          #length of the solenoid in [m]
    step = step0    #step size in [m]
    Bz = []
    z = []
    sym = sym0
    init = -sym      #initial
    end = sym        #end
    for i in range(0,int((end-init)/step)):
        f1 = (init+0.5*L)/(np.sqrt(R*R+(init+0.5*L)**2))
        f2 = (init-0.5*L)/(np.sqrt(R*R+(init-0.5*L)**2))
        b = B0*0.5*(f1-f2)
        if init+end <=1:
            z.append(init+end)
            Bz.append(b+B_low)
        init = init + step
    B0 = Bz[0]           #Magnetic field at postion z=0 in [T]
    B1 = Bz[len(Bz)-1]
    Brho = pz/0.299792458
    vec0 = [1/k(B0,pz), 0, k(B0,pz)]
    vec1 = [1/k(B1,pz), 0, k(B1,pz)]
    #Calculate_lattice_matrix_______________________________________________
    Unit = [[1,0], [0,1]]
    for i in range(0,len(Bz)):
        M = SOL(step, Brho, Bz[i])
        Unit = np.dot(M,Unit)
    #Calculate_Twiss_matrix_________________________________________________
    T = TWISS(Unit)
    return vec1[1] + vec1[0] + vec1[2] - vec0[0]*T[0][0] - vec0[0]*T[1][0] - vec0[0]*T[2][0] - vec0[1]*T[1][1] - vec0[1]*T[0][1] - vec0[1]*T[2][1] - vec0[2]*T[1][2] - vec0[2]*T[0][2] - vec0[2]*T[2][2]
#____________________________________________________________________________
#___Bernd_searching_method___________________________________________________
#____________________________________________________________________________
ff = []
for i in range(0,len(R)):
    ff.append(f(R[i]))
    
fff = []
for i in range(0,len(ff)):
    fff.append(abs(ff[i]))
min = np.amin(fff)
print 'minimal optimization value:', min
for i in range(0,len(ff)):
    if fff[i] == min:
        RR = R[i]
        print 'Solenoid radius=', R[i]
        
#______________________PART_I__________________________________________
#___Plot_beta_function_________________________________________________
#______________________________________________________________________

B0 = B00
R = RR
B_low = B_low0       #low field in [T]
L = L0        #length of the solenoid in [m]
step = step0    #step size in [m]
Bz = []
z = []
sym = sym0
init = -sym      #initial
end = sym        #end
for i in range(0,int((end-init)/step)):
    f1 = (init+0.5*L)/(np.sqrt(R*R+(init+0.5*L)**2))
    f2 = (init-0.5*L)/(np.sqrt(R*R+(init-0.5*L)**2))
    b = B0*0.5*(f1-f2)
    if init+end <=1:
        z.append(init+end)
        Bz.append(b+B_low)
    init = init + step
B0 = Bz[0]           #Magnetic field at postion z=0 in [T]
B1 = Bz[len(Bz)-1]
Brho = pz/0.299792458

b = []
u = []
zzz = []

b.append(1/k(Bz[0],pz))
u.append(0)
zzz.append(0)
L = step
dz = 1e-5

const = 0.1
for i in range(0,int(const/step)):
    z.append(z[len(z)-1]+dz)
    Bz.append(Bz[len(Bz)-1])

z0 = 0
for j in range(0, len(Bz)):
    anfang =  int(L*j/dz)
    ende = int(L*(j+1)/dz)
    k = (0.5*0.29979258)*Bz[j]/pz
    for i in range(anfang,ende):
#        BB.append(B[j])
        bb = b[i] + u[i]*dz
        uu = u[i] + 2*dz/b[i] + u[i]**2*dz/(2*b[i]) - 2*b[i]*k**2*dz
        z0 = z0 + dz
        b.append(bb)
        u.append(uu)
        zzz.append(z0)


plt.subplot(212)
plt.plot(zzz,b, label=r'$\beta$ [m]', lw=2)
plt.plot(zzz,u, label=r'$\beta^{\prime}$', lw=2)
plt.xlabel(r"$z$ [m]", fontsize=15)
plt.legend(loc = 'best', numpoints=1)
plt.axhline(y=0, color='k')

#______________________________________________________________________
#built_a_complete_realistic_solenoidal_field___________________________
#______________________________________________________________________
B0 = B00
R = RR
B_low = B_low0       #low field in [T]
L = L0               #length of the solenoid in [m]
step = step0         #step size in [m]
Bz = []
z = []
sym = sym0
init = -sym      #initial
end = sym        #end
for i in range(0,int((end-init)/step)):
    f1 = (init+0.5*L)/(np.sqrt(R*R+(init+0.5*L)**2))
    f2 = (init-0.5*L)/(np.sqrt(R*R+(init-0.5*L)**2))
    b = B0*0.5*(f1-f2)
    init = init + step
    z.append(init+end)
    Bz.append(b+B_low)
#    if init+end >1:
#        z.append(init+end-1)
#        Bz.append(b+B_low)
#______________________________________________________________________
#______________________________________________________________________

#______________________________________________________________________
#___tanh_solenoidal_field_approximation________________________________
#______________________________________________________________________
def Btanh(z, Bc, Lc, elen, alen, Boff):
#   Bc = field strength [T]
#   Lc = Length of central region [m]
# elne = Length of entrence end rgeion [m]
# alen = End attenuation length [m]
# Boff = Magnetic offset
    a = (Bc*0.5 + Bc*0.5*np.tanh((z - elen)/alen))*np.heaviside((elen+Lc*0.5)-z,1)
    b = (Bc*0.5 + Bc*0.5*np.tanh((elen - z + Lc)/alen))*np.heaviside(z-(elen+Lc*0.5),1)
    return a + b + Boff

popt, pcov = curve_fit(Btanh, z, Bz, p0=(B0, L, sym, L/R, B_low))
print 'High field strength:', popt[0] ,'T'
print 'Central region length:', popt[1], 'm'
print 'End region length:', popt[2], 'm'
print 'End attenuation length:', popt[3], 'm'
print 'Constant offset field strength:', popt[4], 'T'

plt.subplot(211)
plt.plot(z,Bz)
plt.plot(z, Btanh(z,popt[0], popt[1], popt[2], popt[3], popt[4]), 'g--')

#______________________PART_II_________________________________________
#___Create_a_ICOOL_lattice_____________________________________________
#______________________________________________________________________

#for001 file creator
#Material in finge field
FMaterial = 'VAC'
c = 0.29979258      #speed of light 1e-9
#Mean values in [GeV]
x =  0
y = 0
z = 0
px = 0.
py = 0.
#pz = 0.1765
#sigma values in [GeV]
beta0 = pz*2/(0.29979258*Bz[0])
m_mu = 0.1056583755
eisbaer = np.sqrt(pz*norm_emi*m_mu/(beta0))
braunbaer = np.sqrt(norm_emi*beta0*m_mu/pz)
sigmax = braunbaer
sigmay = braunbaer
sigmaz = 0.
sigmapx = eisbaer
sigmapy = eisbaer
sigmapz = 0.
#particle number
f = open('for001.dat','w')

#Header____________________________________________________________

#Control Variables_________________________________________________
A = 0.45
f.write('&SUB apertr  %f' %(A))
f.write('\n')
f.write('&SUB radius  %f' %(R))
f.write('\n')
f.write('&cont npart=%d bgen=.true. dectrk=.true. f9dp=8' %(parnum))
f.write('\n')
f.write('varstep=.true. nprnt=-1 prlevel=1 ntuple=.false.')
f.write('\n')
f.write('spin=.false. spintrk=0 spinmatter=0 bzfldprd=20')
f.write('\n')
f.write('fsav=.true. fsavset=.true. izfile=581 ! rfdiag=89')
f.write('\n')
f.write('phasemodel=3 output1=.true. /')
f.write('\n')
f.write('&bmt nbeamtyp=1 /')
f.write('\n')
#____________________________________________________________________

#Beam parameters_____________________________________________________
f.write('1 2 1. 1')         #Gaussian distribution
f.write('\n')
#________x__y__z__px_py_pz____mean values____________________________
f.write('%f %f %f %f %f %f' %(x, y, z, px, py, pz))
f.write('\n')
#________x__y__z__px_py_pz____sigma values___________________________
f.write('%.9f %.9f %.9f %.9f %.9f %.9f' %(sigmax, sigmay, sigmaz, sigmapx, sigmapy, sigmapz))
f.write('\n')
f.write('0')
f.write('\n')
f.write('0 0 0 0')
f.write('\n')

#____________________________________________________________________
f.write('&ints ldecay=.true. ldedx=.true. lstrag=.true. lscatter=.true. ldray=.false. lspace=F')
f.write('\n')
f.write('delev=2 straglev=4 scatlev=1 declev=1 spacelev=2 parbunsc=4.81e12 /')
f.write('\n')
#____________________________________________________________________

#____________________________________________________________________
f.write('&nhs nhist=2   /')
f.write('\n')
f.write('0.100 0.001 50 10 1')
f.write('\n')
f.write('0.100 0.001 50 6 27')
f.write('\n')
f.write('&nsc nscat=0  sauto=.true. /')
f.write('\n')
f.write('&nzh nzhist=2 /')
f.write('\n')
f.write('2  0. 0.22857 70 0. 0. 33')
f.write('\n')
f.write('2  10. 0.37143 70 0. 0. 33')
f.write('\n')
f.write('&nrh nrhist=0 /')
f.write('\n')
f.write('&nem nemit=10  pxycorr=.true./')
f.write('\n')
f.write('3 6 9 12 15 18 21 24 27 30')
f.write('\n')
f.write('&ncv ncovar=3 /')
f.write('\n')
f.write('1 17 27')
f.write('\n')
#____________________________________________________________________

f.write('SECTION')
f.write('\n')
f.write('\n')


#Begin Cell structure_________________________________________
f.write('CELL')
f.write('\n')
f.write('1')
f.write('\n')
f.write('.FALSE.')
f.write('\n')
f.write('SOL')
f.write('\n')
f.write('2 %f %f %f 1  %f %f 0 0 0  0 0 0 0 0' %(popt[0], popt[1], popt[2], popt[3], popt[4]))
f.write('\n')
f.write('\n')


#Define_SREGION_of_end_length___________________________________
f.write('REPEAT')
f.write('\n')
f.write('%d' %(popt[2]/step))
f.write('\n')
f.write('SREGION')
f.write('\n')
f.write('%f 1 1e-3' %(step))
f.write('\n')
f.write('1 0. &apertr')
f.write('\n')
f.write('NONE')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0  0 0 0 0 0')
f.write('\n')
f.write('%s' % FMaterial)
f.write('\n')
f.write('CBLOCK')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0')
f.write('\n')
#Set output for the for009.dat file____________________________
f.write('OUTPUT')
f.write('\n')
#______________________________________________________________
f.write('ENDREPEAT')
f.write('\n')
f.write('\n')
#______________________________________________________________
#______________________________________________________________
#______________________________________________________________
#Define_SREGION_of_solenoid_enttrence__________________________
f.write('REPEAT')
f.write('\n')
f.write('%d' %(popt[3]/step))
f.write('\n')
f.write('SREGION')
f.write('\n')
f.write('%f 1 1e-3' %(step))
f.write('\n')
f.write('1 0. &radius')
f.write('\n')
f.write('NONE')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0  0 0 0 0 0')
f.write('\n')
f.write('%s' % FMaterial)
f.write('\n')
f.write('CBLOCK')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0')
f.write('\n')
#Set output for the for009.dat file____________________________
f.write('OUTPUT')
f.write('\n')
#______________________________________________________________
f.write('ENDREPEAT')
f.write('\n')
f.write('\n')
#______________________________________________________________
#______________________________________________________________
#______________________________________________________________
#Define_SREGION_of_absorber_posistion__________________________
f.write('REPEAT')
f.write('\n')
f.write('%d' %((popt[1]-2*popt[3])/step))
f.write('\n')
f.write('SREGION')
f.write('\n')
f.write('%f 1 1e-3' %(step))
f.write('\n')
f.write('1 0. &radius')
f.write('\n')
f.write('NONE')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0  0 0 0 0 0')
f.write('\n')
f.write('%s' % AbsMaterial)
f.write('\n')
f.write('CBLOCK')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0')
f.write('\n')
#Set output for the for009.dat file____________________________
f.write('OUTPUT')
f.write('\n')
#______________________________________________________________
f.write('ENDREPEAT')
f.write('\n')
f.write('\n')
#______________________________________________________________
#______________________________________________________________
#______________________________________________________________
#Define_SREGION_of_solenoid_entrence___________________________
f.write('REPEAT')
f.write('\n')
f.write('%d' %(popt[3]/step))
f.write('\n')
f.write('SREGION')
f.write('\n')
f.write('%f 1 1e-3' %(step))
f.write('\n')
f.write('1 0. &radius')
f.write('\n')
f.write('NONE')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0  0 0 0 0 0')
f.write('\n')
f.write('%s' % FMaterial)
f.write('\n')
f.write('CBLOCK')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0')
f.write('\n')
#Set output for the for009.dat file____________________________
f.write('OUTPUT')
f.write('\n')
#______________________________________________________________
f.write('ENDREPEAT')
f.write('\n')
f.write('\n')
#______________________________________________________________
#______________________________________________________________
#______________________________________________________________
#Define_SREGION_of_end_length___________________________________
f.write('REPEAT')
f.write('\n')
f.write('%d' %(popt[2]/step))
f.write('\n')
f.write('SREGION')
f.write('\n')
f.write('%f 1 1e-3' %(step))
f.write('\n')
f.write('1 0. &apertr')
f.write('\n')
f.write('NONE')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0  0 0 0 0 0')
f.write('\n')
f.write('%s' % FMaterial)
f.write('\n')
f.write('CBLOCK')
f.write('\n')
f.write('0 0 0 0 0  0 0 0 0 0')
f.write('\n')
#Set output for the for009.dat file____________________________
f.write('OUTPUT')
f.write('\n')
#______________________________________________________________
f.write('ENDREPEAT')
f.write('\n')

#End cell structure____________________________________________
f.write('\n')
f.write('ENDCELL')
f.write('\n')
f.write('\n')
#______________________________________________________________


f.write('ENDSECTION')
f.close()


plt.show()
