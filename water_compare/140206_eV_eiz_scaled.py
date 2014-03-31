#!/usr/bin/python
import matplotlib               
#import pyreport
import numpy as np                  
from pylab import *
#from pylab import show
from matplotlib import pyplot as pl

x_x,y_x = np.loadtxt('data/CNT6_5_xe2_solid_30.txt',unpack=True, usecols = [0,1])
x_z,y_z = np.loadtxt('data/CNT6_5_ze2_solid_30.txt',unpack=True, usecols = [0,1])
x_w,y_w      = np.loadtxt('data/water-L.txt',unpack=True, usecols = [0,1])
z, RR_eiz_w = np.loadtxt('data/LD_water_pars1.txt', delimiter='     ',usecols=(0,1), unpack=True)

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#------------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s
#coeff = 2.41e14 # in (1 rad)*(1/s)=inverse seconds
T = 297.0
#kb_J = 1.3806488e-23 # in J/K
#hbar = 6.625e-34 # in J/s
#coeff_J = 2.0*np.pi*kb_J*T/hbar#1.602e-19*0.159e15 # in eV #(2.41*1e14) # in rad/s
#n = arange(0,500)
#z = n * coeff
#coeff_J = 1.602e-19*0.159e15 # in eV #(2.41*1e14) # in rad/s

#z = n * coeff
#z = n * coeff_J

eiz_x = empty(len(z))
eiz_y = empty(len(z))
eiz_z = empty(len(z))
eiz_w = empty(len(z))

eiz_x_arg=empty(len(x_x))
eiz_z_arg=empty(len(x_z))
eiz_w_arg=empty(len(x_w))

for j in range(len(z)):
    for i in range(len(x_x)):
        eiz_x_arg[i]=x_x[i]*y_x[i] / (x_x[i]**2 + z[j]**2)
    eiz_x[j] = 1 + (2./pi) * trapz(eiz_x_arg,x_x)

    for m in range(len(x_z)):
        eiz_z_arg[m]=x_z[m]*y_z[m] / (x_z[m]**2 + z[j]**2)
    eiz_z[j] = 1 + (2./pi) * trapz(eiz_z_arg,x_z)    

    for p in range(len(x_w)):
        eiz_w_arg[p]=x_w[p]*y_w[p] / (x_w[p]**2 + z[j]**2)
    eiz_w[j] = 1 + (2./pi) * trapz(eiz_w_arg,x_w)    
#savetxt("data/eiz_x_output.txt", eiz_x)
#savetxt("data/eiz_y_output.txt", eiz_y)
#savetxt("data/eiz_z_output.txt", eiz_z)
#savetxt("data/eiz_w_output.txt", eiz_w)
#
savetxt("data/eiz_x_rr_output_eV.txt", eiz_x)
savetxt("data/eiz_z_rr_output_eV.txt", eiz_z)
savetxt("data/eiz_w_rr_L_output_eV.txt", eiz_w)

#
pl.figure()
pl.plot(z,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
pl.plot(z,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{n})$')
pl.plot(z,eiz_w, color = 'c', label = r'$\varepsilon_{\hat{w}}(i\zeta_{n})$')
pl.plot(z,RR_eiz_w, color = 'm', label = r'$\varepsilon_{RR}(i\zeta_{n})$')
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 24)
pl.legend()
pl.title(r'[6,5] and water eiz')
#pl.savefig('plots/DNA_eiz_x_z.png', dpi = 300 )
pl.savefig('plots/65w65_eiz_with_no_RR.pdf')
pl.show()
#pl.close()


pl.figure()
pl.loglog(z,eiz_x, color = 'b', label = r'$\varepsilon_{\hat{x}}(i\zeta_{N})$')
pl.loglog(z,eiz_z, color = 'r', label = r'$\varepsilon_{\hat{z}}(i\zeta_{n})$')
pl.loglog(z,eiz_w, color = 'c', label = r'$\varepsilon_{\hat{w}}(i\zeta_{n})$')
pl.loglog(z,RR_eiz_w, color = 'm', label = r'$\varepsilon_{RR}(i\zeta_{n})$')
pl.xlabel(r'$N$', size = 24)
pl.ylabel(r'$\varepsilon(i\zeta)$', size = 24)
pl.legend()
pl.title(r'[6,5] and water eiz')
#pl.savefig('plots/DNA_eiz_x_z.png', dpi = 300 )
pl.savefig('plots/loglog_65w65_no_RR_eiz.pdf')
pl.show()
