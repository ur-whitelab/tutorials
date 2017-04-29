import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

colvar_no_eds=np.genfromtxt('COLVAR_RDC_NO_EDS')
#print colvar_no_eds
colvars_w_eds=np.genfromtxt('COLVAR_RDC_EDS')
#print colvars_w_eds
bias=-4*np.ones(len(colvars_w_eds[:,0]))
plt.plot(colvar_no_eds[:,0],colvar_no_eds[:,1],label='No Bias')
plt.plot(colvars_w_eds[:,0],colvars_w_eds[:,1],label='With Bias')
plt.plot(colvars_w_eds[:,0],bias,label='Reference RDC')
plt.xlabel('Time (ps)')
plt.ylabel('Residual Diploar Coupling')
plt.legend()
plt.savefig('distance_woteds_weds_rdc')
