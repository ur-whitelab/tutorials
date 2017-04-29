import numpy as np
import matplotlib.pyplot as plt
def mov_aver(data,window):
    results=np.ones(len(data))
    results[0]=data[0]
    c=2/float(window+1)
    for i in np.arange(len(results)): 
        results[i]=(1-c)*results[i-1]+c*data[i]
    return results


colvar_no_eds=np.genfromtxt('colvars_no_eds.dat')
#print colvar_no_eds
colvars_w_eds=np.genfromtxt('colvars.dat')
#print colvars_w_eds
bias=0.7*np.ones(len(colvars_w_eds[:,0]))
plt.plot(colvar_no_eds[:,0],colvar_no_eds[:,1],label='No Bias')
plt.plot(colvars_w_eds[:,0],colvars_w_eds[:,1],label='With Bias')
plt.plot(colvars_w_eds[:,0],bias,label='Reference')
plt.xlabel('Time (ps)')
plt.ylabel('Distance (nm)')
plt.legend()
plt.savefig('distance_eds_noeds')
