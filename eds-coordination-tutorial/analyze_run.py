#-*- coding: utf-8 -*-


import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import MDAnalysis, MDAnalysis.analysis.rdf
import numpy as np

def plot_traj_rdf(trajs, labels, ax):
    '''Plot the RDF'''
    
    for t,l in zip(trajs, labels):
        u = MDAnalysis.Universe(t)
        g = u.select_atoms('name O')
        rdf = MDAnalysis.analysis.rdf.InterRDF(g, g)
        rdf.run()
        print(rdf.rdf)
        ax.plot(rdf.bins, rdf.rdf, label=l)
        
    ax.legend(loc='best')
    ax.set_xlabel('$r$ [Å]')
    ax.set_ylabel('$g(r)$')

def plot_rdf(rdfs, labels, ax):
    '''Plot the RDF'''
    
    for t,l in zip(rdfs, labels):
        data = np.genfromtxt(t, skip_header=4)
        ax.plot(data[:,1], data[:,2], label=l)        
    ax.legend(loc='best')
    ax.set_xlabel('$r$ [Å]')
    ax.set_ylabel('$g(r)$')
    

def plot_temp(eners, labels, ax, indices=[0,3]):
    '''plot the temperature over time'''
    for e, l in zip(eners, labels):
        data = np.genfromtxt(e)
        ax.plot(data[:,indices[0]], data[:,indices[1]], label=l)
    ax.legend(loc='best')



def plot_cv(cns, labels, set_points):
    '''plt the coordination numbers'''    
    N = len(cns)
    fig, axes = plt.subplots(N,sharex=True, sharey=True, figsize=(5,4 + 2.5 * (N - 1)))
    try:
        iter(axes)
    except TypeError:
        axes = [axes]
    cmap=plt.cm.get_cmap('Accent', N)
    for i,(e,l) in enumerate(zip(cns, labels)):
        data = np.genfromtxt(e, skip_footer=1)
        for j,ax in enumerate(axes):
            ax.plot(data[:,0], data[:,i+1] / set_points[i], label=l, color=cmap(i / N))
    axes[0].legend(loc='best')
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        


def main():
    import sys, os
    names = sys.argv[1:]
    print('Analyzing ', ','.join(names))
    eners = ['{}.ener'.format(n) for n in names]
    cns = ['{}_cn.log'.format(n) for n in names]
    trajs = ['{}_pos.xyz'.format(n) for n in names]
    rdfs = ['{}.rdf'.format(n) for n in names]

    #read in set-points
    with open('set_points.txt', 'r') as f:
        set_points = [float(x[:-1]) for x in f.readlines()]

    mpl.style.use('seaborn-muted')
    mpl.rcParams['savefig.dpi'] = 300


    plt.figure(figsize=(4,3))
    ax = plt.gca()
    ax.set_title('RDF')

    #the reference as well
    data = np.genfromtxt('reference.rdf')
    ax.plot(data[:,0], data[:,1], label='Experiment [Skinner et al.]')
    ax.set_xlim(0,10)

    if os.path.exists(rdfs[0]):        
        plot_rdf(rdfs, names, ax)
    elif os.path.exists(trajs[0]):
        plot_traj_rdf(trajs, names, ax)
    else:
        raise NameError('Could not find trajs or rdfs')
    
    plt.savefig('rdf_comparison.png')
                
    plt.figure(figsize=(4,3))
    ax = plt.gca()    
    ax.set_title('Temperature')
    plot_temp(eners, names, ax)
    plt.savefig('temperature_comparison.png')

    plot_cv(cns, names, set_points)
    plt.savefig('cn_comparison.png')
    
    
if __name__ == '__main__':
    main()
