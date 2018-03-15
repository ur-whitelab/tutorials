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
    plt.tight_layout()

def plot_rdf(rdfs, labels, ax):
    '''Plot the RDF'''

    for t,l in zip(rdfs, labels):
        data = np.genfromtxt(t, skip_header=4)
        ax.plot(data[:,1], data[:,2], label=l)
    ax.legend(loc='best')
    ax.set_xlabel('$r$ [Å]')
    ax.set_ylabel('$g(r)$')
    plt.tight_layout()


def plot_temp(eners, labels, ax, indices=[0,3]):
    '''plot the temperature over time'''
    for e, l in zip(eners, labels):
        data = np.genfromtxt(e, invalid_raise=False)
        ax.plot(data[:,indices[0]], data[:,indices[1]], label=l)
    ax.legend(loc='best')
    ax.set_ylim(250,350)



def plot_cv(cns, labels, set_points):
    '''plt the coordination numbers'''
    N = len(set_points)
    fig, axes = plt.subplots(N,sharex=True, sharey=True, figsize=(5,2 + 1.5 * (N - 1)))
    try:
        iter(axes)
    except TypeError:
        axes = [axes]
    for i,(e,l) in enumerate(zip(cns, labels)):
        data = np.genfromtxt(e, skip_footer=1, invalid_raise=False)
        for j,ax in enumerate(axes):
            ax.plot(data[:,0] / 1000, data[:,j+1] / set_points[j], label=l, color='C{}'.format(i))
            #ax.set_ylim(0,2)
    for p,ax in zip(set_points, axes):
        ax.axhline(y=1.0, linestyle='--', color='blue', label='set-point')
        #ax.get_yaxis().set_ticks([0,0.5,1,1.5])
    axes[0].legend(loc='best')
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)


def plot_bias(eds_files,labels, N, M=9):
    '''plot the biases'''

    #we will plot every other row to avoid to 0s in the mean
    fig, axes = plt.subplots(N,M, sharex=True, figsize=(5 + 3.5 * (M - 1),2 + 1.5 * (N - 1)))

    for i, (e, l) in enumerate(zip(eds_files, labels)):
        with open(e) as f:
            fields_str = f.readline()[2:-1]
            column_labels = fields_str.split()[1:]
            column_dict = {k: i for i,k in enumerate(column_labels)}
        data = np.genfromtxt(e, skip_footer=1, invalid_raise=False)
        for (j,k), ax in np.ndenumerate(axes):
            try:
                #maybe we aren't biasing all moments
                name = column_labels[j * M + k + 1]
            except IndexError:
                continue
            col = data[:, j * M + k + 1]#add one to skip time
            #for means, remove 0s
            if(name.split('_')[-1] == 'mean' or name.split('_')[-1] == 'std'):
                ax.plot(data[col != 0, column_dict['time']] / 1000, col[col != 0], label=l, color='C{}'.format(i))
            else:
                ax.plot(data[:, column_dict['time']] / 1000, col, label=l, color='C{}'.format(i))
            if j != N - 1:
                plt.setp(ax.get_xticklabels(), visible=False)
            if j == 0:
                ax.set_title(name.split('_')[-1])


    axes[0,0].legend(loc='best')
    fig.subplots_adjust(hspace=0)


def main():
    import sys, os
    names = sys.argv[1:]
    print('Analyzing ', ','.join(names))
    eners = ['{}.ener'.format(n) for n in names]
    cns = ['{}_cn.log'.format(n) for n in names]
    eds = ['{}_eds.restart'.format(n) for n in names[1:]]
    trajs = ['{}_pos.xyz'.format(n) for n in names]
    rdfs = ['{}.rdf'.format(n) for n in names]

    #read in set-points
    set_points = []
    with open('set_points.txt', 'r') as f:
        for line in f.readlines():
            if(not line.startswith('#')):
                set_points.append(float(line[:-1]))

    mpl.style.use(['seaborn-muted', 'seaborn-white'])
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

    #skip first (unbiased)
    plot_bias(eds, names[1:], len(set_points))
    plt.savefig('biases_comparison.png')



if __name__ == '__main__':
    main()
