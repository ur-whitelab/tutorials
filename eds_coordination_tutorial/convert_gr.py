#!/usr/bin/env python
 
 
from scipy.integrate import *
from scipy.interpolate import interp1d
from math import pi, exp, sqrt,erf
import numpy as np
 
grtypes = ["gromacs", "lammps", "aimd"]
 
def read_gr(infile, grtype, rmax):
    if(grtype not in grtypes):
        raise ValueError("grtype must be one of {}".format(grtypes))
 
    gr = []
    r = []
    if(grtype == "aimd"):
        with open(infile) as f:
            for line in f:
                if(line[0] in ("#", "@") or len(line) <= 1):
                    continue
                sline = line.split()
                gr.append(float(sline[1]))
                r.append(float(sline[0]))
 
        grtype = "lammps"
         
    if(grtype == "gromacs"):
        with open(infile) as f:
            for line in f:
                if(line[0] in ("#", "@") or len(line) <= 1):
                    continue
                sline = line.split()
                gr.append(float(sline[1]))
                r.append(float(sline[0]) / 10.)
 
    if(grtype == "lammps"):
        with open(infile) as f:
            for line in f:
                sline = line.split()
                if(line[0] in ("#") or len(sline) != 4):
                    continue
                gr.append(float(sline[2]))
                r.append(float(sline[1]))
 
    #truncate to rmax if its given
    if(rmax is not None):
        r = [x for x in r if x <= rmax]
        gr = gr[:len(r)]
    else:
        rmax = r[-1]
 
    if(r[-1] > rmax):
        r.append(rmax)
        gr.append(1.0)
 
    #make sure g(r) is correctly normalized
    grf = interp1d(r, gr)
    grapprox = lambda x: 1. if x > r[-1] else grf(x)
    z = rmax ** 3 / (3 * quad(lambda r : r ** 2 * grapprox(r), r[0], rmax, limit=4 * len(gr))[0])
 
    return r,gr
 
def gr_min(r, gr):
    """ Find the minimum after the first peak
    """
    grmax = gr.index(max(gr))
    r0 = r[gr.index(min(gr[grmax:]))]
    return r0
 
 
def main(infile, grtype, density, r0, m=12, n=6, reval=None, rmax=None, moment=0, w=None):
     
    r,gr = read_gr(infile, grtype, rmax)
     
   #append values to gr
    if(r0 is None):
        #find first minimum
        r0 = gr_min(r,gr)
 
    if(rmax is None):
        rmax = r[-1]
 
    if(w is not None):
        afxn = cn_f2_fxn(gr,r,r0,w,m,n, moment)
        efxn = cn_exact_fxn(gr,r,r0+0.5*w, moment)
    else:
        afxn = cn_approx_fxn(gr,r,r0,m,n, moment)
        efxn = cn_exact_fxn(gr,r,r0, moment)
 
    print ('r0 = {}, n = {}, m = {}, Exact: {}, Approximate: {}'.format(r0, n, m, density * efxn(r[-1]), density * afxn(r[-1])))
    if(reval is not None):
        if(w is not None):
            afxn = cn_approx2_fxn(gr,r,reval,w,m,n, moment)
            efxn = cn_exact_fxn(gr,r,r0 + 0.5*w, moment)
        else:
            afxn = cn_approx_fxn(gr,r,reval,m,n, moment)
            efxn = cn_exact_fxn(gr,r,r0, moment)
 
 
        print('Approx({}) = {}, Exact({}): {}'.format(reval, density * afxn(r[-1]), reval, density * efxn(r[-1])))
     
    #assume rmax = rbox / 2 (cubic box).
    #Get requiv, which is for a sphere but equivalent volume
#    requiv = (3 / (4 * pi) * 8 * (rmax ** 3)) ** (1 / 3.)
#    print 'Check: N(rmax) = {}, assuming cubic box'.format(density * cn_exact_fxn(gr, r, requiv)(requiv))
 
def cn_approx_fxn(gr, r, r0, m, n, moment=0):
    return lambda x : (quad(integrand(gr, r, lambda y : cn_approx_def(y, r0, m, n), moment), r[0], x, points=[r0], limit=4 * len(gr)))[0]
 
def cn_approx2_fxn(gr, r, r0, w, m, n, moment=0):
    return lambda x : (quad(integrand(gr, r, lambda y : cn_approx2_def(y, r0, w, m, n), moment), r[0], x, points=[r0], epsabs = 1e-10,limit=10 * len(gr)))[0]
 
 
def cn_exact_fxn(gr, r, r0, moment=0):
    return lambda x : (simps(*fixed_integrand(gr, r, lambda y : cn_exact_def(y, r0), x, moment), dx=1))
 
def integrand(gr, r, f, moment = 0):
    grf = interp1d(r, gr)    
    grapprox = lambda x: 1. if x > r[-1] else (0 if x <= r[0] else grf(x))
    return lambda x: 4 * pi * x**(2 + moment) * grapprox(x) * f(x)
 
def integrand3(gr, r, f, moment = 0):
    grf = interp1d(r, gr)    
    grapprox = lambda x: 1. if x > r[-1] else (0 if x <= r[0] else grf(x))
    return lambda x: 4 * pi * x**2 * grapprox(x) * f(x)**(1+moment)
 
def fixed_integrand(gr, r, f, r_limit, moment = 0):
    stop = min(1 + np.argmin([(ri - r_limit)**2 for ri in r]), len(r) - 1)
    if(stop == 0):
        return [0,0], r[:2]
    fr = [f(ri) for ri in r[:stop]]
    return 4 * pi * np.array(r[:stop])**(2 + moment) * gr[:stop] * fr, r[:stop]
 
 
def uniform_integrand(gr,r,f, moment = 0):
    return lambda x: 4 * pi * x**(2) * f(x)
 
def cn_approx_def(r, r0, m, n):
    if(r == r0):
        return 1.
    return (1 - (r/r0)**n)  / (1 - (r/r0)**m)
 
def cn_approx2_def(r, r0, w, m, n):
    if(r <= r0):
        return 1.
    if((r - r0)/ w == 1.):
        return float(n) / m
    return (1 - ((r - r0)/w)**n)  / (1 - ((r - r0)/w)**m)
 
 
def cn_plumed_def(r, r0, w, s):
    upper = ((r0 + w) - r) / (sqrt(2.0) * s)
    lower = (r0 - r) / (sqrt(2.0) * s)
    return 0.5 * (erf(upper) - erf(lower))
     
def cn_exact_def(r, r0):
    return 1. if r < r0 else 0.
 
def lj_potential(r, epsilon, sigma):
    return 4 *epsilon* ((sigma / r)**12 - (sigma / r)**6)
 
 
 
if __name__ == '__main__':
    import argparse
 
    parser = argparse.ArgumentParser(description='Given a g(r), find the coordination number of the first minimum and the approximate coordination number for biasing')
    parser.add_argument('gr_file')
    parser.add_argument('gr_file_type')
    parser.add_argument('density', type=float)
    parser.add_argument('-m', default=12, type=float)
    parser.add_argument('-n', default=6, type=float)
    parser.add_argument('-r0', default=None, type=float)
    parser.add_argument('-w', default=None, type=float, help="Coordination number width")
    parser.add_argument('-reval', default=None, type=float)
    parser.add_argument('-rmax', default=None, type=float)
    parser.add_argument('-moment', default=0, type=float, help="Set this to use higher order moments of the g(r) distribution")
    pargs = parser.parse_args()
    main(pargs.gr_file, pargs.gr_file_type, pargs.density, pargs.r0, m=pargs.m, n=pargs.n, reval=pargs.reval, rmax=pargs.rmax, moment=pargs.moment, w=pargs.w)
