# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 10:47:03 2019

@author: rundaji
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#post_processing
def plot_mach(case_no,p):
    f = open('data\\task_2_c\\bump%d_curved_p%d\\bump%d_curved_triangulation.txt' %(case_no,p,case_no), "r");
    nElem_single, nn = [int(string) for string in f.readline().split()];
    tri = [None] * nElem_single;
    for i in range(0,nElem_single):
        tri[i] = [int(string) for string in f.readline().split()]; 
    f.close();
    
    f = open('data\\task_2_c\\bump%d_curved_p%d\\bump%d_curved_Mach_number.txt' %(case_no,p,case_no), "r");
    nElemTot, ns = [int(string) for string in f.readline().split()];
    x = [None] * ns;
    y = [None] * ns;
    M = [None] * ns;
    
    fig = plt.figure(figsize=(15,5));
    for i in range(0,nElemTot):
        for j in range(0,ns):
            x[j], y[j], M[j] = [float(string) for string in f.readline().split()];
        plt.tricontourf(x, y, tri, M, cmap=plt.cm.jet, vmin=0.40, vmax=0.8);
    
    plt.axis('equal');
    plt.title('Mach number, bump%d_curved, p = %d' %(case_no,p));
    
    ax1 = fig.add_axes([0.91, 0.12, 0.025, 0.76])
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=0.40, vmax=0.8)
    mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
    plt.savefig('figure\\bump%d_curved_p%d_mach.pdf' %(case_no,p), dpi=150);
    plt.close(fig);
    return 0;

def plot_history(case_no,p):
    f = open('data\\task_2_c\\bump%d_curved_p%d\\bump%d_curved_history.txt' %(case_no,p,case_no), "r");
    nTime_step = int(f.readline());
    L_inf = [None]*nTime_step;
    for i in range(0,nTime_step):
        L_inf[i] = float(f.readline());
    fig = plt.figure(figsize=(8,6));
    plt.semilogy(L_inf,'k-');
    plt.xlabel('Iteration',fontsize =16);
    plt.ylabel('$L_{\infty}$ error',fontsize =16);
    plt.grid();
    plt.title('Residual norm history, bump%d' %case_no);
    plt.savefig('figure\\bump%d_curved_p%d_L_inf_error.pdf' %(case_no,p), dpi=150);
    plt.close(fig);
    return 0;

def plot_cp(case_no,p):
    f = open('data\\task_2_c\\bump%d_curved_p%d\\bump%d_curved_cp_distribution.txt' %(case_no,p,case_no), "r");
    data = f.readlines();
    global cp;
    cp = [];
    for i in range(0,len(data)):
        cp.append([float(string) for string in data[i].split()]);
    cp = sorted(cp, key=lambda my_list: my_list[0]);
    cp = np.asarray(cp);
    x = cp[:,0];
    y = cp[:,1];
    f1 = plt.figure(figsize=(15,5));
    plt.plot(x,-y,'k-');
    plt.xlabel("x",fontsize =16);
    plt.ylabel('Pressure coefficient $-c_p$',fontsize =16);
    plt.grid();
    plt.title('Pressure coefficient distribution, bump%d' %case_no);
    plt.savefig('figure\\bump%d_curved_p%d_cp_distribution.pdf' %(case_no,p), dpi=150);
    plt.close(f1);
    return 0;

def main():
    global mesh;
    for i in range(0,3):
        for p in range(0,3):
            plot_mach(i,p);
            plot_history(i,p);
            plot_cp(i,p);
    return 0;

if __name__=="__main__":
    main()