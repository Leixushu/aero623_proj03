# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 14:37:40 2019

@author: runda
"""

import numpy as np
import matplotlib.pyplot as plt

cl_exact = 1.537095;
cd_exact = 2.94278 * 10 ** (-6);
Es_exact = 0.0;
    
def find_error(fname):
    global temp;
    f = open(fname, "r");
    f.readline();
    case_no = 9;
    DOF = np.zeros(case_no);
    cl_err = np.zeros(case_no);
    cd_err = np.zeros(case_no);
    Es_err = np.zeros(case_no);
    for i in range(0,case_no):
        title,p,DOF[i],cl,cd,Es = f.readline().split();
        DOF[i] = np.sqrt(int(DOF[i]));
        cl_err[i] = abs(float(cl) - cl_exact);
        cd_err[i] = abs(float(cd) - cd_exact);
        Es_err[i] = abs(float(Es) - Es_exact);
    err = {'DOF':DOF,'e_cl':cl_err,'e_cd':cd_err,'e_Es':Es_err};
    return err;

def plot_error(err):
    f1 = plt.figure(figsize=(10,8));
    plt.loglog(err['DOF'][0:3],err['e_cl'][0:3],'r-^',  label='p = 0, $c_l$');
    plt.loglog(err['DOF'][3:6],err['e_cl'][3:6],'r--x', label='p = 1, $c_l$');
    plt.loglog(err['DOF'][6:9],err['e_cl'][6:9],'r-.s', label='p = 2, $c_l$');
    plt.legend(handlelength=5,loc='upper left',bbox_to_anchor=(1.0, 1.0));
    plt.xlabel('$\sqrt{dof}$',fontsize =20);
    plt.ylabel('Error',fontsize = 20);
    plt.grid();
    plt.tight_layout();
    plt.savefig('conv_study_cl.pdf', dpi=150);
    plt.close(f1);
    
    f2 = plt.figure(figsize=(10,8));
    plt.loglog(err['DOF'][0:3],err['e_cd'][0:3],'g-^',  label='p = 0, $c_d$');
    plt.loglog(err['DOF'][3:6],err['e_cd'][3:6],'g--x', label='p = 1, $c_d$');
    plt.loglog(err['DOF'][6:9],err['e_cd'][6:9],'g-.s', label='p = 2, $c_d$');
    plt.legend(handlelength=5,loc='upper left',bbox_to_anchor=(1.0, 1.0));
    plt.xlabel('$\sqrt{dof}$',fontsize =20);
    plt.ylabel('Error',fontsize = 20);
    plt.grid();
    plt.tight_layout();
    plt.savefig('conv_study_cd.pdf', dpi=150);
    plt.close(f2);
    
    f3 = plt.figure(figsize=(10,8));
    plt.loglog(err['DOF'][0:3],err['e_Es'][0:3],'b-^',  label='p = 0, $E_s$');
    plt.loglog(err['DOF'][3:6],err['e_Es'][3:6],'b--x', label='p = 1, $E_s$');
    plt.loglog(err['DOF'][6:9],err['e_Es'][6:9],'b-.s', label='p = 2, $E_s$');
    plt.legend(handlelength=5,loc='upper left',bbox_to_anchor=(1.0, 1.0));
    plt.xlabel('$\sqrt{dof}$',fontsize =20);
    plt.ylabel('Error',fontsize = 20);
    plt.grid();
    plt.tight_layout();
    plt.savefig('conv_study_Es.pdf', dpi=150);
    plt.close(f3);
    return 0;

def func(x, a, b):
    return a + b*x;

def find_conv_rate(err):
    print('e_cl rates');
    rate0 = np.log2(err['e_cl'][0+2]/err['e_cl'][0+1]);
    rate1 = np.log2(err['e_cl'][3+2]/err['e_cl'][3+1]);
    rate2 = np.log2(err['e_cl'][6+2]/err['e_cl'][6+1]);
    print('p=0 rate: %.2f, p=1 rate: %.2f, p=2 rate: %.2f' %(rate0,rate1,rate2));
    
    print('e_cd rates');
    rate0 = np.log2(err['e_cd'][0+2]/err['e_cd'][0+1]);
    rate1 = np.log2(err['e_cd'][3+2]/err['e_cd'][3+1]);
    rate2 = np.log2(err['e_cd'][6+1]/err['e_cd'][6+0]);
    print('p=0 rate: %.2f, p=1 rate: %.2f, p=2 rate: %.2f' %(rate0,rate1,rate2));
    
    print('e_Es rates');
    rate0 = np.log2(err['e_Es'][0+2]/err['e_Es'][0+1]);
    rate1 = np.log2(err['e_Es'][3+2]/err['e_Es'][3+1]);
    rate2 = np.log2(err['e_Es'][6+2]/err['e_Es'][6+1]);
    print('p=0 rate: %.2f, p=1 rate: %.2f, p=2 rate: %.2f' %(rate0,rate1,rate2));
    return 0;

def main():
    global err_1st, err_2nd;
    err = find_error('cl_cd_Es.txt');
    plot_error(err);
    find_conv_rate(err);
    return 0;

if __name__=="__main__":
    main()