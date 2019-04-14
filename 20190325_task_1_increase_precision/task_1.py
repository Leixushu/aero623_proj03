# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 23:38:47 2019

@author: runda
"""

from my_class import EGROUP
from read_gri import read_gri
from pre_calculation import pre_calculation
import numpy as np
import matplotlib.pyplot as plt

def curve(mesh):
    global new_nodes,new_edges;
    [new_nodes,elem_flag] = split_edge(mesh);
    update_elem(mesh,new_nodes,elem_flag);
    return 0;

def split_edge(mesh):
    edge_flag = np.zeros(mesh['nEdge']);
    elem_flag = np.zeros(mesh['nElemTot']);
    for i in range(0,mesh['nElemTot']):
        for j in range(0,3):
            global_edge_no = mesh['Elems'][i].edges[j];
            edge = mesh['Edges'][global_edge_no];
            if edge.e2 == 'Bottom':
                elem_flag[i] = 1;
        if elem_flag[i] == 1:
            for j in range(0,3):
                global_edge_no = mesh['Elems'][i].edges[j];
                edge_flag[global_edge_no] = 1;
                
    for i in range(0,mesh['nBGroup']):
        if mesh['BGroup'][i].Title == "Bottom":
            mesh['BGroup'][i].nf = 4; # 4 nodes per face on the bottom wall
    #--------------------------------------------------------------------------
    new_nodes = [];
    for i in range(0,mesh['nEdge']):
        if edge_flag[i] == 1:
            edge = mesh['Edges'][i];
            n0 = edge.nodes[0];
            n1 = edge.nodes[1];
            x0,y0 = mesh['node_pos'][n0];
            x1,y1 = mesh['node_pos'][n1];
            trisect_pt_pos_0 = np.array([x0 + (x1-x0)*1.0/3.0, y0 + (y1-y0)*1.0/3.0]);
            trisect_pt_pos_1 = np.array([x0 + (x1-x0)*2.0/3.0, y0 + (y1-y0)*2.0/3.0]);
            #correct the trisect_pt position for lower boundary
            if edge.e2 == 'Bottom':
                trisect_pt_pos_0[1] = 0.0625*np.exp(-25*trisect_pt_pos_0[0]**2);
                trisect_pt_pos_1[1] = 0.0625*np.exp(-25*trisect_pt_pos_1[0]**2);
                #add the trisect_pt to current list
            mesh['node_pos'] = np.append(mesh['node_pos'], [trisect_pt_pos_0], axis = 0);
            new_node_no_0 = mesh['nNode'];
            mesh['nNode'] = mesh['nNode'] + 1;
            new_nodes.append(new_node_no_0);
            mesh['node_pos'] = np.append(mesh['node_pos'], [trisect_pt_pos_1], axis = 0);
            new_node_no_1 = mesh['nNode'];
            mesh['nNode'] = mesh['nNode'] + 1;
            new_nodes.append(new_node_no_1);
            #split the edge            
            mesh['Edges'][i].nodes = [n0,new_node_no_0 ,new_node_no_1,n1];
    return new_nodes, elem_flag;

def update_elem(mesh,new_nodes,elem_flag):
    for i in range(0,mesh['nElemTot']):
        if elem_flag[i] == 1:
            elem = mesh['Elems'][i];
            # add centeroid
            x = np.zeros(3);
            y = np.zeros(3);
            for j in range(0,3):
                global_node_n0 = elem.nodes[j];
                x[j],y[j] = mesh['node_pos'][global_node_n0];
            center_pos = np.array([(x[0]+x[1]+x[2])/3.0,(y[0]+y[1]+y[2])/3.0]);
            mesh['node_pos'] = np.append(mesh['node_pos'], [center_pos], axis = 0);
            new_node_no = mesh['nNode'];
            mesh['nNode'] = mesh['nNode'] + 1;
            n = [None]*10;
            if mesh['Edges'][elem.edges[0]].nodes[0] == elem.nodes[1] and mesh['Edges'][elem.edges[0]].nodes[-1] == elem.nodes[2]:
                n[3],n[6],n[8],n[9] = mesh['Edges'][elem.edges[0]].nodes;
            else: #swap the direction
                n[9],n[8],n[6],n[3] = mesh['Edges'][elem.edges[0]].nodes;
            if mesh['Edges'][elem.edges[1]].nodes[0] == elem.nodes[2] and mesh['Edges'][elem.edges[1]].nodes[-1] == elem.nodes[0]:
                n[9],n[7],n[4],n[0] = mesh['Edges'][elem.edges[1]].nodes;
            else:
                n[0],n[4],n[7],n[9] = mesh['Edges'][elem.edges[1]].nodes;
            if mesh['Edges'][elem.edges[2]].nodes[0] == elem.nodes[0] and mesh['Edges'][elem.edges[2]].nodes[-1] == elem.nodes[1]:
                n[0],n[1],n[2],n[3] = mesh['Edges'][elem.edges[2]].nodes;
            else:
                n[3],n[2],n[1],n[0] = mesh['Edges'][elem.edges[2]].nodes;
            n[5] = new_node_no;
            mesh['Elems'][i].nodes = n;
    #--------------------------------------------------------------------------
    mesh['nEGroup'] = 2;
    mesh['EGroup'] = [None]*mesh['nEGroup'];
    mesh['EGroup'][0] = EGROUP(0, 1, 'TriLagrange', []);
    mesh['EGroup'][1] = EGROUP(0, 3, 'TriLagrange', []);
    for i in range(0,mesh['nEGroup']):
        mesh['EGroup'][i].nElem = 0;
        mesh['EGroup'][i].Basis = 'TriLagrange';
        mesh['EGroup'][i].E = [];
    mesh['EGroup'][0].Order = 1;
    mesh['EGroup'][1].Order = 3;
    for i in range(0,mesh['nElemTot']):
        if elem_flag[i] == 0:
            mesh['EGroup'][0].E.append(i);
            mesh['EGroup'][0].nElem += 1;
        else:
            mesh['EGroup'][1].E.append(i);
            mesh['EGroup'][1].nElem += 1;
    return 0;

def write_gri_file(mesh,folder,title):
    file = open('%s\\%s.gri' %(folder,title), 'w');
    file.write('%d %d %d\n' %(mesh['nNode'],mesh['nElemTot'],2));
    for i in range(0,mesh['nNode']):
        pos = mesh['node_pos'][i];
        file.write('%.16f %.16f\n' %(pos[0], pos[1]));
    file.write('%d\n' %(mesh['nBGroup']));
    for i in range(0,mesh['nBGroup']):
        BGroup = mesh['BGroup'][i];
        file.write('%d %d %s\n' %(BGroup.nBFace, BGroup.nf, BGroup.Title));
        for j in range(0,BGroup.nBFace):
            global_edge_no = BGroup.B[j];
            for k in range(0,BGroup.nf):
                n = mesh['Edges'][global_edge_no].nodes[k] + 1; #the index start from 1
                file.write('%d ' %n);
            file.write('\n');
    for i in range(0,mesh['nEGroup']):
        EGroup = mesh['EGroup'][i];
        file.write('%d %d %s\n' %(EGroup.nElem,EGroup.Order,EGroup.Basis));
        q = EGroup.Order;
        nq = int(0.5*(q + 1)*(q + 2));
        for j in range(0,EGroup.nElem):
            global_elem_no = EGroup.E[j];
            elem = mesh['Elems'][global_elem_no];
            for k in range(0,nq):
                n = elem.nodes[k] + 1;
                file.write('%d ' %n);
            file.write('\n');
    file.close();
    return 0;

def plot_mesh(mesh,folder,title):
    f1 = plt.figure(figsize=([15,4]));
    for i in range(0,mesh['nEdge']):
        Edge = mesh['Edges'][i];
        n = Edge.nodes;
        plot_exterior(mesh,n);
    for i in range(0,mesh['nEGroup']):
        EGroup = mesh['EGroup'][i];
        for j in range(0,EGroup.nElem):
            global_elem_no = EGroup.E[j];
            n = mesh['Elems'][global_elem_no].nodes;
            
            if EGroup.Order == 3:
                plot_interior(mesh,n);
    plt.axis('equal');
    plt.savefig('%s\\%s.pdf' %(folder,title),dpi=150);
    plt.xlim(-0.25, 0.25);
    plt.ylim(-0, 0.2);
    plt.savefig('%s\\%s_zoomed.pdf' %(folder,title),dpi=150);
    plt.close(f1);
    return 0;   

def plot_exterior(mesh,n):
    x = [None]*len(n);
    y = [None]*len(n);
    for i in range(0,len(n)):
        x[i],y[i] = mesh['node_pos'][n[i]];
    plt.plot(x,y,'k-',lw=0.5);

def plot_interior(mesh,n):
    n_order = [1,5,2,6,5,8,7,5,4,1];
    x = [None]*len(n_order);
    y = [None]*len(n_order);
    for i in range(0,len(n_order)):
        local_node_no = n_order[i];
        global_node_no = n[local_node_no];
        x[i],y[i] = mesh['node_pos'][global_node_no];
    plt.plot(x,y,'b-',lw=0.5);
    return 0;

def main():
    global mesh,R_max_bump;
    Ncurve = 4;
    for i in range(0,Ncurve):
        mesh = read_gri('mesh\\bump%d.gri' %i);
        pre_calculation(mesh);
        plot_mesh(mesh,'figure','mesh_bump%d'%i);
        curve(mesh);
        write_gri_file(mesh,'mesh','bump%d_curved'%i);
        plot_mesh(mesh,'figure','mesh_bump%d_curved'%i);
    return 0;

if __name__=="__main__":
    main()