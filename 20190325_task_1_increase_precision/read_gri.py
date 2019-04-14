# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 12:20:59 2019

@author: rundaji
"""

from my_class import CELL, BGROUP, EGROUP
import numpy as np

def read_gri(fname):
    f = open(fname, "r");
    #read some general info
    nNode,nElemTot,Dim = [int(string) for string in f.readline().split()];
    node_pos = [None]*nNode;
    #read the position of nodes
    for i in range(0,nNode):
        x,y = [float(string) for string in f.readline().split()];
        node_pos[i] = [x,y];
    node_pos =  np.asarray(node_pos);
    #read the number of boundary groups
    nBGroup = int(f.readline());
    #read the boundaries
    BGroup = [None]*nBGroup;
    for i in range(0,nBGroup):
        nBFace,nf,Title = f.readline().split();
        nBFace = int(nBFace);
        nf = int(nf);
        B = [None]*nBFace;
        n = [None]*nf;
        for j in range(0,nBFace):
            n = [int(string) for string in f.readline().split()];
            # the given index start from 1, we want the index start from 0
            for k in range(0,nf):
                n[k] -= 1;
            # n.sort(); do not sort! we need the first and last element to determine the adj cell!  
            B[j] = n;
        BGroup[i] = BGROUP(nBFace, nf, Title, B);
    #read cell info
    curTot = 0;
    nEGroup = 0;
    EGroup = [];
    E_global = [None]*nElemTot;
    while curTot != nElemTot:
        nElem,Order,Basis = f.readline().split();
        nElem = int(nElem);
        Order = int(Order);
        q = Order;
        nq = int(0.5*(q + 1)*(q + 2));
        E_local = [None]*nElem;
        for i in range(0,nElem):
            n = [int(string) for string in f.readline().split()];
            # the given index start from 1, we want the index start from 0
            for j in range(0,nq):
                n[j] -= 1;
            if q == 1:
                nodes_edge_0 = [n[1],n[2]];
                nodes_edge_1 = [n[2],n[0]];
                nodes_edge_2 = [n[0],n[1]];
#            elif q == 2:
#                nodes_edge_0 = [n[2],n[4],n[5]];
#                nodes_edge_1 = [n[5],n[3],n[0]];
#                nodes_edge_2 = [n[0],n[1],n[2]];
            elif q == 3:
                nodes_edge_0 = [n[3],n[6],n[8],n[9]];
                nodes_edge_1 = [n[9],n[7],n[4],n[0]];
                nodes_edge_2 = [n[0],n[1],n[2],n[3]];
            else:
                print("Error: q is too high!\n");
            e0 = {'nodes': nodes_edge_0, 't': i, 'e':0};
            e1 = {'nodes': nodes_edge_1, 't': i, 'e':1};
            e2 = {'nodes': nodes_edge_2, 't': i, 'e':2};
            edges = [e0, e1, e2];
            # nodes, edges, no, adj_cell, state, R, dt, A
            E_global[curTot] = CELL(n, edges, i, [], None, None, None, None);
            E_local[i] = curTot;
            curTot += 1;
        EGroup.append(EGROUP(nElem, Order, Basis, E_local));
        nEGroup += 1;
    
    mesh = {'nNode':nNode, 'node_pos':node_pos,'nBGroup':nBGroup, 'BGroup':BGroup, 'nElemTot':nElemTot, 'nEGroup':nEGroup, 'EGroup':EGroup, 'Elems':E_global, 'nEdge': 0, 'Edges':[]};
    f.close();
    return mesh;