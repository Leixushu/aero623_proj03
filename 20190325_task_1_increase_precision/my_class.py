# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 21:13:12 2018

@author: Runda Ji
"""

class EDGE():
    def __init__(self, nodes, t1, e1, t2, e2, norm_vec, length, s):
        self.nodes = nodes;
        self.t1 = t1;
        self.e1 = e1;
        self.t2 = t2;
        self.e2 = e2;
        self.norm_vec = norm_vec;
        self.length = length;
        self.s = s;

class CELL():
    def __init__(self, nodes, edges, no, adj_cell, state, R, dt, A):
        self.nodes = nodes;
        self.edges = edges;
        self.no = no; #the global index of triangle (cell)
        self.adj_cell = adj_cell;
        self.state = state;
        self.R = R;
        self.dt = dt;
        self.A = A;
        
class BGROUP():
    def __init__(self, nBFace, nf, Title, B):
            self.nBFace = nBFace;
            self.nf = nf;
            self.Title = Title;
            self.B = B;
            
class EGROUP():
    def __init__(self, nElem, Order, Basis, E):
            self.nElem = nElem;
            self.Order = Order;
            self.Basis = Basis;
            self.E = E;



