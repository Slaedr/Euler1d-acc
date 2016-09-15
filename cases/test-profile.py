#! /usr/bin/env python
import sys
import numpy as np

fname = "test-profile-99.dat"
m = 1.0;
c = 1.0;
N = 99;
L = 1.0;

pi = 3.14159265359;
dx = L/(float)(N)
x = np.linspace(0.0+dx/2.0,L-dx/2.0,num=N)
S = x + np.ones(N)

np.savetxt(fname,S)
