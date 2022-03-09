from    fenics import *
import  numpy as np
import ufl
import dolfin

# model parameters (non dim)
G1    = 1          # elastic modulus solig
G2    = 0.03       # elastic modulus liquid
G3    = 0.03       # elastic modulus air

gamma1  = 0.06     # surface tension solig
gamma2  = 0.03     # surface tension liquid
gamma3  = 2.24     # surface tension air 

m       = 1e-6     # Cahn-Hilliard mobility
mu      = 1        # Solid Stokes viscosity 

h0      = 1        # substrate height 
h1      = 0.36     # fluid height

H1      = 2

eps     = 0.01   # interface width

interface_factor =  Constant(1.0/sqrt(2))


