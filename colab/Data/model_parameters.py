from    fenics import *
import  numpy as np
import ufl
import dolfin

# model parameters (non dim)
G1    = 1          # elastic modulus solig
G2    = 0.03       # elastic modulus liquid
G3    = 0.03       # elastic modulus air

#gamma1  = 0.06     # surface tension solig
#gamma2  = 0.03     # surface tension liquid
#gamma3  = 2.24     # surface tension air 

gamma1 = 0.21
gamma2 = 0.51
gamma3 = 0.41

m       = 1e-6     # Cahn-Hilliard mobility
mu      = 1        # Solid Stokes viscosity 

h0      = 1        # substrate height 
h1      = 0.36     # fluid height

H1      = 2

eps     = 0.01   # interface width

interface_factor =  Constant(1.0/sqrt(2))

def GshearF(psi1,psi2):
    psi3 = -1-psi1-psi2
    xi1 = (1+psi1)/2
    xi2 = (1+psi2)/2
    xi3 = (1+psi3)/2
    return  (xi1*G1 + xi2*G2 + xi3*G3)
