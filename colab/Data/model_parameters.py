# model parameters (non dim)
G1    = 1          # elastic modulus solig
G2    = 0.03       # elastic modulus liquid
G3    = 0.03       # elastic modulus air

gamma1  = 0.06     # surface tension solig
gamma2  = 0.03     # surface tension liquid
gamma3  = 2.24     # surface tension air 

m       = 1e-6     # Cahn-Hilliard mobility
mu      = 1        # Solid Stokes viscosity 

h0T   = 1          # substrate height 
h1T   = 0.36       # fluid height

eps   = 0.5e-2     # interface width

interface_factor =  Constant(1.0/sqrt(2))


