## all length in m, N, J... SI units

t       = 0
T       = 0.0001
n_steps = 100
dt      = T/(n_steps)

eps     =  6e-9
Mob     =  0.001  # Cahn-Hilliard mobility
mu      =  0.1    # solid mobility /viscosity

G1 = 2.8e4
G2 = 0.3*G1      # E-modulus phase 2
G3 = 0.3*G1      # E-modulus phase 3

# in N/m
sigma1  = 1e-3  #  1.7e-3 #-7.9e-3 #1.7e-3 # 1e-3   # 36.1e-3    # solig
sigma2  = 6e-3  #  0.3e-3 #10.9e-3 #0.3e-3 # 6e-3  # 3.9e-3     # liquid
sigma3  = 24e-3 #  39e-3  #28.4e-3 #39e-3  # 24e-3   # 15.9e-3    # air 

interface_factor = Constant(1.0/sqrt(2))

            
# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True


h0    = 4.0e-6        # 50
h1    = 1.8e-6 #1.87e-6 # das korrespondiert jetzt zu 0.2 unten in den Anfangsdaten... # 1.4e-6 #1.38e-6        # m  

Nx    = 2  
Ny    = 4 

H0    = 0
H1    = 8e-6  # Breite
H2    = 7e-6  #550 # Höhe

H0    = 0/h0
H1    = H1/h0 # Breite
H2    = H2/h0 # Höhe

LL    =  H1 # H1 # 0

noslip    = Constant((0, 0)) 
noslip1d  = Constant(0) 
zero      = Constant(0)

BT      = Constant(0)

mesh = RectangleMesh(Point(H0,H0),Point(H1,H2),Nx,Ny,diagonal="right")

h0T  = h0/h0
h1T  = h1/h0
eps  = eps/h0

P2 = VectorElement("P", mesh.ufl_cell(), 1)
P1 = FiniteElement("P", mesh.ufl_cell(), 1)
R  = FiniteElement("P", mesh.ufl_cell(), 1)     # incompressible


# dfv,dfpsi,dfeta,dv,dpsi,deta
TH = MixedElement([P2,P1,P1,P2,P1,P1,R]) 
W = FunctionSpace(mesh, TH)
S = FunctionSpace(mesh, P1)



initial = Expression(("0",
                      "0",
                      "tanh(ff*(h0-x[1])/eps)",
                      "-1 + 2*(1+tanh(ff*(x[1]/h0-h0)/eps))*(1+tanh(ff/eps*(h1-pow(pow((x[0])/h0-LL/2,2.0)+pow(x[1]/h0-h0+BT,2.0),0.5))))/4",
                      "0","0",
                      "0","0","0"), eps=eps,h0=h0T,h1=h1T,degree=2,ff=interface_factor,LL=LL,BT=BT)

def GshearF(psi1,psi2):
    psi3 = -1-psi1-psi2
    xi1 = (1+psi1)/2
    xi2 = (1+psi2)/2
    xi3 = (1+psi3)/2
    return  (xi1*G1 + xi2*G2 + xi3*G3)

def boundary_bot1(x, on_boundary1):
    tol = 1E-14
    return on_boundary1 and near(x[0], H0, tol)

def boundary_bot2(x, on_boundary2):
    tol = 1E-14
    return on_boundary2 and near(x[1], H0, tol)

def boundary_bot3(x, on_boundary3):
    tol = 1E-14
    return on_boundary3 and near(x[0], H1, tol)

def boundary_bot4(x, on_boundary4):
    tol = 1E-14
    return on_boundary4 and near(x[1], H2, tol)

 
bc1   = DirichletBC(W.sub(0).sub(0), noslip1d, boundary_bot1)
bc2   = DirichletBC(W.sub(0), noslip, boundary_bot2)
bc3   = DirichletBC(W.sub(0).sub(0), noslip1d, boundary_bot3)
bc4   = DirichletBC(W.sub(0), noslip, boundary_bot4)



# Output routine
def output(q, t=0):
    v,psi1,psi2,fv,feta1,feta2 = q.split()
    
    tmp = Function(W)
    ALE.move(mesh,q.sub(0))
  
    tmp.assign(-q)
    ALE.move(mesh,tmp.sub(0))


bc = [bc1,bc2,bc3,bc4]
