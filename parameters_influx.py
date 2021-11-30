from fenics import *
from mshr import *

## Time
t       = 0
T       = 0.1e2
n_steps = 100
dt      = T/(n_steps)

# Constants
Gshear  = Constant(0.01)
kB      = Constant(1.)
T       = Constant(1.)
vg      = Constant(1)
Mob     = Constant(0.5)
chi     = Constant(1.) 
eps     = Constant(0.001) 
B       = Constant((0.01))

K = kB*T/vg

## Mesh
H = Constant(1)
L = Constant(5)
N = Constant(0.)
N = Constant(0.)

Nx = 50
Ny = 10

mesh = RectangleMesh(Point(N,N),Point(L,H),Nx,Ny,diagonal="right")

# Define the function space
P2 = VectorElement("P", mesh.ufl_cell(), 1)
P1 = FiniteElement("P", mesh.ufl_cell(), 1)
R   = FiniteElement("Real", mesh.ufl_cell(), 0)     # incompressible

V  = FunctionSpace(mesh, P2)
S  = FunctionSpace(mesh, P1)
W  = FunctionSpace(mesh, MixedElement([P2,P1,P1,R]))

## Boundary
def boundary_bot1(x, on_boundary1):
    tol = 1E-14
    return on_boundary1 and near(x[0], N, tol)

def boundary_bot2(x, on_boundary2):
    tol = 1E-14
    return on_boundary2 and near(x[1], N, tol)

def boundary_bot3(x, on_boundary3):
    tol = 1E-10
    return on_boundary3 and near(x[0], L, tol)

def boundary_bot4(x, on_boundary4):
    tol = 1E-14
    return on_boundary4 and near(x[1], H, tol)


def up(x, on_boundary):
        tol = 1E-12
        return on_boundary and near(H, x[1], tol)


noslip    = Constant((0, 0)) 
noslip1d  = Constant(0) 


bc1   = DirichletBC(W.sub(0).sub(0), noslip1d, boundary_bot1)
bc2   = DirichletBC(W.sub(0), noslip, boundary_bot2)
bc3   = DirichletBC(W.sub(0).sub(0), noslip1d, boundary_bot3)
#bc4   = DirichletBC(W.sub(0), noslip, boundary_bot4)

bc = [bc1,bc2,bc3] 
