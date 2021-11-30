# 1 phase free boundary incompressible  

t       = 0
T       = 0.1e-2 # 2d T und s_steps proportional hoch für dieses eps = 0.02, T = 0.3e-3, n_steps = 100
n_steps = 1000
dt      = T/(n_steps)

# Constants
Gshear  = Constant(1.)
Mob     = Constant(1.)
eps     = Constant(0.02)
alpha   = Constant(0.4)


L     = Constant(1)
H     = Constant(4)
 
parameters["form_compiler"]["cpp_optimize"] = True


# Boundary conditions, forcing terms, and initial conditions
noslip    = Constant((0, 0)) 
noslip1d  = Constant(0) 
zero      = Constant(0)


Nx      = 100
Ny      = 100

# mesh 
mesh0 = RectangleMesh(Point(L,L),Point(H,H),Nx,Ny,diagonal="right")

domain1 = Circle(Point(0, 0), 0.5)
mesh1   = generate_mesh(domain1-Circle(Point(0,-0.2),0.2),32)

mesh = mesh1


# Define the function space
P2 = VectorElement("P", mesh.ufl_cell(), 1)
P1 = FiniteElement("P", mesh.ufl_cell(), 1)
R  = FiniteElement("Real", mesh.ufl_cell(), 0)

# dfv,dfpsi,dfeta,dv,dpsi,deta
TH = MixedElement([P2,P1,P1,R,R,R,R]) 

W = FunctionSpace(mesh, TH)
X = FunctionSpace(mesh, P2)

xx = interpolate(Expression(('x[0]','x[1]'),degree=2),X)
