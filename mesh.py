# mesh
mesh0 = RectangleMesh(Point(P.L,P.L),Point(P.H,P.H),P.Nx,P.Ny,diagonal="right")

domain1 = Circle(Point(0, 0), 0.5)
mesh1   = generate_mesh(domain1-Circle(Point(0,-0.2),0.2),32)

mesh = mesh0


# Define the function space
P2 = VectorElement("P", mesh.ufl_cell(), 1)
P1 = FiniteElement("P", mesh.ufl_cell(), 1)
R  = FiniteElement("Real", mesh.ufl_cell(), 0)

# dfv,dfpsi,dfeta,dv,dpsi,deta
TH = MixedElement([P2,P1,P1,R,R,R,R]) 

W = FunctionSpace(mesh, TH)
X = FunctionSpace(mesh, P2)

xx = interpolate(Expression(('x[0]','x[1]'),degree=2),X)
