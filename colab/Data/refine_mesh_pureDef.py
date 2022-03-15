from    fenics import *
import  numpy as np
import ufl
import dolfin

def refine_mesh_pureDef(mesh,W,S,q,eps,h0,h1,interface_factor,H1):
    
    ep = 0.12*4
    
    v,  pp1= split(q)

    psi11 = Expression(("tanh(ff*(h0-x[1])/eps)"),eps=eps,h0=h0,h1=h1,ff=interface_factor,degree=2)
    psi22 = Expression(("-1 + 2*(1+tanh(ff*(x[1]/h0-h0)/eps))*(1+tanh(ff/eps*(h1-pow(pow((x[0])/h0-H1/2,2.0)+pow(x[1]/h0-h0,2.0),0.5))))/4"),eps=eps,h0=h0,h1=h1,H1=H1,ff=interface_factor,degree=2)
    
    psi1 = interpolate(psi11, S)
    psi2 = interpolate(psi22, S)

    tmp = project(psi1,S) 
    tmp2 = project(psi2,S)

    cell_markers = MeshFunction("bool", mesh,2)
    cell_markers.set_all(False)

    qTmp = tmp.vector()
    qTmp2 = tmp2.vector()
    
    dm = S.dofmap()
       
    for cell in cells(mesh):
       
        cell_index = cell.index()
        cell_dofs = dm.cell_dofs(cell_index)

        a = qTmp[cell_dofs].min()
        b = qTmp[cell_dofs].max()
    
        a2 = qTmp2[cell_dofs].min()
        b2 = qTmp2[cell_dofs].max()
 
        if (a-ep<0) and (b+ep>0):
            cell_markers[cell] = True
            
        if (a2-ep<0) and (b2+ep>0):
            cell_markers[cell] = True
                            
    mesh = refine(mesh, cell_markers)
    FEu_DEG = -2
    FEp_DEG = -1
    
    #P2        = VectorElement("P", mesh.ufl_cell(), 1)
    Pk      = FiniteElement("Lagrange", mesh.ufl_cell(), -FEu_DEG)
    B       = FiniteElement("Bubble",   mesh.ufl_cell(), 3)
    P2      = VectorElement(NodalEnrichedElement(Pk, B))
    
    P1        = FiniteElement("P", mesh.ufl_cell(), 1)
    #R         = FiniteElement("P", mesh.ufl_cell(), 1)    
    R       = FiniteElement("DG", mesh.ufl_cell(), -FEp_DEG)
    
    TH        = MixedElement([P2,R]) 
    W         = FunctionSpace(mesh, TH)
    S         = FunctionSpace(mesh, P1)
    
    # Boundary
    #bc = [DirichletBC(W.sub(0),  Constant((0, 0)), 'on_boundary')]
    noslip = Constant((0, 0)) 
    noslip1d = Constant((0))

    def boundary_bot1(x, on_boundary1):
        tol = 1E-14
        return on_boundary1 and near(x[0], 0, tol)

    def boundary_bot2(x, on_boundary2):
        tol = 1E-14
        return on_boundary2 and near(x[1], 0, tol)

    def boundary_bot3(x, on_boundary3):
        tol = 1E-14
        return on_boundary3 and near(x[0], 2, tol)

    def boundary_bot4(x, on_boundary4):
        tol = 1E-14
        return on_boundary4 and near(x[1], 2, tol)


    bc1   = DirichletBC(W.sub(0).sub(0),  noslip1d, boundary_bot1)
    bc2   = DirichletBC(W.sub(0),         noslip, boundary_bot2)
    bc3   = DirichletBC(W.sub(0).sub(0),  noslip1d, boundary_bot3)
    bc4   = DirichletBC(W.sub(0),         noslip, boundary_bot4)

    bc = [bc1,bc2,bc3,bc4]
    
    
    psi1 = interpolate(psi11, S)
    psi2 = interpolate(psi22, S)

    return bc, W, S, P2, mesh, psi1, psi2
