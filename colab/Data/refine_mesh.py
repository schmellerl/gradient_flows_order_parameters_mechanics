from    fenics import *
import  numpy as np
import ufl
import dolfin

def refine_mesh(mesh,VxUxR,Vpsi,q):
    
    ep = 0.12*4
    
    v, psi1, psi2, feta1, feta2, pp1= split(q)
    
    tmp   = project(psi1,FunctionSpace(mesh,Vpsi)) 
    tmp2 = project(psi2,FunctionSpace(mesh,Vpsi))
   
    cell_markers = MeshFunction("bool", mesh,2)
    cell_markers.set_all(False)
    
    qTmp = tmp.vector()
    qTmp2 = tmp2.vector()

    dm = FunctionSpace(mesh,Vpsi).dofmap()


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
 
    #Vu      = VectorElement("P", mesh.ufl_cell(), 1) # V_u   = displacements
    Pk      = FiniteElement("Lagrange", mesh.ufl_cell(), -FEu_DEG)
    B       = FiniteElement("Bubble",   mesh.ufl_cell(), 3)
    Vu      = VectorElement(NodalEnrichedElement(Pk, B))
    
    Vpsi  = FiniteElement("P", mesh.ufl_cell(), 1) # V_psi = order-parameter(s)
    Upsi  = FiniteElement("P", mesh.ufl_cell(), 1) # U     = forces
    
    
    #R       = FiniteElement("P", mesh.ufl_cell(), 1) # Lagrange Multiplier
    R       = FiniteElement("DG", mesh.ufl_cell(), -FEp_DEG)
    
    VxUxR   = FunctionSpace(mesh, MixedElement([Vu,Vpsi,Vpsi,Upsi,Upsi,R]))
    #Vu      = FunctionSpace(mesh, Pv)

    #bc = [DirichletBC(VxUxR.sub(0),  Constant((0, 0)), 'on_boundary')]
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


    bc1   = DirichletBC(VxUxR.sub(0).sub(0),  noslip1d, boundary_bot1)
    bc2   = DirichletBC(VxUxR.sub(0),         noslip, boundary_bot2)
    bc3   = DirichletBC(VxUxR.sub(0).sub(0),  noslip1d, boundary_bot3)
    bc4   = DirichletBC(VxUxR.sub(0),         noslip, boundary_bot4)

    bc = [bc1,bc2,bc3,bc4]
    
    return bc, VxUxR, Vpsi, Vu, mesh
