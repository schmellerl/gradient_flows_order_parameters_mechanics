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
     
    Vu    = VectorElement("P", mesh.ufl_cell(), 1) # V_u   = displacements
    Vpsi  = FiniteElement("P", mesh.ufl_cell(), 1) # V_psi = order-parameter(s)
    Uu    = VectorElement("P", mesh.ufl_cell(), 1)
    Upsi  = FiniteElement("P", mesh.ufl_cell(), 1) # U     = forces
    R     = FiniteElement("P", mesh.ufl_cell(), 1) 

    TH    = MixedElement([Vu,Vpsi,Vpsi,Vpsi,Vpsi,R]) # real raum
     
    VxUxR = FunctionSpace(mesh, TH)

    bc = [DirichletBC(VxUxR.sub(0),  Constant((0, 0)), 'on_boundary')]
    
    return bc, VxUxR, Vpsi, mesh
