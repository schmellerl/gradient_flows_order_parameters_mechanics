def refine_mesh_pureDef(mesh,W,S,q,eps,h0,h1,interface_factor,H1):
    
    ep = 0.12*4
    
    v,  pp1= split(q)

    psi11 = Expression(("tanh(ff*(h0-x[1])/eps)"),eps=eps,h0=h0,h1=h1,ff=interface_factor,degree=2)
    psi22 = Expression(("-1 + 2*(1+tanh(ff*(x[1]/h0-h0)/eps))*(1+tanh(ff/eps*(h1-pow(pow((x[0])/h0-LL/2,2.0)+pow(x[1]/h0-h0,2.0),0.5))))/4"),eps=eps,h0=h0,h1=h1,LL=H1,ff=interface_factor,degree=2)
    
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
    
    TH = MixedElement([P2,R])  
    W = FunctionSpace(mesh, TH)
    S = FunctionSpace(mesh, P1)
    
    # Boundary
    bc = [DirichletBC(W.sub(0),  Constant((0, 0)), 'on_boundary')]

    return bc, W, S, mesh 
