from IPython.display import HTML, display
from fenics import *
from mshr import *

def progress(value, max=100):
    return HTML("""
        Iterations: {value} of {max}
        <progress
            value='{value}'
            max='{max}',
            style='width: 100%'
        >
            {value}
        </progress>
    """.format(value=value, max=max))

def output_mesh(mesh,filepath):
    File(filepath + "mesh.xml") << mesh

def output_mesh1(mesh1,filepath):
    File(filepath + "mesh1.xml") << mesh1

def output_solution(q,n,mesh,P2,S,psi1,psi2,filepath,t=0):
    File(filepath + "u"+str(n)+".xml")    << project(q.sub(0),FunctionSpace(mesh,P2)) 
    File(filepath + "psi1"+str(n)+".xml") << project(psi1, S)
    File(filepath + "psi2"+str(n)+".xml") << project(psi2, S)
