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


def output_mesh(mesh):
    File(filepath + "mesh.xml") << mesh

def output_mesh1(mesh1):
    File(filepath + "mesh1.xml") << mesh1

def output_solution(q,n,t=0):
    u, psi1, psi2, *_ = q.split()

    File(filepath + "u"+str(n)+".xml")      << project(q.sub(0),FunctionSpace(mesh,Vu)) 
    File(filepath + "psi1"+str(n)+".xml")   << project(q.sub(1),FunctionSpace(mesh,Vpsi))
    File(filepath + "psi2"+str(n)+".xml")   << project(q.sub(2),FunctionSpace(mesh,Vpsi))
    File(filepath + "solid"+str(n)+".xml")  << project(q.sub(0), V_solid)
