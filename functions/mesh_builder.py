# mesh_builder.py

from dolfin import *
from mshr import *
from numpy import sqrt,cos,sin,double
from mpi4py import MPI

def Fiber_Waveguide_Mesh(r_core,
                        r_clad,
                        single     = True,
                        d_core     = 0.99,
                        ratio_core = 0.02,
                        ratio_clad = 0.01):
    '''

    This is a function to create a circular mesh with a cladding radius of
    r_clad, and one of the two core configurations (toggled by s_or_d):

    1) Single Core of radius r_core. The core is declared as a subdomain
    of the cladding to make sure there are enough cells that take care of the
    discontinuity at the boundary.

    2) Two cores of radius r_core each, separated by a distance d_core.
    Variable y represents the height at which the secondary circles are made
    above the dual core for smoothening the interaction region.

    '''

    if single:

        #Creating Single core
        core = Circle(Point(0.0, 0.0), r_core, int(r_core*100))

    else:
        #Generating Dual core region
        core = Circle(Point(d_core, 0.0), r_core, int(r_core*100)) \
                + Circle(Point(-d_core, 0.0), r_core, int(r_core*100))

    # Creating cladding region
    clad = Circle(Point(0.0, 0.0), r_clad, int(r_clad*50))
    #Labelling the core as a subdomain of cladding
    clad.set_subdomain(1,core)

    #Generating Mesh
    mesh = generate_mesh(clad,1)
    #Generating markers for the core and cladding
    markers = MeshFunction('size_t', mesh, 2, mesh.domains())

    #~ # Refine Core
    mesh, markers = Proper_Refine(mesh, markers, 1, ratio_core*r_core)
    #~ #mesh, markers = Proper_Refine(mesh, markers, 2, ratio_core*r_core)

    #~ # Refine Cladding
    mesh, markers = Proper_Refine(mesh, markers, 0, ratio_clad*r_clad)

    #Returning the refined mesh and markers
    return mesh, markers

def Proper_Refine(mesh, markers, index, limit_radius):

    '''
    This function finds all the cells in subdomain labelled 'index', finds the
    cells within the subdomain with inradius() greater than the limiting radius,
    and keeps refining those cells until the limit is achieved.
    '''
    n_unrefined = 2 # Tag to count the number of unrefined cells

    while (n_unrefined > 0): # loop breaks till all cells are refined

        n_unrefined  = 0 # set to zero
        cell_markers = MeshFunction("bool", mesh, 2, mesh.domains())
        cell_markers.set_all(False)

        for cell in cells(mesh):

            if (markers[cell] == index) and (cell.inradius() >= limit_radius):

                cell_markers[cell] = True
                n_unrefined += 1 # increment n_unrefined if unrefined cell is found

        mesh    = refine(mesh,cell_markers)
        markers = adapt(markers, mesh)

    return mesh, markers


def refine_maxvalue(cf, mesh, expr, threshold):

    W = FunctionSpace(mesh, dolfin.FiniteElement("DG", mesh.ufl_cell(), 0))
    u = project(expr.magnitude, W)

    uv = u.vector()
    dofmap = W.dofmap()
    for cell in dolfin.cells(mesh):
        cell_index = cell.index()
        cell_dofs = dofmap.cell_dofs(cell_index)
        val = np.amax(np.abs(uv[cell_dofs]))
        if (val >= threshold):
            cf[cell] = True

    return cf
