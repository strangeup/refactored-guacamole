//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#include <fenv.h> 
//Generic routines
#include "generic.h" 

// The equations
#include "C1_large_displacement_plate_models.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
// The coupling of the stretching energy
double nu = 0.3;
double h = 0.01;
double eta_u = 1;
double eta_sigma = 1;
double eta_bending = h;
double Theta = Pi / 10. ;

/// Different nondimensionalisations
enum Nondimensional_form {none = 0 , oomph = 1 /*, thickness_scaled= 2*/ };
Nondimensional_form nondimensional_form = none;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// The normal and tangential directions.
void get_normal_and_tangent(const Vector<double>& x, Vector<double>& n, 
 Vector<double>& t, DenseMatrix<double>& Dn, DenseMatrix<double>& Dt)
{
 // Fill in the normal and derivatives of the normal
 n[0] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
 n[1] = x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
 t[0] =-x[1]/(x[0]*x[0]+x[1]*x[1]);
 t[1] = x[0]/(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dt(0,0)=(+2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(0,1)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,0)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,1)=(-2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
}

// Exact solution for constant pressure, circular domain and resting boundary 
// conditions
void get_exact_w(const Vector<double>& xi, Vector<double>& w)
{
 const double /*x = xi[0],*/ y = xi[1];
 w = Vector<double>(18,0.0);

 w[6] =  (cos(Theta)-1)*y;
 w[8] =  (cos(Theta)-1);
  
 w[12] = sin(Theta)*y;
 w[14] = sin(Theta);
}

// Exact solution for constant pressure, circular domain and resting boundary 
// conditions
void get_exact_w_radial(const Vector<double>& xi, Vector<double>& w)
{
 const double x = xi[0], y = xi[1];
 w = Vector<double>(18,0.0);

 const double A = (cos(Theta)-1);
 w[6+0] =  A*y; // A r sin(t)
 w[6+1] =  A*(y / sqrt(x*x+y*y)); // A sin(t)
 w[6+2] =  A*x; // A r cos(t)
 w[6+3] =  0.0;
 w[6+4] = A*(x / sqrt(x*x+y*y)); // A r sin(t)
 w[6+5] =-A*y;
  
 const double B = (sin(Theta));
 w[12+0] =  B*y; // B r sin(t)
 w[12+1] =  B*(y / sqrt(x*x+y*y)); // B sin(t)
 w[12+2] =  B*x; // B r cos(t)
 w[12+3] =  0.0;
 w[12+4] = B*(x / sqrt(x*x+y*y)); // B r sin(t)
 w[12+5] =-B*y;
}

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

/// Constructor
UnstructuredFvKProblem(double element_area = 0.1);

/// Destructor
~UnstructuredFvKProblem()
{
 delete (Surface_mesh_pt);
 delete (Bulk_mesh_pt);
 
};

/// Update after solve (empty)
void actions_after_newton_solve()
{
}

/// Update the problem specs before solve: Re-apply boundary conditions
/// Empty as the boundary conditions stay fixed
void actions_before_newton_solve()
{
apply_boundary_conditions();
}

/// Set the initial values to the exact solution (useful for debugging)
void set_initial_values_to_exact_solution();

/// Doc the solution
void doc_solution(const std::string& comment="");

/// \short Overloaded version of the problem's access function to
/// the mesh. Recasts the pointer to the base Mesh object to
/// the actual mesh type.
TriangleMesh<ELEMENT>* mesh_pt()
{
return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()); 
}

/// Doc info object for labeling output
DocInfo Doc_info;

private:
/// Helper function to apply boundary conditions
void apply_boundary_conditions();

/// \short Helper function to (re-)set boundary condition
/// and complete the build of  all elements
void complete_problem_setup();

/// Trace file to document norm of solution
ofstream Trace_file;

// Keep track of boundary ids
enum
{
 Outer_boundary0 = 0,
 Outer_boundary1 = 1,
 Inner_boundary0 = 2
};

double Element_area;

/// \short Create bending moment elements on the b-th boundary of the
/// problems mesh 
void create_traction_elements(const unsigned &b, Mesh* const & bulk_mesh_py,
                            Mesh* const &surface_mesh_pt);

/// Helper to upgrade the edge elements to curved Bernadou elements
void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

/// Helper to rotate the edge elements' Hermite dofs
void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area)
{
Vector<double> zeta(1);
Vector<double> posn(2);

//Outer boundary
//--------------

double A = 1.0;
double B = 1.0;
Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);

TriangleMeshClosedCurve* outer_boundary_pt = 0;

Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi;
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

outer_boundary_pt =
new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

// Internal bit - this means we can have a boundary which is just the centre
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
vertices[0][0] = 1.0;

unsigned boundary_id = Inner_boundary0;

TriangleMeshPolyLine *boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);

// Total number of open curves in the domain
unsigned n_open_curves = 1;
// We want internal open curves
Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);
// Connect it
boundary2_pt -> connect_initial_vertex_to_curviline(dynamic_cast<TriangleMeshCurviLine*>(outer_curvilinear_boundary_pt[0]),0.0);

// Each internal open curve is defined by a vector of
// TriangleMeshCurveSection,
// on this example we only need one curve section for each internal boundary
Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
internal_curve_section1_pt[0] = boundary2_pt;

// The open curve that define this boundary is composed of just one
// curve section
 inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(internal_curve_section1_pt);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(outer_boundary_pt);

mesh_parameters.element_area() = element_area;

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

// Build an assign bulk mesh
Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

// Create "surface mesh" that will contain only the prescribed-traction
// elements. The constructor creates the mesh without adding any nodes
// elements etc.
Surface_mesh_pt =  new Mesh;

//Add two submeshes to problem
add_sub_mesh(Bulk_mesh_pt);
add_sub_mesh(Surface_mesh_pt);

// Combine submeshes into a single Mesh
build_global_mesh();

// Curved Edge upgrade
upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
 
// Rotate degrees of freedom
rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';

delete outer_boundary_pt;
delete outer_boundary_ellipse_pt;
delete outer_curvilinear_boundary_pt[0];
delete outer_curvilinear_boundary_pt[1];
}


//==start_of_complete======================================================
// Set the initial values to the exact solution (useful for debugging)
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::set_initial_values_to_exact_solution()
{
// Loop over all elements and reset the values
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
 // Upcast from GeneralisedElement to the present element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 // Get the numbers
 const unsigned nnode = el_pt->nnode();
 const unsigned ndof_type = el_pt->nnodal_position_type();
 const unsigned ndisplacement = 3;
 const unsigned dim = 2;
 // Loop over nodes
 for(unsigned i=0;i<nnode;++i)
  {
   // Containers
   Vector<double> x(dim),w(ndisplacement*ndof_type,0.0);
   // Get node pointer
   Node* nod_pt = el_pt->node_pt(i);
   // Get x 
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get test
   TestSoln::get_exact_w(x,w);
   // Set value
   for(unsigned k=0;k<ndisplacement;++k)
    {
    for(unsigned l=0 ; l<ndof_type;++l)
     {
      const unsigned index = el_pt->u_index_koiter_model(l,k);
      nod_pt->set_value(index,w[ndof_type*k+l]);
     }
    }
  }
  // Pin the internal dofs
  const unsigned n_internal_dofs = el_pt->number_of_internal_dofs();
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   // Containers
   Vector<double> x(dim,0.0), s(dim,0.0), w(ndisplacement*ndof_type,0.0);
   // Get internal data
   Data* internal_data_pt = el_pt->internal_data_pt(1);
   el_pt->get_internal_dofs_location(i,s);
   // Get test
   el_pt->get_coordinate_x(s,x); 
   TestSoln::get_exact_w_radial(x,w);
   // HERE need a function that can give us this lookup
   for(unsigned k=0; k<3; ++k)
    {
     // The bubble index
     const unsigned index =i+k*n_internal_dofs;
     internal_data_pt->set_value(index,w[ndof_type*k+0]);
    }
   }  

  //Just loop over the boundary elements
  unsigned nbound = Outer_boundary1 + 1;
  for(unsigned b=0;b<nbound;b++)
   {
   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nb_element;e++)
    {
     // Get pointer to bulk element adjacent to b
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     // Loop over vertices of the element (i={0,1,2} always!)
     for(unsigned i=0;i<nnode;++i)
      {
       Node* nod_pt = el_pt->node_pt(i);
       // If it is on the bth boundary
       if(nod_pt -> is_on_boundary(b))
        {
         // Set the values of the unknowns that aren't pinned
         Vector<double> x(dim,0.0),w(ndof_type*ndisplacement,0.0);
         x[0] = nod_pt->x(0);
         x[1] = nod_pt->x(1);
         TestSoln::get_exact_w_radial(x,w);
         // Just set the nonzero values -> the others are pinned
         // Set value
         for(unsigned k=0;k<ndisplacement;++k)
          {
          for(unsigned l=0 ; l<ndof_type;++l)
           {
            const unsigned index = el_pt->u_index_koiter_model(l,k);
            nod_pt->set_value(index,w[ndof_type*k+l]);
           }
          }
        }
      }
    }
   }
}// End loop over elements
}


//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of 
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{   
unsigned nbound = Outer_boundary1 + 1;
 // Upcast to current element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
//Just loop over outer boundary since inner boundary doesn't have boundary
//conditions
for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get nod
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 Vector<double> x(2,0.0),w(18,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);

 TestSoln::get_exact_w_radial(x,w);

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=0;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal
    if(l != 3)  
     {
     unsigned index=el_pt->u_index_koiter_model(l,i);
     nod_pt->pin(index);
     nod_pt->set_value(index,w[index]);
     }
    }
  }
 }
} // end loop over boundaries 


// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
// Upcast from GeneralisedElement to the present element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

//Set the pressure function pointers and the physical constants
el_pt->nu_pt() = &TestSoln::nu;
el_pt->thickness_pt() = &TestSoln::eta_bending;
el_pt->eta_sigma_pt() = &TestSoln::eta_sigma;
}

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
// apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
 unsigned num_nod=Bulk_mesh_pt->nnode();
 for(unsigned inod=0;inod<num_nod;++inod)
 {
 Node* nod_pt=Bulk_mesh_pt->node_pt(inod);

 // Extract nodal coordinates from node: Vector<double> x(2),w(18);
 Vector<double> x(2),w(18,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);
 // If is at centre
 if(x[0] == 0.0 && x[1] == 0.0)
 {
 // Pin the displacements
 TestSoln::get_exact_w(x,w);
 // Pin Rotation
 nod_pt->pin(12+1);
 nod_pt->set_value(12+1,w[12+1]);
 nod_pt->pin(12+2);
 nod_pt->set_value(12+2,w[12+2]);
 // Pin displacement
 const unsigned n_disp_dofs = 3*6;
 for(unsigned index=0; index<n_disp_dofs;index+=6)
  {
   nod_pt->pin(index+0);
   nod_pt->set_value(index+0,w[index+0]);
  }
 }
 }
} // end set bc

/// A function that upgrades straight sided elements to be curved. This involves
// Setting up the parametric boundary, F(s) and the first derivative F'(s)
// We also need to set the edge number of the upgraded element and the positions
// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
// representation. That is we need to have a continuous 2nd derivative defined 
// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky , 
// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s) 
// as well.
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &bulk_mesh_pt) 
{
 // How many bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt-> nboundary_element(b);
 // These depend on the boundary we are on
 CurvilineGeomObject* parametric_curve_pt; 
 // Define the functions for each part of the boundary
 switch (b)
  {
   // Upper boundary
   case 0:
    parametric_curve_pt = &TestSoln::parametric_curve_top;
   break;

   // Lower boundary
   case 1:
    parametric_curve_pt = &TestSoln::parametric_curve_bottom;
   break;

   default:
    throw OomphLibError(
     "I have encountered a boundary number that I wasn't expecting. This is very\
 peculiar.",
     "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
     OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Loop over nodes
   const unsigned nnode=3;
   unsigned index_of_interior_node=3;

   // The edge that is curved
   MyC1CurvedElements::Edge edge;

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     if(nod_pt->is_on_boundary(0) || nod_pt->is_on_boundary(1))
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar
   double s_ubar, s_obar;

   // s at the next (cyclic) node after interior
   s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
   // s at the previous (cyclic) node before interior
   s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);

   // Assign edge case
   switch(index_of_interior_node)
    {
     case 0: edge= MyC1CurvedElements::zero; 
      break;
     case 1: edge= MyC1CurvedElements::one; 
      break;
     case 2: edge= MyC1CurvedElements::two; 
      break;
     // Should break it here HERE
     default: edge= MyC1CurvedElements::none; 
      throw OomphLibError(
       "The edge number has been set to a value greater than two: either we have\
 quadrilateral elements or more likely the index_of_interior_node was never set\
 and remains at its default value.",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
      break;
     }
   if (s_ubar>s_obar)
    {
     oomph_info <<"Apparent clockwise direction of parametric coordinate."
                <<"This will probably result in an inverted element."
                <<"s_start="<<s_ubar<<"; s_end ="<<s_obar<<std::endl;
     throw OomphLibError(
       "The Edge coordinate appears to be decreasing from s_start to s_end. \
Either the parametric boundary is defined to be clockwise (a no-no) or \
the mesh has returned an inverted element (less likely)",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
    }

   // Upgrade it
   bulk_el_pt->upgrade_to_curved_element(edge,s_ubar,s_obar,
    parametric_curve_pt);     
  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
   // Loop nodes  
   unsigned nnode = el_pt->nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary(0) || 
                       el_pt->node_pt(n)->is_on_boundary(1));}

   // Now if we have nodes on boundary 
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // If it is on the boundary
       if(el_pt->node_pt(n)->is_on_boundary(0) || 
                       el_pt->node_pt(n)->is_on_boundary(1))
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
      }
    // Output that we have found element HERE
    std::cout<<"Element "<<e<<" has "<<bnode<< " nodes on the boundary.\n";

    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
   // Now rotate the nodes
   }
 }
}// end

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                    std::string& comment)
{ 
ofstream some_file;
char filename[100];

// Number of plot points
unsigned npts = 2;

sprintf(filename,"RESLT/coarse_soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Number of plot points
npts = 6;

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/exact_interpolated_soln%i-%f.dat","RESLT",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_exact_w); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Output boundaries
//------------------
sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_boundaries(some_file);
some_file.close();

// Output regions
unsigned n_region = Bulk_mesh_pt->nregion();
if (n_region > 1)
{
for (unsigned r = 0; r < n_region; r++)
{
 //Attempt to output elements in different regions
 sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
   Element_area);
 some_file.open(filename);
 unsigned nel = Bulk_mesh_pt->nregion_element(r);
 for (unsigned e = 0; e < nel; e++)
  {
   Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
  }
 some_file.close();
}
}

// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Error of computed solution: " << sqrt(dummy_error)<< std::endl;
 oomph_info << "Norm of computed solution: "  << sqrt(zero_norm)  << std::endl;

 // Find the solution at r=0
 //   // ----------------------
 // Interestingly the number of bins is miscalculated for element_area = 1
 // whcih causes a segfault
 // Discuss?
 MeshAsGeomObject* Mesh_as_geom_obj_pt=
  new MeshAsGeomObject(Bulk_mesh_pt);
 Vector<double> s(2);
 GeomObject* geom_obj_pt=0;
 Vector<double> r(2,0.0);
 r[1]=0.5;
 Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);
 // Call the interpolated_u function
 Vector<Vector<double> > u_i(3,Vector<double>(6,0.0));
 dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_koiter_plate(s,u_i);

 oomph_info << "r at (0.0,0.5): "  
            <<"("<< 0.0+u_i[0][0] <<","<<0.5+u_i[1][0]<<","<<u_i[2][0]<<")"<< std::endl;
 
 Trace_file << TestSoln::Theta<<" "<< u_i[0][0] <<" "<<u_i[1][0]<<" "<<u_i[2][0]<< std::endl;

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2\n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm;
some_file.close();

// Increment the doc_info number
Doc_info.number()++;

// Clean up
delete Mesh_as_geom_obj_pt;
} // end of doc


//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Validation
 CommandLineArgs::specify_command_line_flag("--validate");

 // Directory for solution
 string output_dir="RSLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 CommandLineArgs::specify_command_line_flag("--theta", &TestSoln::Theta);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,5,KoiterSteigmannPlateEquations> >problem(element_area);
 // Change some tolerances
 problem.max_residuals()=1e10;
 problem.newton_solver_tolerance()=1e-12;
 problem.max_newton_iterations()=30;
 problem.set_initial_values_to_exact_solution();
 problem.newton_solve();
 problem.doc_solution();
} //End of main

