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
double eta = 0;
double p_mag = 0.; 
double radius_of_curvature = Pi; 
const double nu = 0.3;
double h = 0.01;
const double eta_xy = 1.0;//h*h;
double eta_z = 1.0;//h;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Assigns the value of pressure depending on the position (x,y)
inline void get_pressure(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& d_ui_dx_j,  const Vector<double>& ni, 
 Vector<double>& pressure)
{
 // Create the matrix of tangents
 DenseMatrix<double> g = d_ui_dx_j;
 g(0,0) += 1.0;
 g(1,1) += 1.0;
 pressure[0] =p_mag*(-(g(1,1)*g(2,0)) + g(1,0)*g(2,1));
 pressure[1] =p_mag*(  g(0,1)*g(2,0)  - g(0,0)*g(2,1));
 pressure[2] =p_mag*(-(g(0,1)*g(1,0)) + g(0,0)*g(1,1));
}

// Assigns the value of pressure depending on the position (x,y)
inline void get_d_pressure_dn(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj, const Vector<double>& ni, 
 DenseMatrix<double>& d_pressure_dn)
{
}

// Assigns the value of pressure depending on the position (x,y)
void get_d_pressure_dr(const Vector<double>& xi,const Vector<double>& ui,
  const DenseMatrix<double>& dui_dxj,
 const Vector<double>& ni, DenseMatrix<double>& d_pressure_dn)
{
 
}

// Assigns the value of pressure depending on the position (x,y)
inline void get_d_pressure_d_grad_u(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& d_ui_dx_j,
 const Vector<double>& ni, RankThreeTensor<double>& d_pressure_du_grad)
{
 // Create the matrix of tangents
 DenseMatrix<double> g = d_ui_dx_j;
 g(0,0) += 1.0;
 g(1,1) += 1.0;

 // Fill in the pressure components 
 d_pressure_du_grad(0,1,1) =p_mag*(-g(2,0));
 d_pressure_du_grad(0,2,0) =p_mag*(-g(1,1));
 d_pressure_du_grad(0,1,0) =p_mag*( g(2,1));
 d_pressure_du_grad(0,2,1) =p_mag*( g(1,0));

 d_pressure_du_grad(1,0,1) =p_mag*(  g(2,0));
 d_pressure_du_grad(1,2,0) =p_mag*(  g(0,1));
 d_pressure_du_grad(1,0,0) =p_mag*(- g(2,1));
 d_pressure_du_grad(1,2,1) =p_mag*(- g(0,0));

 d_pressure_du_grad(2,0,1) =p_mag*(-g(1,0));
 d_pressure_du_grad(2,1,0) =p_mag*(-g(0,1));
 d_pressure_du_grad(2,0,0) =p_mag*( g(1,1));
 d_pressure_du_grad(2,1,1) =p_mag*( g(0,0));
}

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

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w(const Vector<double>& xi, Vector<double>& w)
{
 const double x = xi[0], y = xi[1];
w[0] =  sin(radius_of_curvature*x)/radius_of_curvature-x;
w[1] =  cos(radius_of_curvature*x)-1.0;
w[2] =  0.0;
w[3] = -sin(radius_of_curvature*x)*radius_of_curvature;
w[4] =  0.0;
w[5] =  0.0;

w[6] =  0.0;
w[7] =  0.0;
w[8] =  0.0;
w[9] =  0.0;
w[10] = 0.0;
w[11] = 0.0;
 
w[12] = (cos(radius_of_curvature*x)-1.0)/radius_of_curvature;
w[13] =-sin(radius_of_curvature*x);
w[14] = 0.0;
w[15] =-cos(radius_of_curvature*x)*radius_of_curvature;
w[16] = 0.0;
w[17] = 0.0;
}

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w_radial(const Vector<double>& xi, Vector<double>& w)
{
 const double x = xi[0], y = xi[1], K=radius_of_curvature;
 const double r = sqrt(x*x+y*y), t = atan2(y,x),r2 = (x*x+y*y);
w[0] =  sin(radius_of_curvature*x)/radius_of_curvature-x;
w[1] =  (cos(K*x)-1)*x/r;
w[2] = -(cos(K*x)-1)*y;
w[3] = -(K*sin(K*x))*x*x/r2;
w[4] =  (K*sin(K*x))*x*y/r - (cos(K*x)-1)*y/r;
w[5] = -(K*sin(K*x))*y*y - (cos(K*x)-1)*x ;

w[6] =  0.0;
w[7] =  0.0;
w[8] =  0.0;
w[9] =  0.0;
w[10] = 0.0;
w[11] = 0.0;
 
w[12] =  (cos(radius_of_curvature*x)-1.0)/radius_of_curvature;
w[13] = -(sin(K*x))*x/r;
w[14] = +(sin(K*x))*y;
w[15] = -(K*cos(K*x))*x*x/r2;
w[16] =  (K*cos(K*x))*x*y/r + (sin(K*x))*y/r;
w[17] = -(K*cos(K*x))*y*y + (sin(K*x))*x ;
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
void set_initial_values_to_exact_solution()
{
// Loop over all elements and reset the values
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
 // Upcast from GeneralisedElement to the present element
 ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 unsigned nnode = el_pt->nnode();
 for(unsigned i=0;i<3;++i)
  {
   // Containers
   Vector<double> x(2),w(18,0.0);
   // Get node pointer
   Node* nod_pt = el_pt->node_pt(i);
   // Get x 
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get test
   TestSoln::get_exact_w(x,w);
   // Set value
   for(unsigned k=0;k<3;++k)
    {
    for(unsigned l=0 ; l<6;++l)
     {
      const unsigned index = el_pt->u_index_koiter_model(l,k);
      nod_pt->set_value(index,w[6*k+l]);
     }
    }
  }
  // Pin the internal dofs
  const unsigned n_internal_dofs = el_pt->number_of_internal_dofs();
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   // Containers
   Vector<double> x(2,0.0), s(2,0.0), w(18,0.0);
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
     internal_data_pt->set_value(index,w[6*k+0]);
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
     for(unsigned i=0;i<3;++i)
      {
       Node* nod_pt = el_pt->node_pt(i);
       // If it is on the bth boundary
       if(nod_pt -> is_on_boundary(b))
        {
         // Set the values of the unknowns that aren't pinned
         Vector<double> x(2,0.0),w(18,0.0);
         x[0] = nod_pt->x(0);
         x[1] = nod_pt->x(1);
         TestSoln::get_exact_w_radial(x,w);
         // Just set the nonzero values -> the others are pinned
         // Set value
         for(unsigned k=0;k<3;++k)
          {
          for(unsigned l=0 ; l<6;++l)
           {
            const unsigned index = el_pt->u_index_koiter_model(l,k);
            nod_pt->set_value(index,w[6*k+l]);
           }
          }
        }
      }
    }
   }
}// End loop over elements
}

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
 Outer_boundary1 = 1
};

double Element_area;

// The extra members required for flux type boundary conditions.
/// \short Number of "bulk" elements (We're attaching flux elements
/// to bulk mesh --> only the first Nkirchhoff_elements elements in
/// the mesh are bulk elements!)
// unsigned Nkirchhoff_elements;

/// \short Create bending moment elements on the b-th boundary of the
/// problems mesh 
void create_traction_elements(const unsigned &b, Mesh* const & bulk_mesh_py,
                            Mesh* const &surface_mesh_pt);

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);


/// \short Delete traction elements and wipe the surface mesh
void delete_traction_elements(Mesh* const &surface_mesh_pt);

/// \short Set pointer to prescribed-flux function for all elements
/// in the surface mesh
void set_prescribed_traction_pt();

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

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(outer_boundary_pt);

mesh_parameters.element_area() = element_area;

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
el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
el_pt->nu_pt() = &TestSoln::nu;
el_pt->d_pressure_dn_fct_pt() = &TestSoln::get_d_pressure_dn;
el_pt->d_pressure_dr_fct_pt() = &TestSoln::get_d_pressure_dr;
el_pt->d_pressure_d_grad_u_fct_pt() = &TestSoln::get_d_pressure_d_grad_u;
el_pt->thickness_pt() = &TestSoln::h;
// el_pt->eta_xy_pt() = &TestSoln::eta_xy;
// el_pt->eta_z_pt() = &TestSoln::eta_z;

//// For each node pin to the solution
//unsigned num_nod=el_pt->nnode();
//for (unsigned inod=0;inod<num_nod;inod++)
// {
// // Get nod
// Node* nod_pt=el_pt->node_pt(inod);
// Vector<double> x(2,0.0),w(18,0.0);
// x[0]=nod_pt->x(0);
// x[1]=nod_pt->x(1);
// TestSoln::get_exact_w(x,w);
//
// // Pin unknown values (everything except for the second normal derivative)
// for(unsigned i=0;i<3;++i)
//  {
//   for(unsigned l=0;l<6;++l)
//    {
//     unsigned index=el_pt->u_index_koiter_model(l,i);
//     nod_pt->set_value(index,w[index]);
//    }
// }
// }
}

// Loop over flux elements to pass pointer to prescribed traction function

/// Set pointer to prescribed traction function for traction elements
//set_prescribed_traction_pt();

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{

// Loop over all boundary nodes
//Just loop over outer boundary conditions
unsigned nbound = Outer_boundary1 + 1;
// Upcast to first element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));

for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get node
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 
 // Extract nodal coordinates from node:
 Vector<double> x(2),w(18);
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
} 

} // end set bc

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
create_traction_elements(const unsigned &b, Mesh* const &bulk_mesh_pt, 
                             Mesh* const &surface_mesh_pt)
{
}// end create traction elements

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
     {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary());}

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
       if(el_pt->node_pt(n)->is_on_boundary())
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
}// end create traction elements

//==start_of_set_prescribed_traction_pt===================================
/// Set pointer to prescribed traction function for all elements in the 
/// surface mesh
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::set_prescribed_traction_pt()
{
}// end of set prescribed flux pt

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
 
 Trace_file << TestSoln::p_mag << " " << "\n ";

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2 log(err/norm) \n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm <<" ";

// Only divide by norm if its nonzero
//some_file<<0.5*(log10(fabs(dummy_error))-log10(zero_norm))<<"\n";
some_file.close();

// Increment the doc_info number
Doc_info.number()++;

} // end of doc

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>
::delete_traction_elements(Mesh* const &surface_mesh_pt)
{
// How many surface elements are in the surface mesh
unsigned n_element = surface_mesh_pt->nelement();

// Loop over the surface elements
for(unsigned e=0;e<n_element;e++)
{
// Kill surface element
delete surface_mesh_pt->element_pt(e);
}

// Wipe the mesh
surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements


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

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Directory for solution
 string output_dir="RSLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 //CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);
 
 // P_step
 unsigned n_step=25;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);
 
 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,5,KoiterSteigmannPlateEquations> >problem(element_area);

 // Change some tolerances
 problem.max_residuals()=1e10;
 problem.max_newton_iterations()=30;
 // Initial guess `nearby'
 TestSoln::h=0.1;
 TestSoln::radius_of_curvature=Pi/2.0001;
 problem.set_initial_values_to_exact_solution();
 TestSoln::p_mag=pow(TestSoln::radius_of_curvature,3)*pow(TestSoln::h,2)/(12.*(1-pow(TestSoln::nu,2)));
 problem.newton_solve();

 // Test Parameters 
 TestSoln::h=0.1;
 TestSoln::radius_of_curvature=Pi/2.0;
 TestSoln::p_mag=pow(TestSoln::radius_of_curvature,3)*pow(TestSoln::h,2)/(12.*(1-pow(TestSoln::nu,2)));
 problem.newton_solve();
 problem.doc_solution();
 // Short circuit
 exit(0);

 // Set to zero before initial solve
 TestSoln::radius_of_curvature=0.0;
 // Loop Curvatures
 for(unsigned i=0;i<n_step;++i)
  {
  //Increment
  TestSoln::radius_of_curvature+=Pi/(n_step);
  oomph_info<<"Solving for curvature=" << TestSoln::radius_of_curvature<<"\n";
  try {
  problem.newton_solve();
  }
  catch(...)
   {
  // Dump the data 
  char filename[100];
  std::ofstream filestream;
  filestream.precision(15);
  sprintf(filename,"%s/fvk_circle_data%i-%f.dump",
          problem.Doc_info.directory().c_str(),
          problem.Doc_info.number(),
          TestSoln::p_mag
         );
  filestream.open(filename);
  problem.dump(filestream);
  filestream.close();
  exit(0); 
  }
  // Document
  problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << " kappa (" << TestSoln::radius_of_curvature << ")" << std::endl;
  oomph_info << "Current n_step  (" << n_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
}
 const unsigned n_h_step =10; 
 const double h_target = 0.1, h_initial = TestSoln::h;
 for(unsigned i=0;i<n_h_step;++i)
  {
  //Increment
  TestSoln::h+=(h_target - h_initial)/(n_h_step);
  oomph_info<<"Solving for curvature=" << TestSoln::radius_of_curvature<<"\n";
  try {
  problem.newton_solve();
  }
  catch(...)
   {
  // Dump the data 
  char filename[100];
  std::ofstream filestream;
  filestream.precision(15);
  sprintf(filename,"%s/fvk_circle_data%i-%f.dump",
          problem.Doc_info.directory().c_str(),
          problem.Doc_info.number(),
          TestSoln::p_mag
         );
  filestream.open(filename);
  problem.dump(filestream);
  filestream.close();
  exit(0); 
  }
  // Document
  problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Kappa (" << TestSoln::radius_of_curvature << ")" << std::endl;
  oomph_info << "Thickness (" << TestSoln::h << ")" << std::endl;
  oomph_info << "Current n_step  (" << n_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
}
 // Debug Jacobian
// problem.describe_dofs();
// oomph_info << "----------PRINTING GLOBAL JACOBIAN-----------" << std::endl;
// DoubleVector residuals;
// CRDoubleMatrix jac;
// DenseMatrix<double> jac_fd(problem.ndof(),problem.ndof(),0.0);
// problem.get_jacobian(residuals,jac);
// 
// // Open filestream, high precision
// std::ofstream filestream;
// filestream.precision(15);
// filestream.open("RESLT/Jacobian.csv");
// 
// // Output Jacobian in csv format
// for (unsigned n=0; n<problem.ndof(); ++n)
// {
//   for(unsigned m=0; m<problem.ndof();++m)
//     filestream<<jac(n,m)<<((m==problem.ndof()-1) ? "\n":",");
//     // if we are at the end of a row output "\n" else ","
// }
// filestream.close();
//
// oomph_info << "----------PRINTING FD JACOBIAN-----------" << std::endl;
// problem.get_fd_jacobian(residuals,jac_fd);
//// oomph::BlackBoxFDNewtonSolver::FD_step=1e-9;
// 
// // Open filestream, high precision
// filestream.precision(15);
// filestream.open("RESLT/FD_Jacobian.csv");
// // Output Jacobian in csv format
// for (unsigned n=0; n<problem.ndof(); ++n)
// {
//   for(unsigned m=0; m<problem.ndof();++m)
//     filestream<<jac_fd(n,m)<<((m==problem.ndof()-1) ? "\n":",");
//     // if we are at the end of a row output "\n" else ","
// }
// filestream.close();

} //End of main

