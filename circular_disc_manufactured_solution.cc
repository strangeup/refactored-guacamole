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
#include <boost/math/special_functions/ellint_2.hpp>

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
double h = 0.2;
double eta_u = 1;
double eta_sigma = 1;
double eta_bending = h;

/// Helper function to update parameters that could potentially depend on nu and
// h
void update_problem_parameters()
 {
  eta_u = 1;
  eta_sigma =1;
  eta_bending = h;
 }

/// Different nondimensionalisations
/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Use a Mathematica friendly alias for std::pow
template<typename T1, typename T2>
double Power(const T1& arg, const T2& exp)
{ return std::pow(arg,exp); }

// Assigns the value of pressure depending on the position (x,y)
inline void get_koiter_steigmann_pressure(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& d_ui_dx_j,  const Vector<double>& ni, 
 Vector<double>& pressure)
{
 const double x=xi[0];
 // Create aliases
 double (*Sin)(double theta);
 double (*Cos)(double theta);
 double (*Sqrt)(double theta);
 Sin = & std::sin;
 Cos = & std::cos;
 Sqrt = & std::sqrt;

 // Full forcing for cheeky ellipse solution
 pressure[0] = -(Sqrt(2)*Power(h,4)*Cos(x)*Sin(x)*(4 - 3*Power(h,2)*
  Power(Cos(x),2) - 5*Power(h,2)*Power(Sin(x),2) + Power(h,4)*Power(Sin(x),4)))/
  (3.*(-1 + Power(nu,2))*Power(2 - Power(h,2) + Power(h,2)*Cos(2*x),2.5));

 pressure[1] =0.0;

 pressure[2] =-(Power(h,3)*Cos(x))/(12.*(-1 + Power(nu,2)));
}

inline void get_foppl_correction_pressure(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& d_ui_dx_j,  const Vector<double>& ni, 
 Vector<double>& pressure)
{
 // Create aliases
 double (*Sin)(double theta);
 double (*Cos)(double theta);
 double (*Sqrt)(double theta);
 Sin = & std::sin;
 Cos = & std::cos;
 Sqrt = & std::sqrt;
 // Basic ordinates
 const double x=xi[0], Chi = 1.0 - Power(h*Sin(x),2);

 // Pressure for cheeky bend to ellipse case
 pressure[0] = -(Power(h,4)*Cos(x)*(Chi*(3 + Chi*(4 - Sqrt(Chi) - 7*Chi + 
  5*Power(Chi,1.5))) + (-3 + Chi*(-1 - 2*Chi + 3*Power(Chi,1.5)))*Power(h,2)*
  Power(Cos(x),2))*Sin(x))/(12.*Power(Chi,2.5)*(-1 + Power(nu,2)));

 pressure[1] = 0.0;

 pressure[2] =-(Power(h,3)*Cos(x)*(Chi*(6 + Chi*(-4 - 11*Sqrt(Chi) + 18*Chi 
  + 2*Power(Chi,1.5) - 20*Power(Chi,2) + 13*Power(Chi,2.5))) + 2*(-1 + Sqrt(Chi))
  *(13 - Sqrt(Chi) - 13*Chi + Power(Chi,1.5) - 2*Power(Chi,2) + 6*Power(Chi,3))
  *Power(h,2)*Power(Cos(x),2) + (5 - 7*Sqrt(Chi))*Power(h,4)*Power(Sin(2*x),2)))
  /(48.*Power(Chi,2.5)*(-1 + Power(nu,2)));
 }

 void pressure(const Vector<double>& xi, Vector<double>& pressure)
  {
   get_foppl_correction_pressure(xi,Vector<double>(3),DenseMatrix<double>(2,2),Vector<double>(3),pressure);
  }

// Alias for pressurre function 
typedef void (*PressureFctPt)(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj, const Vector<double>& ni, 
 Vector<double>& pressure);

// Assigns the value of pressure depending on the position (x,y)
inline void get_empty_d_pressure_dn(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& dui_dxj, const Vector<double>& ni, 
 DenseMatrix<double>& d_pressure_dn)
{/* Empty */}

// Assigns the value of pressure depending on the position (x,y)
void get_empty_d_pressure_dr(const Vector<double>& xi,const Vector<double>& ui,
  const DenseMatrix<double>& dui_dxj,
 const Vector<double>& ni, DenseMatrix<double>& d_pressure_dn)
{/* Empty */}

// Assigns the value of pressure depending on the position (x,y)
inline void get_empty_d_pressure_d_grad_u(const Vector<double>& xi,const Vector<double>& ui,
 const DenseMatrix<double>& d_ui_dx_j,
 const Vector<double>& ni, RankThreeTensor<double>& d_pressure_du_grad)
{/* Empty */}

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
 // Basic ordinates 
 const double x = xi[0]/*, y = xi[1] */;

 // Create aliases
 double (*EllipticE)(double m, double theta);
 double (*Power)(double base, int exponent);
 double (*Sin)(double theta);
 double (*Cos)(double theta);
 double (*Sqrt)(double theta);
 EllipticE = boost::math::ellint_2;
 Power = & std::pow;
 Sin = & std::sin;
 Cos = & std::cos;
 Sqrt = & std::sqrt;
 // Now fill in
 const unsigned ndisp = 3, ndoftype = 6;
 w=Vector<double>(ndisp*ndoftype,0.0);
 // In--plane displacement
 w[0] = -x + EllipticE(h,x);
 w[1] = -1 + Sqrt(1 - Power(h,2)*Power(Sin(x),2));
 w[3] = -((Power(h,2)*Cos(x)*Sin(x))/Sqrt(1 - Power(h,2)*Power(Sin(x),2)));
 // Out-of-plane displacement
 w[2*ndoftype + 0] = h*(-1 + Cos(x));
 w[2*ndoftype + 1] = -(h*Sin(x));
 w[2*ndoftype + 3] = -(h*Cos(x));;
}

// Bool for high resolution
bool High_resolution = false;

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
class ProblemWithDocInfo : public virtual Problem
{
public:
 ProblemWithDocInfo(){}
 ~ProblemWithDocInfo(){}
 /// Doc info object for labeling output
 DocInfo Doc_info;
 /// Set initial values to exact solution
 virtual void set_initial_values_to_exact_solution()=0;
};

//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual ProblemWithDocInfo
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
 // Doc the solution
 doc_solution();
}

/// Update the problem specs before solve: Re-apply boundary conditions
/// Empty as the boundary conditions stay fixed
void actions_before_newton_solve()
{
// Update any global problem parameters
TestSoln::update_problem_parameters();
// Update any problem paramters in elemnts
complete_problem_setup();
}

/// Set the initial values to the exact solution (useful for debugging)
void set_initial_values_to_exact_solution();

/// Enable finite element Jacobian
void enable_finite_difference_jacobian()
 {
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to the present element
  ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
  el_pt->enable_finite_difference_jacobian();
 }
}
/// Disable finite element Jacobian
void disable_finite_difference_jacobian()
 {
 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to the present element
  ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
  el_pt->disable_finite_difference_jacobian();
 }
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

void print_jacobian();

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

/// Pointer to Pressure function
TestSoln::PressureFctPt Pressure_fct_pt;

/// Get pressure fct pt
public:
TestSoln::PressureFctPt& pressure_fct_pt() {return Pressure_fct_pt;}
const TestSoln::PressureFctPt pressure_fct_pt() const  {return Pressure_fct_pt;}
}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area), Pressure_fct_pt(0)
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
//rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

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
   TestSoln::get_exact_w(x,w);
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
         TestSoln::get_exact_w/*_radial*/(x,w);
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

 TestSoln::get_exact_w(x,w);

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=0;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal
    // if(l != 3)  
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
el_pt->pressure_fct_pt() = Pressure_fct_pt;
el_pt->nu_pt() = &TestSoln::nu;
el_pt->d_pressure_dn_fct_pt() = &TestSoln::get_empty_d_pressure_dn;
el_pt->d_pressure_dr_fct_pt() = &TestSoln::get_empty_d_pressure_dr;
el_pt->d_pressure_d_grad_u_fct_pt() = &TestSoln::get_empty_d_pressure_d_grad_u;
el_pt->thickness_pt() = &TestSoln::eta_bending;
el_pt->eta_sigma_pt() = &TestSoln::eta_sigma;
//// Test case for thickness scaling (DISABLED)
// el_pt->eta_u_pt() = &TestSoln::eta_u;
// el_pt->enable_finite_difference_jacobian();
}

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

 TestSoln::get_exact_w(x,w);

 // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i=0;i<3;++i)
  {
   for(unsigned l=0;l<6;++l)
    {
    // Don't pin the second normal
    //if(l != 3)  
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
if (TestSoln::High_resolution)
 {npts = 50;}
else
 {npts = 6;}

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

//  Output exact solution
sprintf(filename,"%s/pressure_soln%i-%f.dat","RESLT",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::pressure); 
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
 r[0]=0.5;
 Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);
 // Call the interpolated_u function
 Vector<Vector<double> > u_i(3,Vector<double>(6,0.0));
 dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_koiter_plate(s,u_i);

 oomph_info << "r at (0.5,0.0): "  
            <<"("<< 0.5+u_i[0][0] <<","<<0.0+u_i[1][0]<<","<<u_i[2][0]<<")"<< std::endl;
 
 Trace_file << TestSoln::h <<" "<< u_i[0][0] <<" "<<u_i[1][0]<<" "<<u_i[2][0]<< std::endl;

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

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::print_jacobian()
{
 // Print the Jacobian
 describe_dofs();
 oomph_info << "----------PRINTING JACOBIAN-----------" << std::endl;
 DoubleVector residuals;
 CRDoubleMatrix jac;
 get_jacobian(residuals,jac);
 
 // Open filestream, high precision
 std::ofstream filestream;
 filestream.precision(15);
 // Filename
 char jacobian_filename[100];
 sprintf(jacobian_filename, "%s/Jacobian.csv",Doc_info.directory().c_str());
 filestream.open(jacobian_filename);
 
 // Output Jacobian in csv format
 for (unsigned n=0; n<ndof(); ++n)
 {
   for(unsigned m=0; m<ndof();++m)
     filestream<<jac(n,m)<<((m==ndof()-1) ? "\n":",");
     // if we are at the end of a row output "\n" else ","
 }
 filestream.close();
 exit(0);
}

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

 // High resolution solutions
 CommandLineArgs::specify_command_line_flag("--high_resolution");

 // High resolution solutions
 CommandLineArgs::specify_command_line_flag("--do_fvk_correction");

 // Directory for solution
 string output_dir="RSLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--use_fd_jacobian");
 
 // P_step
 unsigned n_step=24;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);
 
 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Update any problem parameters
 TestSoln::update_problem_parameters();

 // If validate flag provided
 const bool validate=CommandLineArgs::
   command_line_flag_has_been_set("--validate");

 // Set the flag 
 TestSoln::High_resolution=CommandLineArgs::
   command_line_flag_has_been_set("--high_resolution");

 // Set the flag 
 const bool do_fvk_correction=CommandLineArgs::
   command_line_flag_has_been_set("--do_fvk_correction");

 // Set the flag 
 const bool use_fd_jacobian=CommandLineArgs::
   command_line_flag_has_been_set("--use_fd_jacobian");

 // Validate block
 if(validate)
  {
  // Test case parameter
  // Problem instance
  UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,5,KoiterSteigmannPlateEquations> >problem_1(element_area);
  problem_1.set_initial_values_to_exact_solution();
  problem_1.pressure_fct_pt() = &TestSoln::get_koiter_steigmann_pressure;
  // Problem instance
  UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement<2,2,5,FoepplVonKarmanCorrectionEquations> >problem_2(element_area);
  problem_2.set_initial_values_to_exact_solution();
  problem_2.pressure_fct_pt() = &TestSoln::get_foppl_correction_pressure;

  // The problem we are interested in solving based on flag
  Problem* problem_pt;
  // Use Koiter Steigmann
  if(!do_fvk_correction)
   {problem_pt = &problem_1; }
  // Use FvK correction
  else
   {problem_pt = &problem_2; }

  // Change some tolerances
  problem_pt->max_residuals()=1e10;
  problem_pt->max_newton_iterations()=30;
  // Newton Solve
  problem_pt->newton_solve();
  exit(0);
  }

 // Problem instance
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement
   <2,2,5,KoiterSteigmannPlateEquations> >problem_1(element_area);
 UnstructuredFvKProblem<LargeDisplacementPlateC1CurvedBellElement
   <2,2,5,FoepplVonKarmanCorrectionEquations> >problem_2(element_area);

 // Use FD jacobian
 if(use_fd_jacobian)
  {
   problem_1.enable_finite_difference_jacobian();
   problem_2.enable_finite_difference_jacobian();
  }

 // The problem we are interested in solving based on flag
 ProblemWithDocInfo* problem_pt;
 // Use Koiter Steigmann
 if(!do_fvk_correction)
  {problem_pt = &problem_1; }
 // Use FvK correction
 else
  {problem_pt = &problem_2; }


 // Change some tolerances
 problem_pt->max_residuals()=1e10;
 problem_pt->max_newton_iterations()=30;
 problem_pt->set_initial_values_to_exact_solution();

 // h max
 const double h_max = 0.25, h_start = 0.01, h_step = (h_max-h_start)/n_step;
 TestSoln::h =  h_start;
 TestSoln::update_problem_parameters();

 // Loop Curvatures
 for(unsigned i=0;i<n_step+1;++i)
  {
  oomph_info<<"Thickness =" << TestSoln::h<<"\n";
  try {
   problem_pt->newton_solve();
  }
  catch(...)
   {
   // Dump the data 
   char filename[100];
   std::ofstream filestream;
   filestream.precision(15);
   sprintf(filename,"%s/circle_data%i-%f.dump",
           problem_pt->Doc_info.directory().c_str(),
           problem_pt->Doc_info.number(),
           TestSoln::h
          );
   filestream.open(filename);
   problem_pt->dump(filestream);
   filestream.close();
   exit(0); 
  }
  // Document
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Thickness (" << TestSoln::h << ")" << std::endl;
  oomph_info << "Current n_step  (" << n_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem_pt->Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
  //Increment
  TestSoln::h+=h_step;
  // Update any problem parameters
  TestSoln::update_problem_parameters();
}
} //End of main

