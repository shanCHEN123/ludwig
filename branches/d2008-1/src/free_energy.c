/*****************************************************************************
 *
 *  free_energy
 *
 *  This is the symmetric phi^4 free energy:
 *
 *  F[\phi] = (1/2) A \phi^2 + (1/4) B \phi^4 + (1/2) \kappa (\nabla\phi)^2
 *
 *  The first two terms represent the bulk free energy, while the
 *  final term penalises curvature in the interface. For a complete
 *  description see Kendon et al., J. Fluid Mech., 440, 147 (2001).
 *
 *  $Id: free_energy.c,v 1.4.2.2 2008-03-21 09:23:53 kevin Exp $
 *
 *  Kevin Stratford (kevin@epcc.ed.ac.uk)
 *  (c) 2007 The University of Edinburgh
 *
 *****************************************************************************/

#include <assert.h>
#include <math.h>

#include "pe.h"
#include "phi.h"
#include "runtime.h"
#include "utilities.h"
#include "free_energy.h"

static double A_     = -0.003125;
static double B_     = +0.003125;
static double C_     =  0.0;
static double kappa_ = +0.002;

/*****************************************************************************
 *
 *  free_energy_init
 *
 *  Get the user's parameters, if present. 
 *
 *****************************************************************************/

void init_free_energy() {

  int n;

  n = RUN_get_double_parameter("A", &A_);
  n = RUN_get_double_parameter("B", &B_);
  n = RUN_get_double_parameter("K", &kappa_);

#ifdef _SINGLE_FLUID_
#else
  if (A_ > 0.0) {
    fatal("The free energy parameter A must be negative\n");
  }
  if (B_ < 0.0) {
    fatal("The free energy parameter B must be positive\n");
  }
  if (kappa_ < 0.0) {
    fatal("The free energy parameter kappa must be positive\n");
  }

  info("\nSymmetric phi^4 free energy:\n");
  info("Bulk parameter A      = %f\n", A_);
  info("Bulk parameter B      = %f\n", B_);
  info("Surface penalty kappa = %f\n", kappa_);
  info("Surface tension       = %f\n", surface_tension());
  info("Interfacial width     = %f\n", interfacial_width());
#endif

  return;
}

/*****************************************************************************
 *
 *  free_energy_A
 *  free_energy_B
 *  free_energy_K
 *
 *****************************************************************************/

double free_energy_A() {
  return A_;
}
double free_energy_B() {
  return B_;
}
double free_energy_K() {
  return kappa_;
}

void free_energy_set_A(double a) {
  assert(a < 0.0);
  A_ = a;
  return;
}

void free_energy_set_B(double b) {
  assert(b > 0.0);
  B_ = b;
  return;
}

void free_energy_set_kappa(double k) {
  assert(k > 0.0);
  kappa_ = k;
  return;
}

/*****************************************************************************
 *
 *  surface_tension
 *
 *  Return the theoretical surface tension for the model.
 *
 *****************************************************************************/

double surface_tension() {

  return sqrt(-8.0*kappa_*A_*A_*A_ / (9.0*B_*B_));
}

/*****************************************************************************
 *
 *  interfacial_width
 *
 *  Return the theoretical interfacial width. Note that there is a
 *  typo in 6.17 of Kendon et al (2001); the factor of 2 should be in
 *  numerator.
 *
 *****************************************************************************/

double interfacial_width() {

  return sqrt(-2.0*kappa_ / A_);
}

/*****************************************************************************
 *
 *  free_energy_get_chemical_potential
 *
 *  Return the chemical potential at given position index.
 *
 *****************************************************************************/

double free_energy_get_chemical_potential(const int index) {

  double phi, delsq_phi, delsq_sq_phi, mu;

  phi = phi_get_phi_site(index);
  delsq_phi = phi_get_delsq_phi_site(index);
  delsq_sq_phi = phi_get_delsq_sq_phi_site(index);

  mu = phi*(A_ + B_*phi*phi) - kappa_*delsq_phi + C_*delsq_sq_phi;

  return mu;
}

/*****************************************************************************
 *
 *  free_energy_get_chemical_stress
 *
 *  Return the chemical stress tensor for given position index.
 *  P_ab = [1/2 A phi^2 + 3/4 B phi^4 - kappa phi \nabla^2 phi
 *       -  1/2 kappa (\nbla phi)^2] \delta_ab
 *       +  kappa \nalba_a phi \nabla_b phi
 *
 *****************************************************************************/

void free_energy_get_chemical_stress(const int index, double p[3][3]) {

  int ia, ib;
  double phi, bulk, delsq_phi, grad_phi_sq;
  double grad_phi[3];
  extern const double d_[3][3]; /* Pending Refactor util etc. */ 

  phi = phi_get_phi_site(index);
  phi_get_grad_phi_site(index, grad_phi);
  delsq_phi = phi_get_delsq_phi_site(index);

  bulk = 0.5*phi*phi*(A_ + 1.5*B_*phi*phi);
  grad_phi_sq = dot_product(grad_phi, grad_phi);

  for (ia = 0; ia < 3; ia++) {
    for (ib = 0; ib < 3; ib++) {
      p[ia][ib] = (bulk - kappa_*(phi*delsq_phi + 0.5*grad_phi_sq))*d_[ia][ib]
	+ kappa_*grad_phi[ia]*grad_phi[ib];
    }
  }

  return;
}

/*****************************************************************************
 *
 *  free_energy_get_isotropic_pressure
 *
 *  Return the isotrpoic part of the pressure tensor.
 *  P_0 = [1/2 A phi^2 + 3/4 B phi^4 - kappa phi \nabla^2 phi
 *       -  1/2 kappa (\nabla phi)^2]
 *
 *****************************************************************************/

double free_energy_get_isotropic_pressure(const int index) {

  double p0, phi, bulk, delsq_phi, grad_phi_sq;
  double grad_phi[3];

  phi = phi_get_phi_site(index);
  phi_get_grad_phi_site(index, grad_phi);
  delsq_phi = phi_get_delsq_phi_site(index);

  bulk = 0.5*phi*phi*(A_ + 1.5*B_*phi*phi);
  grad_phi_sq = dot_product(grad_phi, grad_phi);
  p0 = bulk - kappa_*(phi*delsq_phi + 0.5*grad_phi_sq);

  return p0;
}

/*****************************************************************************
 *
 *  free_energy_density
 *
 *  Return the free energy density
 *  E = (1/2) A phi^2 + (1/4) B phi^4 + (1/2) kappa (\nabla phi)^2
 *    + (1/2) C (\nabla^2 phi)^2
 *
 *****************************************************************************/

double free_energy_density(const int index) {

  double e, bulk;
  double phi, dphi[3], delsq;

  phi = phi_get_phi_site(index);
  phi_get_grad_phi_site(index, dphi);
  delsq = phi_get_delsq_phi_site(index);

  bulk = phi*phi*(A_ + 0.5*B_*phi*phi);
  e = 0.5*(bulk + kappa_*dot_product(dphi, dphi) + C_*delsq*delsq);

  return e;
}
