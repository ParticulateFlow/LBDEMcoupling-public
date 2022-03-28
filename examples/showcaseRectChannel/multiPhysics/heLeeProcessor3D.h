/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* This code was written with help of a Fortran code kindly provided
 * by Prof. Taehun Lee, and contains important contributions by
 * Andrea Parmigiani.
 */

#ifndef HE_LEE_PROCESSOR_3D_H
#define HE_LEE_PROCESSOR_3D_H

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/globalDefs.h"
#include "multiPhysics/interparticlePotential.h"

namespace plb {

/** GENERAL CONCEPTS
 *  ================
 *  In this implementation of the He/Lee binary fluid, the macroscopic
 *  variables like velocity and concentration are stored in separate
 *  fields (unlike the Shan/Chen implementation in Palabos, in which
 *  this kind of variables are stored in external scalars). The reason
 *  is that the model makes a massive use of finite difference schemes
 *  to compute derivatives, and derivatives of derivaties (up to fourth
 *  derivative). Derivatives are non-local, and so, before each deriva-
 *  tive some inter-block communication is needed to fill the envelopes.
 *  This is most easily achieved by putting each macroscopic variable in
 *  a separate field which then spontaneously updates its envelope when
 *  it has been modified.
 *
 *  LITERATURE
 *  ==========
 *  - T. Lee and C.-L. Lin, J Comp Phys 206 (2005), 16-47.
 *  - T. Lee, Comp Math App 58 (2009), 987-994.
 *
 *  VARIABLES
 *  =========
 *  The variables, each stored in a separate lattice or field, are:
 *  Block-Lattices:
 *  - g:   "The fluid" --> velocity u and pressure p1.
 *  - f:   "The concentration" --> Conc. C, density rho and
 *                                 potential mu.
 *  Scalar-Fields:
 *  - C:    Concentration of the heavy fluid.
 *  - rho:  Total density, summed over both components.
 *  - mu, laplaceMu:  Chemical potential and its Laplacian.
 *  - p1:             Flow pressure.
 *  Tensor-Fields:
 *  - gradC, gradMu:  Gradients of C and mu.
 *  - u:              Flow velocity.
 *
 *  Note that not all terms used in the He/Lee model are not represented
 *  as fields. For example, although the gradient of the pressure is used
 *  in the model, it is not stored in a tensor-field; instead, it is
 *  computed on-the-fly during the collision to save memory and to improve
 *  efficiency. This would not have been possible for the gradient of C,
 *  for example, because this quantity is used not only directly in the
 *  collision, but also for the computation of the velocity u.
 *
 *  DATA PROCESSORS FOR THE IMPLEMENTATION OF THE MODEL
 *  ===================================================
 *  The He/Lee algorithm is split over 4 data processors; the three first
 *  are needed to compute the macroscopic variables and their space deri-
 *  vatives from the populations f and g; the fourth implements the
 *  collision step.
 *
 *  - Compute_C_processor
 *      In: f, laplaceMu(t-1)
 *      Out: C
 *  - Compute_gradC_rho_mu_processor
 *      In: C
 *      Out: gradC, rho, mu
 *  - Compute_gradMu_laplaceMu_u_p1_processor
 *      In: g, C, rho, gradC, mu
 *      Out: gradMu, laplaceMu, u, p1
 *  - HeLeeCollisionProcessor
 *      In: f, g, C, rho, gradC, mu, gradMu, laplaceMu, u, p1
 *      Out: f, g
 *
 *  Note that all processors take the field f as their first argument,
 *  even if they don't need it, for a purely technical reason: in this
 *  way they can all be added as internal processors to f, to guarantee
 *  that they are executed automatically.
 *
 *  DATA PROCESSORS FOR SIMULATION SETUP
 *  ====================================
 *  To create an initial condition, the user must assign a value to C,
 *  u, and p1. Then, derivative variables can be computed by invoking
 *  Compute_gradC_rho_mu_processor. A special data processor is available
 *  to compute the derivatives of mu without overwriting u and p1:
 *  - Compute_gradMu_laplaceMu_processor
 *      In: mu
 *      Out: gradMu, laplaceMu
 *  Finally, the HeLeeCollisionProcessor can be used to initialize f and g
 *  to their equilibrium, by constructing this processor with an optional
 *  boolean argument initialize=true.
 *
 */

/// Compute the concentration C of the heavy phase.
/**
 *  In: f, laplaceMu(t-1)
 *  Out: C
 *
 *  Comment: laplaceMu is taken from the previous time step
 *  to preserve an explicit scheme.
 **/
template <typename T, template <typename U> class Descriptor>
class Compute_C_processor : public BoxProcessingFunctional3D {
public:
    Compute_C_processor(T M_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual Compute_C_processor<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T M;
};

/// Compute quantities derived from C.
/** Can be used both to implement the model and to set up the
 *  simulation.
 *  In: C
 *  Out: gradC, rho, mu
 */
template <typename T>
class Compute_gradC_rho_mu_processor : public BoxProcessingFunctional3D {
public:
    Compute_gradC_rho_mu_processor(T beta_, T kappa_, T rho_h_, T rho_l_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual Compute_gradC_rho_mu_processor<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T beta, kappa;
    T rho_h, rho_l;
};

/// Compute derivatives of the chemical potential mu, to be used during
/// the simulation setup.
/**
 *  In: mu
 *  Out: gradMu, laplaceMu
 **/
template <typename T, template <typename U> class Descriptor>
class Compute_gradMu_laplaceMu_processor : public BoxProcessingFunctional3D {
public:
    Compute_gradMu_laplaceMu_processor();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual Compute_gradMu_laplaceMu_processor<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Compute derivatives of the chemical potential mu, and the flow velocity
/// and pressure.
/**
 *  In: mu
 *  Out: gradMu, laplaceMu, u, p1.
 **/
template <typename T, template <typename U> class Descriptor>
class Compute_gradMu_laplaceMu_u_p1_processor : public BoxProcessingFunctional3D {
public:
    Compute_gradMu_laplaceMu_u_p1_processor(T rho_h_, T rho_l_, T RT_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T rho_h, rho_l, RT;
};

/// Implementation of the collision step. If the optional argument
/// initialize is set to true, the populations are initialized to
/// equilibrium.
/**
 *  In: f, g, C, rho, gradC, mu, gradMu, laplaceMu, u, p1
 *  Out: f, g
 **/
template <typename T, template <typename U> class Descriptor>
class HeLeeCollisionProcessor : public BoxProcessingFunctional3D {
public:
    HeLeeCollisionProcessor(
        T rho_h_, T rho_l_, T tau_h_, T tau_l_, T M_, T RT_, bool initialize_ = false);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual HeLeeCollisionProcessor<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    void computeAdvectionTerms(
        ScalarField3D<T> const &C, T &adv_gradC, T &bias_adv_gradC, plint iX, plint iY, plint iZ,
        plint iPop);

private:
    T rho_h, rho_l, tau_h, tau_l, M, RT;
    bool initialize;
};

}  // namespace plb

#endif  // HE_LEE_LATTICES_3D_H
