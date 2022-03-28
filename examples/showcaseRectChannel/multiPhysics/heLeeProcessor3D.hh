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
 * by Prof. Taehun Lee, and is massively co-authored by
 * Andrea Parmigiani.
 */

#ifndef HE_LEE_PROCESSOR_3D_HH
#define HE_LEE_PROCESSOR_3D_HH

#include "heLeeProcessor3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T>
class FiniteDifference {
public:
    FiniteDifference(ScalarField3D<T> const &field_, plint x_, plint y_, plint z_) :
        field(field_), x(x_), y(y_), z(z_)
    { }
    T val(plint iX, plint iY, plint iZ) const
    {
        return field.get(iX, iY, iZ);
    }
    void centralGradient(Array<T, 3> &gradient) const
    {
        const T c1 = (T)2. / (T)9.;
        const T c2 = (T)1. / (T)18.;
        const T c3 = (T)1. / (T)72.;
        gradient[0] =
            (val(x + 1, y, z) - val(x - 1, y, z)) * c1
            + (val(x + 1, y + 1, z) - val(x - 1, y - 1, z) + val(x + 1, y - 1, z)
               - val(x - 1, y + 1, z) + val(x + 1, y, z + 1) - val(x - 1, y, z - 1)
               + val(x + 1, y, z - 1) - val(x - 1, y, z + 1))
                  * c2
            + (val(x + 1, y + 1, z + 1) - val(x - 1, y - 1, z - 1) + val(x + 1, y - 1, z + 1)
               - val(x - 1, y + 1, z - 1) + val(x + 1, y + 1, z - 1) - val(x - 1, y - 1, z + 1)
               + val(x + 1, y - 1, z - 1) - val(x - 1, y + 1, z + 1))
                  * c3;
        gradient[1] =
            (val(x, y + 1, z) - val(x, y - 1, z)) * c1
            + (val(x + 1, y + 1, z) - val(x - 1, y - 1, z) + val(x - 1, y + 1, z)
               - val(x + 1, y - 1, z) + val(x, y + 1, z + 1) - val(x, y - 1, z - 1)
               + val(x, y + 1, z - 1) - val(x, y - 1, z + 1))
                  * c2
            + (val(x + 1, y + 1, z + 1) - val(x - 1, y - 1, z - 1) + val(x + 1, y + 1, z - 1)
               - val(x - 1, y - 1, z + 1) + val(x - 1, y + 1, z + 1) - val(x + 1, y - 1, z - 1)
               + val(x - 1, y + 1, z - 1) - val(x + 1, y - 1, z + 1))
                  * c3;
        gradient[2] =
            (val(x, y, z + 1) - val(x, y, z - 1)) * c1
            + (val(x + 1, y, z + 1) - val(x - 1, y, z - 1) + val(x - 1, y, z + 1)
               - val(x + 1, y, z - 1) + val(x, y + 1, z + 1) - val(x, y - 1, z - 1)
               + val(x, y - 1, z + 1) - val(x, y + 1, z - 1))
                  * c2
            + (val(x + 1, y + 1, z + 1) - val(x - 1, y - 1, z - 1) + val(x + 1, y - 1, z + 1)
               - val(x - 1, y + 1, z - 1) + val(x - 1, y + 1, z + 1) - val(x + 1, y - 1, z - 1)
               + val(x - 1, y - 1, z + 1) - val(x + 1, y + 1, z - 1))
                  * c3;
    }
    void biasedGradient(Array<T, 3> &gradient) const
    {
        const T c1 = (T)1. / (T)9.;
        const T c2 = (T)1. / (T)36.;
        const T c3 = (T)1. / (T)144.;
        gradient[0] =
            (-val(x + 2, y, z) + 4. * val(x + 1, y, z) - 4. * val(x - 1, y, z) + val(x - 2, y, z))
                * c1
            + (-val(x + 2, y + 2, z) + 4. * val(x + 1, y + 1, z) - 4. * val(x - 1, y - 1, z)
               + val(x - 2, y - 2, z) - val(x + 2, y - 2, z) + 4. * val(x + 1, y - 1, z)
               - 4. * val(x - 1, y + 1, z) + val(x - 2, y + 2, z) - val(x + 2, y, z + 2)
               + 4. * val(x + 1, y, z + 1) - 4. * val(x - 1, y, z - 1) + val(x - 2, y, z - 2)
               - val(x + 2, y, z - 2) + 4. * val(x + 1, y, z - 1) - 4. * val(x - 1, y, z + 1)
               + val(x - 2, y, z + 2))
                  * c2
            + (-val(x + 2, y + 2, z + 2) + 4. * val(x + 1, y + 1, z + 1)
               - 4. * val(x - 1, y - 1, z - 1) + val(x - 2, y - 2, z - 2) - val(x + 2, y - 2, z + 2)
               + 4. * val(x + 1, y - 1, z + 1) - 4. * val(x - 1, y + 1, z - 1)
               + val(x - 2, y + 2, z - 2) - val(x + 2, y + 2, z - 2) + 4. * val(x + 1, y + 1, z - 1)
               - 4. * val(x - 1, y - 1, z + 1) + val(x - 2, y - 2, z + 2) - val(x + 2, y - 2, z - 2)
               + 4. * val(x + 1, y - 1, z - 1) - 4. * val(x - 1, y + 1, z + 1)
               + val(x - 2, y + 2, z + 2))
                  * c3;
        gradient[1] =
            (-val(x, y + 2, z) + 4. * val(x, y + 1, z) - 4. * val(x, y - 1, z) + val(x, y - 2, z))
                * c1
            + (-val(x + 2, y + 2, z) + 4. * val(x + 1, y + 1, z) - 4. * val(x - 1, y - 1, z)
               + val(x - 2, y - 2, z) - val(x - 2, y + 2, z) + 4. * val(x - 1, y + 1, z)
               - 4. * val(x + 1, y - 1, z) + val(x + 2, y - 2, z) - val(x, y + 2, z + 2)
               + 4. * val(x, y + 1, z + 1) - 4. * val(x, y - 1, z - 1) + val(x, y - 2, z - 2)
               - val(x, y + 2, z - 2) + 4. * val(x, y + 1, z - 1) - 4. * val(x, y - 1, z + 1)
               + val(x, y - 2, z + 2))
                  * c2
            + (-val(x + 2, y + 2, z + 2) + 4. * val(x + 1, y + 1, z + 1)
               - 4. * val(x - 1, y - 1, z - 1) + val(x - 2, y - 2, z - 2) - val(x + 2, y + 2, z - 2)
               + 4. * val(x + 1, y + 1, z - 1) - 4. * val(x - 1, y - 1, z + 1)
               + val(x - 2, y - 2, z + 2) - val(x - 2, y + 2, z + 2) + 4. * val(x - 1, y + 1, z + 1)
               - 4. * val(x + 1, y - 1, z - 1) + val(x + 2, y - 2, z - 2) - val(x - 2, y + 2, z - 2)
               + 4. * val(x - 1, y + 1, z - 1) - 4. * val(x + 1, y - 1, z + 1)
               + val(x + 2, y - 2, z + 2))
                  * c3;
        gradient[2] =
            (-val(x, y, z + 2) + 4. * val(x, y, z + 1) - 4. * val(x, y, z - 1) + val(x, y, z - 2))
                * c1
            + (-val(x + 2, y, z + 2) + 4. * val(x + 1, y, z + 1) - 4. * val(x - 1, y, z - 1)
               + val(x - 2, y, z - 2) - val(x - 2, y, z + 2) + 4. * val(x - 1, y, z + 1)
               - 4. * val(x + 1, y, z - 1) + val(x + 2, y, z - 2) - val(x, y + 2, z + 2)
               + 4. * val(x, y + 1, z + 1) - 4. * val(x, y - 1, z - 1) + val(x, y - 2, z - 2)
               - val(x, y - 2, z + 2) + 4. * val(x, y - 1, z + 1) - 4. * val(x, y + 1, z - 1)
               + val(x, y + 2, z - 2))
                  * c2
            + (-val(x + 2, y + 2, z + 2) + 4. * val(x + 1, y + 1, z + 1)
               - 4. * val(x - 1, y - 1, z - 1) + val(x - 2, y - 2, z - 2) - val(x + 2, y - 2, z + 2)
               + 4. * val(x + 1, y - 1, z + 1) - 4. * val(x - 1, y + 1, z - 1)
               + val(x - 2, y + 2, z - 2) - val(x - 2, y + 2, z + 2) + 4. * val(x - 1, y + 1, z + 1)
               - 4. * val(x + 1, y - 1, z - 1) + val(x + 2, y - 2, z - 2) - val(x - 2, y - 2, z + 2)
               + 4. * val(x - 1, y - 1, z + 1) - 4. * val(x + 1, y + 1, z - 1)
               + val(x + 2, y + 2, z - 2))
                  * c3;
    }

    T laplacian() const
    {
        static const T c1 = (T)4. / (T)9.;
        static const T c2 = (T)1. / (T)9.;
        static const T c3 = (T)1. / (T)36.;
        static const T c4 = (T)38 / (T)9.;

        return (val(x + 1, y, z) + val(x - 1, y, z) + val(x, y + 1, z) + val(x, y - 1, z)
                + val(x, y, z + 1) + val(x, y, z - 1))
                   * c1
               + (val(x + 1, y + 1, z) + val(x - 1, y - 1, z) + val(x + 1, y - 1, z)
                  + val(x - 1, y + 1, z) + val(x + 1, y, z + 1) + val(x - 1, y, z - 1)
                  + val(x + 1, y, z - 1) + val(x - 1, y, z + 1) + val(x, y + 1, z + 1)
                  + val(x, y - 1, z - 1) + val(x, y + 1, z - 1) + val(x, y - 1, z + 1))
                     * c2
               + (val(x + 1, y + 1, z + 1) + val(x + 1, y + 1, z - 1) + val(x + 1, y - 1, z + 1)
                  + val(x + 1, y - 1, z - 1) + val(x - 1, y + 1, z + 1) + val(x - 1, y + 1, z - 1)
                  + val(x - 1, y - 1, z + 1) + val(x - 1, y - 1, z - 1))
                     * c3
               - val(x, y, z) * c4;
    }

private:
    ScalarField3D<T> const &field;
    plint x, y, z;
};

/* *************** Compute_C_processor ***************** */

template <typename T, template <typename U> class Descriptor>
Compute_C_processor<T, Descriptor>::Compute_C_processor(T M_) : M(M_)
{ }

template <typename T, template <typename U> class Descriptor>
void Compute_C_processor<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    BlockLattice3D<T, Descriptor> &f = *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    // Laplacian of chemical potential at previous time step t-1.
    ScalarField3D<T> &laplaceMu = *dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    ScalarField3D<T> &C = *dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Use templates to compute order-0 moment of f.
                C.get(iX, iY, iZ) = momentTemplates<T, Descriptor>::get_rhoBar(f.get(iX, iY, iZ))
                                    + 0.5 * M * laplaceMu.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Compute_C_processor<T, Descriptor> *Compute_C_processor<T, Descriptor>::clone() const
{
    return new Compute_C_processor<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void Compute_C_processor<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // f
    modified[1] = modif::nothing;          // laplaceMu(t-1)
    modified[2] = modif::staticVariables;  // C
}

/* *************** Compute_gradC_rho_mu_processor ***************** */

template <typename T>
Compute_gradC_rho_mu_processor<T>::Compute_gradC_rho_mu_processor(
    T beta_, T kappa_, T rho_h_, T rho_l_) :
    beta(beta_), kappa(kappa_), rho_h(rho_h_), rho_l(rho_l_)
{ }

template <typename T>
void Compute_gradC_rho_mu_processor<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    // blocks[0] stands for the unused f.
    ScalarField3D<T> &C = *dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    TensorField3D<T, 3> &gradC = *dynamic_cast<TensorField3D<T, 3> *>(blocks[2]);
    ScalarField3D<T> &rho = *dynamic_cast<ScalarField3D<T> *>(blocks[3]);
    ScalarField3D<T> &mu = *dynamic_cast<ScalarField3D<T> *>(blocks[4]);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                FiniteDifference<T>(C, iX, iY, iZ).centralGradient(gradC.get(iX, iY, iZ));
                T C_ = C.get(iX, iY, iZ);
                rho.get(iX, iY, iZ) = C_ * (rho_h - rho_l) + rho_l;
                mu.get(iX, iY, iZ) = 4. * beta * C_ * (C_ - 0.5) * (C_ - 1.)
                                     - kappa * FiniteDifference<T>(C, iX, iY, iZ).laplacian();
            }
        }
    }
}

template <typename T>
Compute_gradC_rho_mu_processor<T> *Compute_gradC_rho_mu_processor<T>::clone() const
{
    return new Compute_gradC_rho_mu_processor<T>(*this);
}

template <typename T>
void Compute_gradC_rho_mu_processor<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // f
    modified[1] = modif::nothing;          // C
    modified[2] = modif::staticVariables;  // gradC
    modified[3] = modif::staticVariables;  // rho
    modified[4] = modif::staticVariables;  // mu
}

/* *************** Compute_gradMu_laplaceMu_processor ***************** */

template <typename T, template <typename U> class Descriptor>
Compute_gradMu_laplaceMu_processor<T, Descriptor>::Compute_gradMu_laplaceMu_processor()
{ }

template <typename T, template <typename U> class Descriptor>
void Compute_gradMu_laplaceMu_processor<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    // blocks[0] stands for the unused f.
    ScalarField3D<T> &mu = *dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    TensorField3D<T, 3> &gradMu = *dynamic_cast<TensorField3D<T, 3> *>(blocks[2]);
    ScalarField3D<T> &laplaceMu = *dynamic_cast<ScalarField3D<T> *>(blocks[3]);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                FiniteDifference<T>(mu, iX, iY, iZ).centralGradient(gradMu.get(iX, iY, iZ));
                laplaceMu.get(iX, iY, iZ) = FiniteDifference<T>(mu, iX, iY, iZ).laplacian();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Compute_gradMu_laplaceMu_processor<T, Descriptor>
    *Compute_gradMu_laplaceMu_processor<T, Descriptor>::clone() const
{
    return new Compute_gradMu_laplaceMu_processor<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void Compute_gradMu_laplaceMu_processor<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // f
    modified[1] = modif::nothing;          // mu
    modified[2] = modif::staticVariables;  // gradMu
    modified[3] = modif::staticVariables;  // laplaceMu
}

/* *************** Compute_gradMu_laplaceMu_u_p1_processor ***************** */

template <typename T, template <typename U> class Descriptor>
Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>::Compute_gradMu_laplaceMu_u_p1_processor(
    T rho_h_, T rho_l_, T RT_) :
    rho_h(rho_h_), rho_l(rho_l_), RT(RT_)
{ }

template <typename T, template <typename U> class Descriptor>
void Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    // blocks[0] stands for the unused f.
    BlockLattice3D<T, Descriptor> &g = *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[1]);
    ScalarField3D<T> &C = *dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    ScalarField3D<T> &rho = *dynamic_cast<ScalarField3D<T> *>(blocks[3]);
    TensorField3D<T, 3> &gradC = *dynamic_cast<TensorField3D<T, 3> *>(blocks[4]);
    ScalarField3D<T> &mu = *dynamic_cast<ScalarField3D<T> *>(blocks[5]);
    TensorField3D<T, 3> &gradMu = *dynamic_cast<TensorField3D<T, 3> *>(blocks[6]);
    ScalarField3D<T> &laplaceMu = *dynamic_cast<ScalarField3D<T> *>(blocks[7]);
    TensorField3D<T, 3> &u = *dynamic_cast<TensorField3D<T, 3> *>(blocks[8]);
    ScalarField3D<T> &p1 = *dynamic_cast<ScalarField3D<T> *>(blocks[9]);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> &u_ = u.get(iX, iY, iZ);
                Array<T, 3> &gradMu_ = gradMu.get(iX, iY, iZ);
                Array<T, 3> &gradC_ = gradC.get(iX, iY, iZ);
                T invRho = (T)1 / rho.get(iX, iY, iZ);

                FiniteDifference<T>(mu, iX, iY, iZ).centralGradient(gradMu_);
                laplaceMu.get(iX, iY, iZ) = FiniteDifference<T>(mu, iX, iY, iZ).laplacian();
                // Use templates to compute order-0 and order-1 moment of g.
                momentTemplates<T, Descriptor>::get_rhoBar_j(
                    g.get(iX, iY, iZ), p1.get(iX, iY, iZ), u_);
                u_ *= invRho / RT;
                u_ -= gradMu_ * ((T)0.5 * invRho * C.get(iX, iY, iZ));
                p1.get(iX, iY, iZ) += 0.5 * RT * (rho_h - rho_l)
                                      * VectorTemplateImpl<T, 3>::scalarProduct(u_, gradC_);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>
    *Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>::clone() const
{
    return new Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void Compute_gradMu_laplaceMu_u_p1_processor<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // f
    modified[1] = modif::nothing;          // g
    modified[2] = modif::nothing;          // C
    modified[3] = modif::nothing;          // rho
    modified[4] = modif::nothing;          // gradC
    modified[5] = modif::nothing;          // mu
    modified[6] = modif::staticVariables;  // gradMu
    modified[7] = modif::staticVariables;  // laplaceMu
    modified[8] = modif::staticVariables;  // u
    modified[9] = modif::staticVariables;  // p1
}

/* *************** HeLeeCollisionProcessor ***************** */

template <typename T, template <typename U> class Descriptor>
HeLeeCollisionProcessor<T, Descriptor>::HeLeeCollisionProcessor(
    T rho_h_, T rho_l_, T tau_h_, T tau_l_, T M_, T RT_, bool initialize_) :
    rho_h(rho_h_),
    rho_l(rho_l_),
    tau_h(tau_h_),
    tau_l(tau_l_),
    M(M_),
    RT(RT_),
    initialize(initialize_)
{ }

template <typename T, template <typename U> class Descriptor>
void HeLeeCollisionProcessor<T, Descriptor>::computeAdvectionTerms(
    ScalarField3D<T> const &C, T &adv_gradC, T &bias_adv_gradC, plint iX, plint iY, plint iZ,
    plint iPop)
{
    typedef Descriptor<T> D;
    T diff_p2 = C.get(iX + 2 * D::c[iPop][0], iY + 2 * D::c[iPop][1], iZ + 2 * D::c[iPop][2])
                - C.get(iX + D::c[iPop][0], iY + D::c[iPop][1], iZ + D::c[iPop][2]);
    T diff_p1 =
        C.get(iX + D::c[iPop][0], iY + D::c[iPop][1], iZ + D::c[iPop][2]) - C.get(iX, iY, iZ);
    T diff_p0 =
        C.get(iX, iY, iZ) - C.get(iX - D::c[iPop][0], iY - D::c[iPop][1], iZ - D::c[iPop][2]);
    adv_gradC = diff_p1 - 0.5 * (diff_p1 - diff_p0);
    bias_adv_gradC = diff_p1 - 0.25 * (diff_p2 - diff_p0);
}

template <typename T, template <typename U> class Descriptor>
void HeLeeCollisionProcessor<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    typedef Descriptor<T> D;
    BlockLattice3D<T, Descriptor> &f = *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    BlockLattice3D<T, Descriptor> &g = *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[1]);
    ScalarField3D<T> &C = *dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    ScalarField3D<T> &rho = *dynamic_cast<ScalarField3D<T> *>(blocks[3]);
    TensorField3D<T, 3> &gradC = *dynamic_cast<TensorField3D<T, 3> *>(blocks[4]);
    ScalarField3D<T> &mu = *dynamic_cast<ScalarField3D<T> *>(blocks[5]);
    TensorField3D<T, 3> &gradMu = *dynamic_cast<TensorField3D<T, 3> *>(blocks[6]);
    ScalarField3D<T> &laplaceMu = *dynamic_cast<ScalarField3D<T> *>(blocks[7]);
    TensorField3D<T, 3> &u = *dynamic_cast<TensorField3D<T, 3> *>(blocks[8]);
    ScalarField3D<T> &p1 = *dynamic_cast<ScalarField3D<T> *>(blocks[9]);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> &u_ = u.get(iX, iY, iZ);
                T &C_ = C.get(iX, iY, iZ);
                T &rho_ = rho.get(iX, iY, iZ);
                T &laplaceMu_ = laplaceMu.get(iX, iY, iZ);
                T &p1_ = p1.get(iX, iY, iZ);
                Cell<T, Descriptor> &f_ = f.get(iX, iY, iZ);
                Cell<T, Descriptor> &g_ = g.get(iX, iY, iZ);

                T tau = C_ * (tau_h - tau_l) + tau_l;

                Array<T, 3> biasGradC, biasGradMu, gradP, biasGradP;
                FiniteDifference<T>(C, iX, iY, iZ).biasedGradient(biasGradC);
                FiniteDifference<T>(mu, iX, iY, iZ).biasedGradient(biasGradMu);
                FiniteDifference<T>(p1, iX, iY, iZ).centralGradient(gradP);
                FiniteDifference<T>(p1, iX, iY, iZ).biasedGradient(biasGradP);

                T uGradC = VectorTemplateImpl<T, 3>::scalarProduct(u_, gradC.get(iX, iY, iZ));
                T uBiasGradC = VectorTemplateImpl<T, 3>::scalarProduct(u_, biasGradC);
                T uGradMu = VectorTemplateImpl<T, 3>::scalarProduct(u_, gradMu.get(iX, iY, iZ));
                T uBiasGradMu = VectorTemplateImpl<T, 3>::scalarProduct(u_, biasGradMu);
                T uGradP = VectorTemplateImpl<T, 3>::scalarProduct(u_, gradP);
                T uBiasGradP = VectorTemplateImpl<T, 3>::scalarProduct(u_, biasGradP);
                T uSqr = VectorTemplateImpl<T, 3>::normSqr(u_);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    T ci_u = D::c[iPop][0] * u_[0] + D::c[iPop][1] * u_[1] + D::c[iPop][2] * u_[2];
                    T ti = Descriptor<T>::t[iPop];
                    T gamma0 = ti;
                    T gammaBar = 3. * ci_u + 4.5 * ci_u * ci_u - 1.5 * uSqr;
                    T gamma = ti * (1. + gammaBar);

                    T adv_gradC, bias_adv_gradC;
                    computeAdvectionTerms(C, adv_gradC, bias_adv_gradC, iX, iY, iZ, iPop);
                    adv_gradC -= uGradC;
                    bias_adv_gradC -= uBiasGradC;

                    T adv_gradP, bias_adv_gradP;
                    computeAdvectionTerms(p1, adv_gradP, bias_adv_gradP, iX, iY, iZ, iPop);
                    adv_gradP -= uGradP;
                    bias_adv_gradP -= uBiasGradP;

                    T adv_gradMu, bias_adv_gradMu;
                    computeAdvectionTerms(mu, adv_gradMu, bias_adv_gradMu, iX, iY, iZ, iPop);
                    adv_gradMu -= uGradMu;
                    adv_gradMu *= C_;
                    bias_adv_gradMu -= uBiasGradMu;
                    bias_adv_gradMu *= C_;

                    T fieq =
                        ti * C_ * (1 + gammaBar)
                        - 0.5 * gamma * (adv_gradC - 1. / RT * (adv_gradP + adv_gradMu) * C_ / rho_)
                        - 0.5 * gamma * M * laplaceMu_;

                    T gieq = ti * (p1_ + rho_ * RT * gammaBar)
                             - 0.5 * RT * adv_gradC * (rho_h - rho_l) * (gamma - gamma0)
                             + 0.5 * gamma * adv_gradMu;
                    if (initialize) {
                        f_[iPop] = fieq;
                        g_[iPop] = gieq;
                    } else {
                        f_[iPop] +=
                            -(f_[iPop] - fieq) * 1. / (tau + 0.5)
                            + gamma
                                  * (bias_adv_gradC
                                     - (bias_adv_gradP + bias_adv_gradMu) * 1. / RT * C_ / rho_)
                            + M * laplaceMu_ * gamma;
                        g_[iPop] += -(g_[iPop] - gieq) * 1. / (tau + 0.5)
                                    + bias_adv_gradC * (rho_h - rho_l) * RT * (gamma - gamma0)
                                    - bias_adv_gradMu * gamma;
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
HeLeeCollisionProcessor<T, Descriptor> *HeLeeCollisionProcessor<T, Descriptor>::clone() const
{
    return new HeLeeCollisionProcessor<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void HeLeeCollisionProcessor<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // f
    modified[1] = modif::staticVariables;  // g
    modified[2] = modif::nothing;          // C
    modified[3] = modif::nothing;          // rho
    modified[4] = modif::nothing;          // gradC
    modified[5] = modif::nothing;          // mu
    modified[6] = modif::nothing;          // gradMu
    modified[7] = modif::nothing;          // laplaceMu
    modified[8] = modif::nothing;          // u
    modified[9] = modif::nothing;          // p1
}

}  // namespace plb

#endif  // HE_LEE_PROCESSOR_3D_HH
