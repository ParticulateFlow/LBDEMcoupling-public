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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef SMAGORINSKY_GENERICS_3D_H
#define SMAGORINSKY_GENERICS_3D_H

namespace plb {

template <typename T, template <typename U> class Descriptor, class SmagoFunction>
class StaticSmagorinskyFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    StaticSmagorinskyFunctional3D(SmagoFunction smagoFunction_, T cSmago0_) :
        smagoFunction(smagoFunction_), cSmago0(cSmago0_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
    {
        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T cSmago = smagoFunction(
                        iX + relativeOffset.x, iY + relativeOffset.y, iZ + relativeOffset.z,
                        cSmago0);
                    T omega0 = lattice.get(iX, iY, iZ).getDynamics().getOmega();
                    lattice.attributeDynamics(
                        iX, iY, iZ,
                        cloneAndInsertAtTopDynamics(
                            lattice.get(iX, iY, iZ).getDynamics(),
                            new SmagorinskyDynamics<T, Descriptor>(
                                new NoDynamics<T, Descriptor>, omega0, cSmago)));
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        // Composite dynamics needs to be instantiated everywhere, including envelope.
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual StaticSmagorinskyFunctional3D<T, Descriptor, SmagoFunction> *clone() const
    {
        return new StaticSmagorinskyFunctional3D<T, Descriptor, SmagoFunction>(*this);
    }

private:
    SmagoFunction smagoFunction;
    T cSmago0;
};

template <typename T>
T constOmegaFromOmega0(plint iX, plint iY, plint iZ, T omega0)
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
void instantiateStaticSmagorinsky(BlockLattice3D<T, Descriptor> &lattice, Box3D domain, T cSmago)
{
    instantiateStaticSmagorinsky(lattice, domain, constOmegaFromOmega0<T>, cSmago);
}

template <typename T, template <typename U> class Descriptor>
void instantiateStaticSmagorinsky(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, T cSmago)
{
    instantiateStaticSmagorinsky(lattice, domain, constOmegaFromOmega0<T>, cSmago);
}

template <typename T, template <typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(
    BlockLattice3D<T, Descriptor> &lattice, Box3D domain, SmagoFunction smagoFunction, T cSmago0)
{
    applyProcessingFunctional(
        new StaticSmagorinskyFunctional3D<T, Descriptor, SmagoFunction>(smagoFunction, cSmago0),
        domain, lattice);
}

template <typename T, template <typename U> class Descriptor, class SmagoFunction>
void instantiateStaticSmagorinsky(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, SmagoFunction smagoFunction,
    T cSmago0)
{
    applyProcessingFunctional(
        new StaticSmagorinskyFunctional3D<T, Descriptor, SmagoFunction>(smagoFunction, cSmago0),
        domain, lattice);
}

}  // namespace plb

#endif  // SMAGORINSKY_GENERICS_3D_H
