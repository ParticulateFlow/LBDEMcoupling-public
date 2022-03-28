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

/** \file
 * Interface for dataProcessing steps -- implementation file.
 */

#include "atomicBlock/dataProcessor3D.h"

#include "core/util.h"

namespace plb {

////////////////////// Class DataProcessor3D /////////////////

/** This method returns the maximum extent of the processor, over any
 *  direction. By extent one means the the size of the neighborhood
 *  for non-local accesses. For example, if the processor implements
 *  a second-order accurate laplace operator, the value returned
 *  by extent() is 1, because the second derivative is evaluated
 *  with help of the -1 2 -1 stencil, requiring one left and one
 *  right neighbor.
 *  The default implementation of this method returns an extent of 1,
 *  based on the assumption that most LB stuff is somehow based on
 *  nearest-neighbor interaction.
 *  This is a bit dangerous though, as one easily forgets to override
 *  the method in case of larger-than-nearest-neighbor Processors.
 *  Still, LatticeProcessors are much easier to write when there's only
 *  one method to override, so we'll let it be this way.
 */
plint DataProcessor3D::extent() const
{
    return 1;
}

/** By default, this method assumes a symmetric neighborhood relation
 *  and refers to the non-directed version of extent().
 */
plint DataProcessor3D::extent(int direction) const
{
    return extent();
}

/** Return -1 as default to help transition period as some
 *  data processors have no ID.
 **/
int DataProcessor3D::getStaticId() const
{
    return -1;
}

////////////////////// Class DataProcessorGenerator3D /////////////////

DataProcessorGenerator3D::~DataProcessorGenerator3D() { }

BlockDomain::DomainT DataProcessorGenerator3D::appliesTo() const
{
    return BlockDomain::bulk;
}

void DataProcessorGenerator3D::rescale(double dxScale_, double dtScale_) { }

/** \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 */
void DataProcessorGenerator3D::setscale(int dxScale_, int dtScale_) { }

void DataProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/** Return -1 as default to help transition period as some
 *  data processors have no ID.
 **/
int DataProcessorGenerator3D::getStaticId() const
{
    return -1;
}

/** Default action does nothing to help transition period. **/
void DataProcessorGenerator3D::serialize(Box3D &domain, std::string &data) const { }

////////////////////// Class BoxedDataProcessorGenerator3D /////////////////

BoxedDataProcessorGenerator3D::BoxedDataProcessorGenerator3D(Box3D domain_) : domain(domain_) { }

void BoxedDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    domain = domain.shift(deltaX, deltaY, deltaZ);
}

void BoxedDataProcessorGenerator3D::multiply(plint scale)
{
    domain = domain.multiply(scale);
}

void BoxedDataProcessorGenerator3D::divide(plint scale)
{
    domain = domain.divide(scale);
}

bool BoxedDataProcessorGenerator3D::extract(Box3D subDomain)
{
    Box3D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    } else {
        return false;
    }
}

Box3D BoxedDataProcessorGenerator3D::getDomain() const
{
    return domain;
}

void BoxedDataProcessorGenerator3D::serialize(Box3D &domain_, std::string &data) const
{
    domain_ = domain;
}

////////////////////// Class MultiBoxedDataProcessorGenerator3D /////////////////

MultiBoxedDataProcessorGenerator3D::MultiBoxedDataProcessorGenerator3D(
    std::vector<Box3D> const &domains_) :
    domains(domains_)
{ }

void MultiBoxedDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].shift(deltaX, deltaY, deltaZ);
    }
}

void MultiBoxedDataProcessorGenerator3D::multiply(plint scale)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].multiply(scale);
    }
}

void MultiBoxedDataProcessorGenerator3D::divide(plint scale)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].divide(scale);
    }
}

bool MultiBoxedDataProcessorGenerator3D::extract(Box3D subDomain)
{
    std::vector<Box3D> intersections;
    for (pluint i = 0; i < domains.size(); ++i) {
        Box3D intersection;
        if (intersect(domains[i], subDomain, intersection)) {
            intersections.push_back(intersection);
        }
    }
    if (intersections.empty()) {
        return false;
    } else {
        intersections.swap(domains);
        return true;
    }
}

std::vector<Box3D> const &MultiBoxedDataProcessorGenerator3D::getDomains() const
{
    return domains;
}

void MultiBoxedDataProcessorGenerator3D::serialize(Box3D &domain_, std::string &data) const { }

////////////////////// Class ReductiveDataProcessorGenerator3D /////////////////

ReductiveDataProcessorGenerator3D::ReductiveDataProcessorGenerator3D() : dxScale(0), dtScale(0) { }

ReductiveDataProcessorGenerator3D::~ReductiveDataProcessorGenerator3D() { }

BlockDomain::DomainT ReductiveDataProcessorGenerator3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** Default implementation does nothing, to reflect the fact that many
 *  processor generators are invariant under change of scale. It is
 *  however important to remember to redefine this method for data
 *  processors which dependent on physical units, because otherwise
 *  the processor yields erroneous results in a MultiGrid environment.
 *
 *  \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 */
void ReductiveDataProcessorGenerator3D::rescale(double dxScale, double dtScale) { }

/** \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 *
 *  If you override this function, you should always remember to invoke the
 *  original implementation as well, so the scale factors are properly
 *  set in the generator.
 */
void ReductiveDataProcessorGenerator3D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int ReductiveDataProcessorGenerator3D::getDxScale() const
{
    return dxScale;
}

int ReductiveDataProcessorGenerator3D::getDtScale() const
{
    return dtScale;
}

/** Default implementation: constant dimensions.
 **/
void ReductiveDataProcessorGenerator3D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

/** Default implementation: constant dimensions.
 **/
void ReductiveDataProcessorGenerator3D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDataProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

void ReductiveDataProcessorGenerator3D::serialize(Box3D &domain, std::string &data) const { }

////////////////////// Class BoxedReductiveDataProcessorGenerator3D /////////////////

BoxedReductiveDataProcessorGenerator3D::BoxedReductiveDataProcessorGenerator3D(Box3D domain_) :
    domain(domain_)
{ }

void BoxedReductiveDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    domain = domain.shift(deltaX, deltaY, deltaZ);
}

void BoxedReductiveDataProcessorGenerator3D::multiply(plint scale)
{
    domain = domain.multiply(scale);
}

void BoxedReductiveDataProcessorGenerator3D::divide(plint scale)
{
    domain = domain.divide(scale);
}

bool BoxedReductiveDataProcessorGenerator3D::extract(Box3D subDomain)
{
    Box3D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    } else {
        return false;
    }
}

Box3D BoxedReductiveDataProcessorGenerator3D::getDomain() const
{
    return domain;
}

void BoxedReductiveDataProcessorGenerator3D::serialize(Box3D &domain_, std::string &data) const
{
    domain_ = domain;
}

////////////////////// Class MultiBoxedReductiveDataProcessorGenerator3D /////////////////

MultiBoxedReductiveDataProcessorGenerator3D::MultiBoxedReductiveDataProcessorGenerator3D(
    std::vector<Box3D> const &domains_) :
    domains(domains_)
{ }

void MultiBoxedReductiveDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].shift(deltaX, deltaY, deltaZ);
    }
}

void MultiBoxedReductiveDataProcessorGenerator3D::multiply(plint scale)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].multiply(scale);
    }
}

void MultiBoxedReductiveDataProcessorGenerator3D::divide(plint scale)
{
    for (pluint i = 0; i < domains.size(); ++i) {
        domains[i] = domains[i].divide(scale);
    }
}

bool MultiBoxedReductiveDataProcessorGenerator3D::extract(Box3D subDomain)
{
    std::vector<Box3D> intersections;
    for (pluint i = 0; i < domains.size(); ++i) {
        Box3D intersection;
        if (intersect(domains[i], subDomain, intersection)) {
            intersections.push_back(intersection);
        }
    }
    if (intersections.empty()) {
        return false;
    } else {
        intersections.swap(domains);
        return true;
    }
}

std::vector<Box3D> const &MultiBoxedReductiveDataProcessorGenerator3D::getDomains() const
{
    return domains;
}

void MultiBoxedReductiveDataProcessorGenerator3D::serialize(Box3D &domain_, std::string &data) const
{ }

////////////////////// Class DottedDataProcessorGenerator3D /////////////////

DottedDataProcessorGenerator3D::DottedDataProcessorGenerator3D(DotList3D const &dots_) : dots(dots_)
{ }

void DottedDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    dots = dots.shift(deltaX, deltaY, deltaZ);
}

void DottedDataProcessorGenerator3D::multiply(plint scale)
{
    dots = dots.multiply(scale);
}

void DottedDataProcessorGenerator3D::divide(plint scale)
{
    dots = dots.divide(scale);
}

bool DottedDataProcessorGenerator3D::extract(Box3D subDomain)
{
    DotList3D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    } else {
        return false;
    }
}

DotList3D const &DottedDataProcessorGenerator3D::getDotList() const
{
    return dots;
}

////////////////////// Class DottedReductiveDataProcessorGenerator3D /////////////////

DottedReductiveDataProcessorGenerator3D::DottedReductiveDataProcessorGenerator3D(
    DotList3D const &dots_) :
    dots(dots_)
{ }

void DottedReductiveDataProcessorGenerator3D::shift(plint deltaX, plint deltaY, plint deltaZ)
{
    dots = dots.shift(deltaX, deltaY, deltaZ);
}

void DottedReductiveDataProcessorGenerator3D::multiply(plint scale)
{
    dots = dots.multiply(scale);
}

void DottedReductiveDataProcessorGenerator3D::divide(plint scale)
{
    dots = dots.divide(scale);
}

bool DottedReductiveDataProcessorGenerator3D::extract(Box3D subDomain)
{
    DotList3D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    } else {
        return false;
    }
}

DotList3D const &DottedReductiveDataProcessorGenerator3D::getDotList() const
{
    return dots;
}

}  // namespace plb
