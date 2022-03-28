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
#include "atomicBlock/dataProcessor2D.h"

#include "core/util.h"

namespace plb {

////////////////////// Class DataProcessor2D /////////////////

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
plint DataProcessor2D::extent() const
{
    return 1;
}

/** By default, this method assumes a symmetric neighborhood relation
 *  and refers to the non-directed version of extent().
 */
plint DataProcessor2D::extent(int direction) const
{
    return extent();
}

/** Return -1 as default to help transition period as some
 *  data processors have no ID.
 **/
int DataProcessor2D::getStaticId() const
{
    return -1;
}

////////////////////// Class DataProcessorGenerator2D /////////////////

DataProcessorGenerator2D::~DataProcessorGenerator2D() { }

BlockDomain::DomainT DataProcessorGenerator2D::appliesTo() const
{
    return BlockDomain::bulk;
}

void DataProcessorGenerator2D::rescale(double dxScale_, double dtScale_) { }

/** \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 */
void DataProcessorGenerator2D::setscale(int dxScale_, int dtScale_) { }

void DataProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
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
int DataProcessorGenerator2D::getStaticId() const
{
    return -1;
}

/** Default action does nothing to help transition period. **/
void DataProcessorGenerator2D::serialize(Box2D &domain, std::string &data) const { }

////////////////////// Class BoxedDataProcessorGenerator2D /////////////////

BoxedDataProcessorGenerator2D::BoxedDataProcessorGenerator2D(Box2D domain_) : domain(domain_) { }

void BoxedDataProcessorGenerator2D::shift(plint deltaX, plint deltaY)
{
    domain = domain.shift(deltaX, deltaY);
}

void BoxedDataProcessorGenerator2D::multiply(plint scale)
{
    domain = domain.multiply(scale);
}

void BoxedDataProcessorGenerator2D::divide(plint scale)
{
    domain = domain.divide(scale);
}

bool BoxedDataProcessorGenerator2D::extract(Box2D subDomain)
{
    Box2D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    } else {
        return false;
    }
}

Box2D BoxedDataProcessorGenerator2D::getDomain() const
{
    return domain;
}

void BoxedDataProcessorGenerator2D::serialize(Box2D &domain_, std::string &data) const
{
    domain_ = domain;
}

////////////////////// Class ReductiveDataProcessorGenerator2D /////////////////

ReductiveDataProcessorGenerator2D::ReductiveDataProcessorGenerator2D() : dxScale(0), dtScale(0) { }

ReductiveDataProcessorGenerator2D::~ReductiveDataProcessorGenerator2D() { }

BlockDomain::DomainT ReductiveDataProcessorGenerator2D::appliesTo() const
{
    return BlockDomain::bulk;
}

void ReductiveDataProcessorGenerator2D::rescale(double dxScale, double dtScale) { }

/** \param dxScale Scale factor for space scale dx.
 *  \param dtScale Scale factor for time scale dt.
 *
 *  If you override this function, you should always remember to invoke the
 *  original implementation as well, so the scale factors are properly
 *  set in the generator.
 */
void ReductiveDataProcessorGenerator2D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int ReductiveDataProcessorGenerator2D::getDxScale() const
{
    return dxScale;
}

int ReductiveDataProcessorGenerator2D::getDtScale() const
{
    return dtScale;
}

/** Default implementation: constant dimensions.
 **/
void ReductiveDataProcessorGenerator2D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

/** Default implementation: constant dimensions.
 **/
void ReductiveDataProcessorGenerator2D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDataProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

void ReductiveDataProcessorGenerator2D::serialize(Box2D &domain, std::string &data) const { }

////////////////////// Class BoxedReductiveDataProcessorGenerator2D /////////////////

BoxedReductiveDataProcessorGenerator2D::BoxedReductiveDataProcessorGenerator2D(Box2D domain_) :
    domain(domain_)
{ }

void BoxedReductiveDataProcessorGenerator2D::shift(plint deltaX, plint deltaY)
{
    domain = domain.shift(deltaX, deltaY);
}

void BoxedReductiveDataProcessorGenerator2D::multiply(plint scale)
{
    domain = domain.multiply(scale);
}

void BoxedReductiveDataProcessorGenerator2D::divide(plint scale)
{
    domain = domain.divide(scale);
}

bool BoxedReductiveDataProcessorGenerator2D::extract(Box2D subDomain)
{
    Box2D intersection;
    if (intersect(domain, subDomain, intersection)) {
        domain = intersection;
        return true;
    } else {
        return false;
    }
}

Box2D BoxedReductiveDataProcessorGenerator2D::getDomain() const
{
    return domain;
}

void BoxedReductiveDataProcessorGenerator2D::serialize(Box2D &domain_, std::string &data) const
{
    domain_ = domain;
}

////////////////////// Class DottedDataProcessorGenerator2D /////////////////

DottedDataProcessorGenerator2D::DottedDataProcessorGenerator2D(DotList2D const &dots_) : dots(dots_)
{ }

void DottedDataProcessorGenerator2D::shift(plint deltaX, plint deltaY)
{
    dots = dots.shift(deltaX, deltaY);
}

void DottedDataProcessorGenerator2D::multiply(plint scale)
{
    dots = dots.multiply(scale);
}

void DottedDataProcessorGenerator2D::divide(plint scale)
{
    dots = dots.divide(scale);
}

bool DottedDataProcessorGenerator2D::extract(Box2D subDomain)
{
    DotList2D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    } else {
        return false;
    }
}

DotList2D const &DottedDataProcessorGenerator2D::getDotList() const
{
    return dots;
}

////////////////////// Class DottedReductiveDataProcessorGenerator2D /////////////////

DottedReductiveDataProcessorGenerator2D::DottedReductiveDataProcessorGenerator2D(
    DotList2D const &dots_) :
    dots(dots_)
{ }

void DottedReductiveDataProcessorGenerator2D::shift(plint deltaX, plint deltaY)
{
    dots = dots.shift(deltaX, deltaY);
}

void DottedReductiveDataProcessorGenerator2D::multiply(plint scale)
{
    dots = dots.multiply(scale);
}

void DottedReductiveDataProcessorGenerator2D::divide(plint scale)
{
    dots = dots.divide(scale);
}

bool DottedReductiveDataProcessorGenerator2D::extract(Box2D subDomain)
{
    DotList2D intersection;
    if (intersect(subDomain, dots, intersection)) {
        dots = intersection;
        return true;
    } else {
        return false;
    }
}

DotList2D const &DottedReductiveDataProcessorGenerator2D::getDotList() const
{
    return dots;
}

}  // namespace plb
