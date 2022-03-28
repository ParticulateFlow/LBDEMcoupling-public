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

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class ReductiveBoxProcessingFunctional2D ************************* */

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT ReductiveBoxProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void ReductiveBoxProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void ReductiveBoxProcessingFunctional2D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveBoxProcessingFunctional2D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveBoxProcessingFunctional2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/** Return a default value of -1. This is to help transition for legacy
 *  code that has not yet implemented this method.
 ***/
int ReductiveBoxProcessingFunctional2D::getStaticId() const
{
    return -1;
}

void ReductiveBoxProcessingFunctional2D::serialize(std::string &data) const { }

void ReductiveBoxProcessingFunctional2D::unserialize(std::string &data) { }

/* *************** Class PlainReductiveBoxProcessingFunctional2D ************************* */

BlockStatistics const &PlainReductiveBoxProcessingFunctional2D::getStatistics() const
{
    return statistics;
}

BlockStatistics &PlainReductiveBoxProcessingFunctional2D::getStatistics()
{
    return statistics;
}

/* *************** Class ReductiveBoxProcessor2D ************************************ */

ReductiveBoxProcessor2D::ReductiveBoxProcessor2D(
    ReductiveBoxProcessingFunctional2D *functional_, Box2D domain_,
    std::vector<AtomicBlock2D *> atomicBlocks_) :
    functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

Box2D ReductiveBoxProcessor2D::getDomain() const
{
    return domain;
}

void ReductiveBoxProcessor2D::process()
{
    functional->processGenericBlocks(domain, atomicBlocks);
    functional->getStatistics().evaluate();
}

ReductiveBoxProcessor2D *ReductiveBoxProcessor2D::clone() const
{
    return new ReductiveBoxProcessor2D(*this);
}

int ReductiveBoxProcessor2D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class ReductiveBoxProcessorGenerator2D *************************** */

ReductiveBoxProcessorGenerator2D::ReductiveBoxProcessorGenerator2D(
    ReductiveBoxProcessingFunctional2D *functional_, Box2D domain) :
    BoxedReductiveDataProcessorGenerator2D(domain), functional(functional_)
{
    // Must be non-null, because it is then used without further checks.
    PLB_ASSERT(functional);
}

ReductiveBoxProcessorGenerator2D::~ReductiveBoxProcessorGenerator2D()
{
    delete functional;
}

ReductiveBoxProcessorGenerator2D::ReductiveBoxProcessorGenerator2D(
    ReductiveBoxProcessorGenerator2D const &rhs) :
    BoxedReductiveDataProcessorGenerator2D(rhs), functional(rhs.functional->clone())
{ }

ReductiveBoxProcessorGenerator2D &ReductiveBoxProcessorGenerator2D::operator=(
    ReductiveBoxProcessorGenerator2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT ReductiveBoxProcessorGenerator2D::appliesTo() const
{
    return functional->appliesTo();
}

void ReductiveBoxProcessorGenerator2D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void ReductiveBoxProcessorGenerator2D::getDimensionsX(std::vector<int> &dimensions) const
{
    functional->getDimensionsX(dimensions);
}

void ReductiveBoxProcessorGenerator2D::getDimensionsT(std::vector<int> &dimensions) const
{
    functional->getDimensionsT(dimensions);
}

void ReductiveBoxProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void ReductiveBoxProcessorGenerator2D::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor2D *ReductiveBoxProcessorGenerator2D::generate(
    std::vector<AtomicBlock2D *> atomicBlocks)
{
    // Don't clone functional. Given that the functional contains the BlockStatistics object,
    //   everybody must point to the same instance.
    return new ReductiveBoxProcessor2D(functional, this->getDomain(), atomicBlocks);
}

ReductiveBoxProcessorGenerator2D *ReductiveBoxProcessorGenerator2D::clone() const
{
    return new ReductiveBoxProcessorGenerator2D(*this);
}

BlockStatistics const &ReductiveBoxProcessorGenerator2D::getStatistics() const
{
    return functional->getStatistics();
}

BlockStatistics &ReductiveBoxProcessorGenerator2D::getStatistics()
{
    return functional->getStatistics();
}

ReductiveBoxProcessingFunctional2D const &ReductiveBoxProcessorGenerator2D::getFunctional() const
{
    return *functional;
}

void ReductiveBoxProcessorGenerator2D::serialize(Box2D &domain, std::string &data) const
{
    BoxedReductiveDataProcessorGenerator2D::serialize(domain, data);
    functional->serialize(data);
}

int ReductiveBoxProcessorGenerator2D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class ReductiveDotProcessingFunctional2D ************************* */

/** Operation is not executed on envelope by default. **/
BlockDomain::DomainT ReductiveDotProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void ReductiveDotProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void ReductiveDotProcessingFunctional2D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDotProcessingFunctional2D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDotProcessingFunctional2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class PlainReductiveDotProcessingFunctional2D ************************* */

BlockStatistics const &PlainReductiveDotProcessingFunctional2D::getStatistics() const
{
    return statistics;
}

BlockStatistics &PlainReductiveDotProcessingFunctional2D::getStatistics()
{
    return statistics;
}

/* *************** Class ReductiveDotProcessor2D ************************************ */

ReductiveDotProcessor2D::ReductiveDotProcessor2D(
    ReductiveDotProcessingFunctional2D *functional_, DotList2D const &dotList_,
    std::vector<AtomicBlock2D *> atomicBlocks_) :
    functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

DotList2D const &ReductiveDotProcessor2D::getDotList() const
{
    return dotList;
}

void ReductiveDotProcessor2D::process()
{
    functional->processGenericBlocks(dotList, atomicBlocks);
}

ReductiveDotProcessor2D *ReductiveDotProcessor2D::clone() const
{
    return new ReductiveDotProcessor2D(*this);
}

/* *************** Class ReductiveDotProcessorGenerator2D *************************** */

ReductiveDotProcessorGenerator2D::ReductiveDotProcessorGenerator2D(
    ReductiveDotProcessingFunctional2D *functional_, DotList2D const &dotList) :
    DottedReductiveDataProcessorGenerator2D(dotList), functional(functional_)
{ }

ReductiveDotProcessorGenerator2D::~ReductiveDotProcessorGenerator2D()
{
    delete functional;
}

ReductiveDotProcessorGenerator2D::ReductiveDotProcessorGenerator2D(
    ReductiveDotProcessorGenerator2D const &rhs) :
    DottedReductiveDataProcessorGenerator2D(rhs), functional(rhs.functional->clone())
{ }

ReductiveDotProcessorGenerator2D &ReductiveDotProcessorGenerator2D::operator=(
    ReductiveDotProcessorGenerator2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT ReductiveDotProcessorGenerator2D::appliesTo() const
{
    return functional->appliesTo();
}

void ReductiveDotProcessorGenerator2D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void ReductiveDotProcessorGenerator2D::getDimensionsX(std::vector<int> &dimensions) const
{
    functional->getDimensionsX(dimensions);
}

void ReductiveDotProcessorGenerator2D::getDimensionsT(std::vector<int> &dimensions) const
{
    functional->getDimensionsT(dimensions);
}

void ReductiveDotProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void ReductiveDotProcessorGenerator2D::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor2D *ReductiveDotProcessorGenerator2D::generate(
    std::vector<AtomicBlock2D *> atomicBlocks)
{
    return new ReductiveDotProcessor2D(functional, this->getDotList(), atomicBlocks);
}

ReductiveDotProcessorGenerator2D *ReductiveDotProcessorGenerator2D::clone() const
{
    return new ReductiveDotProcessorGenerator2D(*this);
}

BlockStatistics const &ReductiveDotProcessorGenerator2D::getStatistics() const
{
    return functional->getStatistics();
}

BlockStatistics &ReductiveDotProcessorGenerator2D::getStatistics()
{
    return functional->getStatistics();
}

ReductiveDotProcessingFunctional2D const &ReductiveDotProcessorGenerator2D::getFunctional() const
{
    return *functional;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional2D ********** */

/** Operation is not executed on envelope by default. **/
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void BoundedReductiveBoxProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void BoundedReductiveBoxProcessingFunctional2D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void BoundedReductiveBoxProcessingFunctional2D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void BoundedReductiveBoxProcessingFunctional2D::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

ReductiveBoxProcessingFunctional2D *BoundedReductiveBoxProcessingFunctional2D::getBulkProcessor()
    const
{
    return new BulkWrapperFunctional(this->clone());
}

ReductiveBoxProcessingFunctional2D *BoundedReductiveBoxProcessingFunctional2D::getEdgeProcessor(
    int direction, int orientation) const
{
    return new EdgeWrapperFunctional(this->clone(), direction, orientation);
}

ReductiveBoxProcessingFunctional2D *BoundedReductiveBoxProcessingFunctional2D::getCornerProcessor(
    int normalX, int normalY) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY);
}

BlockStatistics const &BoundedReductiveBoxProcessingFunctional2D::getStatistics() const
{
    return statistics;
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional2D::getStatistics()
{
    return statistics;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::BulkWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_) :
    boundedFunctional(boundedFunctional_)
{ }

BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::BulkWrapperFunctional(
    BulkWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone())
{ }

BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::~BulkWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::operator=(
        BulkWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional2D::BulkWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

/* *************** Class BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_, int direction_,
    int orientation_) :
    boundedFunctional(boundedFunctional_), direction(direction_), orientation(orientation_)
{ }

BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    EdgeWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::~EdgeWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::operator=(
        EdgeWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processEdgeGeneric(direction, orientation, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional2D::EdgeWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

/* *************** Class BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::CornerWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional2D *boundedFunctional_, int normalX_, int normalY_) :
    boundedFunctional(boundedFunctional_), normalX(normalX_), normalY(normalY_)
{ }

BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::CornerWrapperFunctional(
    CornerWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()), normalX(rhs.normalX), normalY(rhs.normalY)
{ }

BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::~CornerWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::operator=(
        CornerWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    return *this;
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processCornerGeneric(normalX, normalY, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional2D::CornerWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

void BoundedReductiveBoxProcessingFunctional2D::getGenerators(
    Box2D const &fullDomain, plint boundaryWidth,
    std::vector<ReductiveBoxProcessorGenerator2D *> &generators)
{
    generators.resize(9);
    BlockSurface2D surf(fullDomain, boundaryWidth);
    generators[0] = new ReductiveBoxProcessorGenerator2D(getBulkProcessor(), surf.bulk());

    generators[1] = new ReductiveBoxProcessorGenerator2D(getEdgeProcessor(0, -1), surf.edge0N());
    generators[2] = new ReductiveBoxProcessorGenerator2D(getEdgeProcessor(0, +1), surf.edge0P());
    generators[3] = new ReductiveBoxProcessorGenerator2D(getEdgeProcessor(1, -1), surf.edge1N());
    generators[4] = new ReductiveBoxProcessorGenerator2D(getEdgeProcessor(1, +1), surf.edge1P());

    generators[5] =
        new ReductiveBoxProcessorGenerator2D(getCornerProcessor(-1, -1), surf.cornerNN());
    generators[6] =
        new ReductiveBoxProcessorGenerator2D(getCornerProcessor(+1, -1), surf.cornerPN());
    generators[7] =
        new ReductiveBoxProcessorGenerator2D(getCornerProcessor(-1, +1), surf.cornerNP());
    generators[8] =
        new ReductiveBoxProcessorGenerator2D(getCornerProcessor(+1, +1), surf.cornerPP());
}

}  // namespace plb
