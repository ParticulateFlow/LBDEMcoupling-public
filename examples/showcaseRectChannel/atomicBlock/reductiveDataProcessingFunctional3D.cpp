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

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class ReductiveBoxProcessingFunctional3D ************************* */

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT ReductiveBoxProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void ReductiveBoxProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void ReductiveBoxProcessingFunctional3D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveBoxProcessingFunctional3D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveBoxProcessingFunctional3D::getModificationPattern(std::vector<bool> &isWritten) const
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
int ReductiveBoxProcessingFunctional3D::getStaticId() const
{
    return -1;
}

void ReductiveBoxProcessingFunctional3D::serialize(std::string &data) const { }

void ReductiveBoxProcessingFunctional3D::unserialize(std::string &data) { }

/* *************** Class PlainReductiveBoxProcessingFunctional3D ************************* */

BlockStatistics const &PlainReductiveBoxProcessingFunctional3D::getStatistics() const
{
    return statistics;
}

BlockStatistics &PlainReductiveBoxProcessingFunctional3D::getStatistics()
{
    return statistics;
}

/* *************** Class ReductiveBoxProcessor3D ************************************ */

ReductiveBoxProcessor3D::ReductiveBoxProcessor3D(
    ReductiveBoxProcessingFunctional3D *functional_, Box3D domain_,
    std::vector<AtomicBlock3D *> atomicBlocks_) :
    functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

Box3D ReductiveBoxProcessor3D::getDomain() const
{
    return domain;
}

void ReductiveBoxProcessor3D::process()
{
    functional->processGenericBlocks(domain, atomicBlocks);
    functional->getStatistics().evaluate();
}

ReductiveBoxProcessor3D *ReductiveBoxProcessor3D::clone() const
{
    return new ReductiveBoxProcessor3D(*this);
}

int ReductiveBoxProcessor3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class ReductiveBoxProcessorGenerator3D *************************** */

ReductiveBoxProcessorGenerator3D::ReductiveBoxProcessorGenerator3D(
    ReductiveBoxProcessingFunctional3D *functional_, Box3D domain) :
    BoxedReductiveDataProcessorGenerator3D(domain), functional(functional_)
{
    // Must be non-null, because it is then used without further checks.
    PLB_ASSERT(functional);
}

ReductiveBoxProcessorGenerator3D::~ReductiveBoxProcessorGenerator3D()
{
    delete functional;
}

ReductiveBoxProcessorGenerator3D::ReductiveBoxProcessorGenerator3D(
    ReductiveBoxProcessorGenerator3D const &rhs) :
    BoxedReductiveDataProcessorGenerator3D(rhs), functional(rhs.functional->clone())
{ }

ReductiveBoxProcessorGenerator3D &ReductiveBoxProcessorGenerator3D::operator=(
    ReductiveBoxProcessorGenerator3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT ReductiveBoxProcessorGenerator3D::appliesTo() const
{
    return functional->appliesTo();
}

void ReductiveBoxProcessorGenerator3D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void ReductiveBoxProcessorGenerator3D::getDimensionsX(std::vector<int> &dimensions) const
{
    functional->getDimensionsX(dimensions);
}

void ReductiveBoxProcessorGenerator3D::getDimensionsT(std::vector<int> &dimensions) const
{
    functional->getDimensionsT(dimensions);
}
void ReductiveBoxProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void ReductiveBoxProcessorGenerator3D::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor3D *ReductiveBoxProcessorGenerator3D::generate(
    std::vector<AtomicBlock3D *> atomicBlocks)
{
    // Don't clone functional. Given that the functional contains the BlockStatistics object,
    //   everybody must point to the same instance.
    return new ReductiveBoxProcessor3D(functional, this->getDomain(), atomicBlocks);
}

ReductiveBoxProcessorGenerator3D *ReductiveBoxProcessorGenerator3D::clone() const
{
    return new ReductiveBoxProcessorGenerator3D(*this);
}

BlockStatistics const &ReductiveBoxProcessorGenerator3D::getStatistics() const
{
    return functional->getStatistics();
}

BlockStatistics &ReductiveBoxProcessorGenerator3D::getStatistics()
{
    return functional->getStatistics();
}

ReductiveBoxProcessingFunctional3D const &ReductiveBoxProcessorGenerator3D::getFunctional() const
{
    return *functional;
}

void ReductiveBoxProcessorGenerator3D::serialize(Box3D &domain, std::string &data) const
{
    BoxedReductiveDataProcessorGenerator3D::serialize(domain, data);
    functional->serialize(data);
}

/** Return a default value of -1. This is to help transition for legacy
 *  code that has not yet implemented this method.
 ***/
int ReductiveBoxProcessorGenerator3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class ReductiveDotProcessingFunctional3D ************************* */

/** Operation is not executed on envelope by default. **/
BlockDomain::DomainT ReductiveDotProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void ReductiveDotProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void ReductiveDotProcessingFunctional3D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDotProcessingFunctional3D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void ReductiveDotProcessingFunctional3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class PlainReductiveDotProcessingFunctional3D ************************* */

BlockStatistics const &PlainReductiveDotProcessingFunctional3D::getStatistics() const
{
    return statistics;
}

BlockStatistics &PlainReductiveDotProcessingFunctional3D::getStatistics()
{
    return statistics;
}

/* *************** Class ReductiveDotProcessor3D ************************************ */

ReductiveDotProcessor3D::ReductiveDotProcessor3D(
    ReductiveDotProcessingFunctional3D *functional_, DotList3D const &dotList_,
    std::vector<AtomicBlock3D *> atomicBlocks_) :
    functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

DotList3D const &ReductiveDotProcessor3D::getDotList() const
{
    return dotList;
}

void ReductiveDotProcessor3D::process()
{
    functional->processGenericBlocks(dotList, atomicBlocks);
}

ReductiveDotProcessor3D *ReductiveDotProcessor3D::clone() const
{
    return new ReductiveDotProcessor3D(*this);
}

/* *************** Class ReductiveDotProcessorGenerator3D *************************** */

ReductiveDotProcessorGenerator3D::ReductiveDotProcessorGenerator3D(
    ReductiveDotProcessingFunctional3D *functional_, DotList3D const &dotList) :
    DottedReductiveDataProcessorGenerator3D(dotList), functional(functional_)
{ }

ReductiveDotProcessorGenerator3D::~ReductiveDotProcessorGenerator3D()
{
    delete functional;
}

ReductiveDotProcessorGenerator3D::ReductiveDotProcessorGenerator3D(
    ReductiveDotProcessorGenerator3D const &rhs) :
    DottedReductiveDataProcessorGenerator3D(rhs), functional(rhs.functional->clone())
{ }

ReductiveDotProcessorGenerator3D &ReductiveDotProcessorGenerator3D::operator=(
    ReductiveDotProcessorGenerator3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT ReductiveDotProcessorGenerator3D::appliesTo() const
{
    return functional->appliesTo();
}

void ReductiveDotProcessorGenerator3D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void ReductiveDotProcessorGenerator3D::getDimensionsX(std::vector<int> &dimensions) const
{
    functional->getDimensionsX(dimensions);
}

void ReductiveDotProcessorGenerator3D::getDimensionsT(std::vector<int> &dimensions) const
{
    functional->getDimensionsT(dimensions);
}

void ReductiveDotProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void ReductiveDotProcessorGenerator3D::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor3D *ReductiveDotProcessorGenerator3D::generate(
    std::vector<AtomicBlock3D *> atomicBlocks)
{
    return new ReductiveDotProcessor3D(functional, this->getDotList(), atomicBlocks);
}

ReductiveDotProcessorGenerator3D *ReductiveDotProcessorGenerator3D::clone() const
{
    return new ReductiveDotProcessorGenerator3D(*this);
}

BlockStatistics const &ReductiveDotProcessorGenerator3D::getStatistics() const
{
    return functional->getStatistics();
}

BlockStatistics &ReductiveDotProcessorGenerator3D::getStatistics()
{
    return functional->getStatistics();
}

ReductiveDotProcessingFunctional3D const &ReductiveDotProcessorGenerator3D::getFunctional() const
{
    return *functional;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D ************************* */

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void BoundedReductiveBoxProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void BoundedReductiveBoxProcessingFunctional3D::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void BoundedReductiveBoxProcessingFunctional3D::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.clear();
}

void BoundedReductiveBoxProcessingFunctional3D::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

ReductiveBoxProcessingFunctional3D *BoundedReductiveBoxProcessingFunctional3D::getBulkProcessor()
    const
{
    return new BulkWrapperFunctional(this->clone());
}

ReductiveBoxProcessingFunctional3D *BoundedReductiveBoxProcessingFunctional3D::getPlaneProcessor(
    int direction, int orientation) const
{
    return new PlaneWrapperFunctional(this->clone(), direction, orientation);
}

ReductiveBoxProcessingFunctional3D *BoundedReductiveBoxProcessingFunctional3D::getEdgeProcessor(
    int plane, int normal1, int normal2) const
{
    return new EdgeWrapperFunctional(this->clone(), plane, normal1, normal2);
}

ReductiveBoxProcessingFunctional3D *BoundedReductiveBoxProcessingFunctional3D::getCornerProcessor(
    int normalX, int normalY, int normalZ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY, normalZ);
}

BlockStatistics const &BoundedReductiveBoxProcessingFunctional3D::getStatistics() const
{
    return statistics;
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional3D::getStatistics()
{
    return statistics;
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::BulkWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_) :
    boundedFunctional(boundedFunctional_)
{ }

BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::BulkWrapperFunctional(
    BulkWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone())
{ }

BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::~BulkWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::operator=(
        BulkWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional3D::BulkWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::PlaneWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int direction_,
    int orientation_) :
    boundedFunctional(boundedFunctional_), direction(direction_), orientation(orientation_)
{ }

BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::PlaneWrapperFunctional(
    PlaneWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::~PlaneWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::operator=(
        PlaneWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processPlaneGeneric(direction, orientation, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::clone() const
{
    return new PlaneWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional3D::PlaneWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int plane_, int normal1_,
    int normal2_) :
    boundedFunctional(boundedFunctional_), plane(plane_), normal1(normal1_), normal2(normal2_)
{ }

BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    EdgeWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    plane(rhs.plane),
    normal1(rhs.normal1),
    normal2(rhs.normal2)
{ }

BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::~EdgeWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::operator=(
        EdgeWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    plane = rhs.plane;
    normal1 = rhs.normal1;
    normal2 = rhs.normal2;
    return *this;
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processEdgeGeneric(plane, normal1, normal2, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional3D::EdgeWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

/* *************** Class BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional ** */

BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::CornerWrapperFunctional(
    BoundedReductiveBoxProcessingFunctional3D *boundedFunctional_, int normalX_, int normalY_,
    int normalZ_) :
    boundedFunctional(boundedFunctional_), normalX(normalX_), normalY(normalY_), normalZ(normalZ_)
{ }

BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::CornerWrapperFunctional(
    CornerWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    normalX(rhs.normalX),
    normalY(rhs.normalY),
    normalZ(rhs.normalZ)
{ }

BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::~CornerWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional &
    BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::operator=(
        CornerWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    normalZ = rhs.normalZ;
    return *this;
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processCornerGeneric(normalX, normalY, normalZ, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::appliesTo()
    const
{
    return boundedFunctional->appliesTo();
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getDimensionsX(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsX(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getDimensionsT(
    std::vector<int> &dimensions) const
{
    boundedFunctional->getDimensionsT(dimensions);
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional *
    BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

BlockStatistics const &
    BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getStatistics() const
{
    return boundedFunctional->getStatistics();
}

BlockStatistics &BoundedReductiveBoxProcessingFunctional3D::CornerWrapperFunctional::getStatistics()
{
    return boundedFunctional->getStatistics();
}

void BoundedReductiveBoxProcessingFunctional3D::getGenerators(
    Box3D const &fullDomain, plint boundaryWidth,
    std::vector<ReductiveBoxProcessorGenerator3D *> &generators)
{
    generators.resize(27);
    BlockSurface3D surf(fullDomain, boundaryWidth);
    generators[0] = new ReductiveBoxProcessorGenerator3D(getBulkProcessor(), surf.bulk());

    generators[1] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(0, -1), surf.surface0N());
    generators[2] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(0, +1), surf.surface0P());
    generators[3] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(1, -1), surf.surface1N());
    generators[4] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(1, +1), surf.surface1P());
    generators[5] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(2, -1), surf.surface2N());
    generators[6] =
        new ReductiveBoxProcessorGenerator3D(getPlaneProcessor(2, +1), surf.surface2P());

    generators[7] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(0, -1, -1), surf.edge0NN());
    generators[8] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(0, -1, 1), surf.edge0NP());
    generators[9] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(0, 1, -1), surf.edge0PN());
    generators[10] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(0, 1, 1), surf.edge0PP());
    generators[11] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(1, -1, -1), surf.edge1NN());
    generators[12] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(1, -1, 1), surf.edge1NP());
    generators[13] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(1, 1, -1), surf.edge1PN());
    generators[14] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(1, 1, 1), surf.edge1PP());
    generators[15] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(2, -1, -1), surf.edge2NN());
    generators[16] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(2, -1, 1), surf.edge2NP());
    generators[17] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(2, 1, -1), surf.edge2PN());
    generators[18] =
        new ReductiveBoxProcessorGenerator3D(getEdgeProcessor(2, 1, 1), surf.edge2PP());

    generators[19] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(-1, -1, -1), surf.cornerNNN());
    generators[20] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(-1, -1, 1), surf.cornerNNP());
    generators[21] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(-1, 1, -1), surf.cornerNPN());
    generators[22] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(-1, 1, 1), surf.cornerNPP());
    generators[23] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(1, -1, -1), surf.cornerPNN());
    generators[24] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(1, -1, 1), surf.cornerPNP());
    generators[25] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(1, 1, -1), surf.cornerPPN());
    generators[26] =
        new ReductiveBoxProcessorGenerator3D(getCornerProcessor(1, 1, 1), surf.cornerPPP());
}

}  // namespace plb
