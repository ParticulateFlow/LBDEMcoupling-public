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
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class BoxProcessingFunctional2D ************************* */

BoxProcessingFunctional2D::BoxProcessingFunctional2D() : dxScale(0), dtScale(0) { }

/** Operations are not executed on envelope by default. **/
BlockDomain::DomainT BoxProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

void BoxProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void BoxProcessingFunctional2D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int BoxProcessingFunctional2D::getDxScale() const
{
    return dxScale;
}

int BoxProcessingFunctional2D::getDtScale() const
{
    return dtScale;
}

void BoxProcessingFunctional2D::serialize(std::string &data) const
{
    std::string newData = util::val2str(dxScale, dtScale);
    data += newData + " ";
}

void BoxProcessingFunctional2D::unserialize(std::string &data)
{
    data = util::consume(data, dxScale, dtScale);
}

/** Return a default value of -1. This is to help transition for legacy
 *  code that has not yet implemented this method.
 ***/
int BoxProcessingFunctional2D::getStaticId() const
{
    return -1;
}

void BoxProcessingFunctional2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class BoxProcessor2D ************************************ */

BoxProcessor2D::BoxProcessor2D(
    BoxProcessingFunctional2D *functional_, Box2D domain_,
    std::vector<AtomicBlock2D *> atomicBlocks_) :
    functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

BoxProcessor2D::BoxProcessor2D(BoxProcessor2D const &rhs) :
    functional(rhs.functional->clone()), domain(rhs.domain), atomicBlocks(rhs.atomicBlocks)
{ }

BoxProcessor2D &BoxProcessor2D::operator=(BoxProcessor2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    domain = rhs.domain;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

BoxProcessor2D::~BoxProcessor2D()
{
    delete functional;
}

Box2D BoxProcessor2D::getDomain() const
{
    return domain;
}

void BoxProcessor2D::process()
{
    functional->processGenericBlocks(domain, atomicBlocks);
}

BoxProcessor2D *BoxProcessor2D::clone() const
{
    return new BoxProcessor2D(*this);
}

int BoxProcessor2D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class BoxProcessorGenerator2D *************************** */

BoxProcessorGenerator2D::BoxProcessorGenerator2D(
    BoxProcessingFunctional2D *functional_, Box2D domain) :
    BoxedDataProcessorGenerator2D(domain), functional(functional_)
{ }

BoxProcessorGenerator2D::~BoxProcessorGenerator2D()
{
    delete functional;
}

BoxProcessorGenerator2D::BoxProcessorGenerator2D(BoxProcessorGenerator2D const &rhs) :
    BoxedDataProcessorGenerator2D(rhs), functional(rhs.functional->clone())
{ }

BoxProcessorGenerator2D &BoxProcessorGenerator2D::operator=(BoxProcessorGenerator2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT BoxProcessorGenerator2D::appliesTo() const
{
    return functional->appliesTo();
}

void BoxProcessorGenerator2D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void BoxProcessorGenerator2D::setscale(int dxScale, int dtScale)
{
    functional->setscale(dxScale, dtScale);
}

void BoxProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void BoxProcessorGenerator2D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor2D *BoxProcessorGenerator2D::generate(std::vector<AtomicBlock2D *> atomicBlocks) const
{
    return new BoxProcessor2D(functional->clone(), this->getDomain(), atomicBlocks);
}

BoxProcessorGenerator2D *BoxProcessorGenerator2D::clone() const
{
    return new BoxProcessorGenerator2D(*this);
}

void BoxProcessorGenerator2D::serialize(Box2D &domain, std::string &data) const
{
    BoxedDataProcessorGenerator2D::serialize(domain, data);
    functional->serialize(data);
}

int BoxProcessorGenerator2D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class DotProcessor2D ************************************ */

DotProcessor2D::DotProcessor2D(
    DotProcessingFunctional2D *functional_, DotList2D const &dotList_,
    std::vector<AtomicBlock2D *> atomicBlocks_) :
    functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

DotProcessor2D::DotProcessor2D(DotProcessor2D const &rhs) :
    functional(rhs.functional->clone()), dotList(rhs.dotList), atomicBlocks(rhs.atomicBlocks)
{ }

DotProcessor2D &DotProcessor2D::operator=(DotProcessor2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    dotList = rhs.dotList;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

DotProcessor2D::~DotProcessor2D()
{
    delete functional;
}

void DotProcessor2D::process()
{
    functional->processGenericBlocks(dotList, atomicBlocks);
}

DotProcessor2D *DotProcessor2D::clone() const
{
    return new DotProcessor2D(*this);
}

DotList2D const &DotProcessor2D::getDotList() const
{
    return dotList;
}

/* *************** Class DotProcessorGenerator2D *************************** */

DotProcessorGenerator2D::DotProcessorGenerator2D(
    DotProcessingFunctional2D *functional_, DotList2D const &dotList) :
    DottedDataProcessorGenerator2D(dotList), functional(functional_)
{ }

DotProcessorGenerator2D::~DotProcessorGenerator2D()
{
    delete functional;
}

DotProcessorGenerator2D::DotProcessorGenerator2D(DotProcessorGenerator2D const &rhs) :
    DottedDataProcessorGenerator2D(rhs), functional(rhs.functional->clone())
{ }

DotProcessorGenerator2D &DotProcessorGenerator2D::operator=(DotProcessorGenerator2D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT DotProcessorGenerator2D::appliesTo() const
{
    return functional->appliesTo();
}

void DotProcessorGenerator2D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void DotProcessorGenerator2D::setscale(int dxScale, int dtScale)
{
    functional->setscale(dxScale, dtScale);
}

void DotProcessorGenerator2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void DotProcessorGenerator2D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor2D *DotProcessorGenerator2D::generate(std::vector<AtomicBlock2D *> atomicBlocks) const
{
    return new DotProcessor2D(functional->clone(), this->getDotList(), atomicBlocks);
}

DotProcessorGenerator2D *DotProcessorGenerator2D::clone() const
{
    return new DotProcessorGenerator2D(*this);
}

/* *************** Class DotProcessingFunctional2D ************************* */

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT DotProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void DotProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void DotProcessingFunctional2D::setscale(int dxScale, int dtScale) { }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
void DotProcessingFunctional2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class BoundedBoxProcessingFunctional2D ************************* */

BoundedBoxProcessingFunctional2D::BoundedBoxProcessingFunctional2D() : dxScale(0), dtScale(0) { }

/** Operation is not executed on envelope by default. **/
BlockDomain::DomainT BoundedBoxProcessingFunctional2D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void BoundedBoxProcessingFunctional2D::rescale(double dxScale, double dtScale) { }

void BoundedBoxProcessingFunctional2D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int BoundedBoxProcessingFunctional2D::getDxScale() const
{
    return dxScale;
}

int BoundedBoxProcessingFunctional2D::getDtScale() const
{
    return dtScale;
}

void BoundedBoxProcessingFunctional2D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

BoxProcessingFunctional2D *BoundedBoxProcessingFunctional2D::getBulkProcessor() const
{
    return new BulkWrapperFunctional(this->clone());
}

BoxProcessingFunctional2D *BoundedBoxProcessingFunctional2D::getEdgeProcessor(
    int direction, int orientation) const
{
    return new EdgeWrapperFunctional(this->clone(), direction, orientation);
}

BoxProcessingFunctional2D *BoundedBoxProcessingFunctional2D::getCornerProcessor(
    int normalX, int normalY) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY);
}

/* *************** Class BoundedBoxProcessingFunctional2D::BulkWrapperFunctional ** */

BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::BulkWrapperFunctional(
    BoundedBoxProcessingFunctional2D *boundedFunctional_) :
    boundedFunctional(boundedFunctional_)
{ }

BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::BulkWrapperFunctional(
    BulkWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone())
{ }

BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::~BulkWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional2D::BulkWrapperFunctional &
    BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::operator=(
        BulkWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

void BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional2D::BulkWrapperFunctional *
    BoundedBoxProcessingFunctional2D::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional ** */

BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    BoundedBoxProcessingFunctional2D *boundedFunctional_, int direction_, int orientation_) :
    boundedFunctional(boundedFunctional_), direction(direction_), orientation(orientation_)
{ }

BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    EdgeWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::~EdgeWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional &
    BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::operator=(
        EdgeWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

void BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processEdgeGeneric(direction, orientation, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional *
    BoundedBoxProcessingFunctional2D::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional2D::CornerWrapperFunctional ** */

BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::CornerWrapperFunctional(
    BoundedBoxProcessingFunctional2D *boundedFunctional_, int normalX_, int normalY_) :
    boundedFunctional(boundedFunctional_), normalX(normalX_), normalY(normalY_)
{ }

BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::CornerWrapperFunctional(
    CornerWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()), normalX(rhs.normalX), normalY(rhs.normalY)
{ }

BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::~CornerWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional2D::CornerWrapperFunctional &
    BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::operator=(
        CornerWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    return *this;
}

void BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> atomicBlocks)
{
    boundedFunctional->processCornerGeneric(normalX, normalY, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional2D::CornerWrapperFunctional *
    BoundedBoxProcessingFunctional2D::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

void BoundedBoxProcessingFunctional2D::getGenerators(
    Box2D const &fullDomain, plint boundaryWidth,
    std::vector<BoxProcessorGenerator2D *> &generators)
{
    generators.resize(9);
    BlockSurface2D surf(fullDomain, boundaryWidth);
    generators[0] = new BoxProcessorGenerator2D(getBulkProcessor(), surf.bulk());

    generators[1] = new BoxProcessorGenerator2D(getEdgeProcessor(0, -1), surf.edge0N());
    generators[2] = new BoxProcessorGenerator2D(getEdgeProcessor(0, +1), surf.edge0P());
    generators[3] = new BoxProcessorGenerator2D(getEdgeProcessor(1, -1), surf.edge1N());
    generators[4] = new BoxProcessorGenerator2D(getEdgeProcessor(1, +1), surf.edge1P());

    generators[5] = new BoxProcessorGenerator2D(getCornerProcessor(-1, -1), surf.cornerNN());
    generators[6] = new BoxProcessorGenerator2D(getCornerProcessor(+1, -1), surf.cornerPN());
    generators[7] = new BoxProcessorGenerator2D(getCornerProcessor(-1, +1), surf.cornerNP());
    generators[8] = new BoxProcessorGenerator2D(getCornerProcessor(+1, +1), surf.cornerPP());
}

}  // namespace plb
