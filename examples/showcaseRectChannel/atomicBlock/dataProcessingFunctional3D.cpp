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
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** Class BoxProcessingFunctional3D ************************* */

BoxProcessingFunctional3D::BoxProcessingFunctional3D() : dxScale(0), dtScale(0) { }

/** Operations are not executed on envelope by default. **/
BlockDomain::DomainT BoxProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

void BoxProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void BoxProcessingFunctional3D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int BoxProcessingFunctional3D::getDxScale() const
{
    return dxScale;
}

int BoxProcessingFunctional3D::getDtScale() const
{
    return dtScale;
}

void BoxProcessingFunctional3D::serialize(std::string &data) const
{
    std::string newData = util::val2str(dxScale, dtScale);
    data += newData + " ";
}

void BoxProcessingFunctional3D::unserialize(std::string &data)
{
    data = util::consume(data, dxScale, dtScale);
}

/** Return a default value of -1. This is to help transition for legacy
 *  code that has not yet implemented this method.
 ***/
int BoxProcessingFunctional3D::getStaticId() const
{
    return -1;
}

void BoxProcessingFunctional3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class BoxProcessor3D ************************************ */

BoxProcessor3D::BoxProcessor3D(
    BoxProcessingFunctional3D *functional_, Box3D domain_,
    std::vector<AtomicBlock3D *> atomicBlocks_) :
    functional(functional_), domain(domain_), atomicBlocks(atomicBlocks_)
{ }

BoxProcessor3D::BoxProcessor3D(BoxProcessor3D const &rhs) :
    functional(rhs.functional->clone()), domain(rhs.domain), atomicBlocks(rhs.atomicBlocks)
{ }

BoxProcessor3D &BoxProcessor3D::operator=(BoxProcessor3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    domain = rhs.domain;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

BoxProcessor3D::~BoxProcessor3D()
{
    delete functional;
}

Box3D BoxProcessor3D::getDomain() const
{
    return domain;
}

void BoxProcessor3D::process()
{
    functional->processGenericBlocks(domain, atomicBlocks);
}

BoxProcessor3D *BoxProcessor3D::clone() const
{
    return new BoxProcessor3D(*this);
}

int BoxProcessor3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class MultiBoxProcessor3D ************************************ */

MultiBoxProcessor3D::MultiBoxProcessor3D(
    BoxProcessingFunctional3D *functional_, std::vector<Box3D> const &domains_,
    std::vector<AtomicBlock3D *> atomicBlocks_) :
    functional(functional_), domains(domains_), atomicBlocks(atomicBlocks_)
{ }

MultiBoxProcessor3D::MultiBoxProcessor3D(MultiBoxProcessor3D const &rhs) :
    functional(rhs.functional->clone()), domains(rhs.domains), atomicBlocks(rhs.atomicBlocks)
{ }

MultiBoxProcessor3D &MultiBoxProcessor3D::operator=(MultiBoxProcessor3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    domains = rhs.domains;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

MultiBoxProcessor3D::~MultiBoxProcessor3D()
{
    delete functional;
}

std::vector<Box3D> const &MultiBoxProcessor3D::getDomains() const
{
    return domains;
}

void MultiBoxProcessor3D::process()
{
    for (pluint i = 0; i < domains.size(); ++i) {
        functional->processGenericBlocks(domains[i], atomicBlocks);
    }
}

MultiBoxProcessor3D *MultiBoxProcessor3D::clone() const
{
    return new MultiBoxProcessor3D(*this);
}

int MultiBoxProcessor3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class BoxProcessorGenerator3D *************************** */

BoxProcessorGenerator3D::BoxProcessorGenerator3D(
    BoxProcessingFunctional3D *functional_, Box3D domain) :
    BoxedDataProcessorGenerator3D(domain), functional(functional_)
{ }

BoxProcessorGenerator3D::~BoxProcessorGenerator3D()
{
    delete functional;
}

BoxProcessorGenerator3D::BoxProcessorGenerator3D(BoxProcessorGenerator3D const &rhs) :
    BoxedDataProcessorGenerator3D(rhs), functional(rhs.functional->clone())
{ }

BoxProcessorGenerator3D &BoxProcessorGenerator3D::operator=(BoxProcessorGenerator3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT BoxProcessorGenerator3D::appliesTo() const
{
    return functional->appliesTo();
}

void BoxProcessorGenerator3D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void BoxProcessorGenerator3D::setscale(int dxScale, int dtScale)
{
    functional->setscale(dxScale, dtScale);
}

void BoxProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void BoxProcessorGenerator3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor3D *BoxProcessorGenerator3D::generate(std::vector<AtomicBlock3D *> atomicBlocks) const
{
    return new BoxProcessor3D(functional->clone(), this->getDomain(), atomicBlocks);
}

BoxProcessorGenerator3D *BoxProcessorGenerator3D::clone() const
{
    return new BoxProcessorGenerator3D(*this);
}

void BoxProcessorGenerator3D::serialize(Box3D &domain, std::string &data) const
{
    BoxedDataProcessorGenerator3D::serialize(domain, data);
    functional->serialize(data);
}

int BoxProcessorGenerator3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class MultiBoxProcessorGenerator3D *************************** */

MultiBoxProcessorGenerator3D::MultiBoxProcessorGenerator3D(
    BoxProcessingFunctional3D *functional_, std::vector<Box3D> const &domains) :
    MultiBoxedDataProcessorGenerator3D(domains), functional(functional_)
{ }

MultiBoxProcessorGenerator3D::~MultiBoxProcessorGenerator3D()
{
    delete functional;
}

MultiBoxProcessorGenerator3D::MultiBoxProcessorGenerator3D(
    MultiBoxProcessorGenerator3D const &rhs) :
    MultiBoxedDataProcessorGenerator3D(rhs), functional(rhs.functional->clone())
{ }

MultiBoxProcessorGenerator3D &MultiBoxProcessorGenerator3D::operator=(
    MultiBoxProcessorGenerator3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT MultiBoxProcessorGenerator3D::appliesTo() const
{
    return functional->appliesTo();
}

void MultiBoxProcessorGenerator3D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void MultiBoxProcessorGenerator3D::setscale(int dxScale, int dtScale)
{
    functional->setscale(dxScale, dtScale);
}

void MultiBoxProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void MultiBoxProcessorGenerator3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor3D *MultiBoxProcessorGenerator3D::generate(
    std::vector<AtomicBlock3D *> atomicBlocks) const
{
    return new MultiBoxProcessor3D(functional->clone(), this->getDomains(), atomicBlocks);
}

MultiBoxProcessorGenerator3D *MultiBoxProcessorGenerator3D::clone() const
{
    return new MultiBoxProcessorGenerator3D(*this);
}

void MultiBoxProcessorGenerator3D::serialize(Box3D &domain, std::string &data) const
{
    MultiBoxedDataProcessorGenerator3D::serialize(domain, data);
    functional->serialize(data);
}

int MultiBoxProcessorGenerator3D::getStaticId() const
{
    return functional->getStaticId();
}

/* *************** Class DotProcessingFunctional3D ************************* */

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT DotProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void DotProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void DotProcessingFunctional3D::setscale(int dxScale, int dtScale) { }

/** The default assumption is conservative: all blocks have potentially been modified.
 */
void DotProcessingFunctional3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

/* *************** Class DotProcessor3D ************************************ */

DotProcessor3D::DotProcessor3D(
    DotProcessingFunctional3D *functional_, DotList3D const &dotList_,
    std::vector<AtomicBlock3D *> atomicBlocks_) :
    functional(functional_), dotList(dotList_), atomicBlocks(atomicBlocks_)
{ }

DotProcessor3D::DotProcessor3D(DotProcessor3D const &rhs) :
    functional(rhs.functional->clone()), dotList(rhs.dotList), atomicBlocks(rhs.atomicBlocks)
{ }

DotProcessor3D &DotProcessor3D::operator=(DotProcessor3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    dotList = rhs.dotList;
    atomicBlocks = rhs.atomicBlocks;
    return *this;
}

DotProcessor3D::~DotProcessor3D()
{
    delete functional;
}

void DotProcessor3D::process()
{
    functional->processGenericBlocks(dotList, atomicBlocks);
}

DotProcessor3D *DotProcessor3D::clone() const
{
    return new DotProcessor3D(*this);
}

DotList3D const &DotProcessor3D::getDotList() const
{
    return dotList;
}

/* *************** Class DotProcessorGenerator3D *************************** */

DotProcessorGenerator3D::DotProcessorGenerator3D(
    DotProcessingFunctional3D *functional_, DotList3D const &dotList) :
    DottedDataProcessorGenerator3D(dotList), functional(functional_)
{ }

DotProcessorGenerator3D::~DotProcessorGenerator3D()
{
    delete functional;
}

DotProcessorGenerator3D::DotProcessorGenerator3D(DotProcessorGenerator3D const &rhs) :
    DottedDataProcessorGenerator3D(rhs), functional(rhs.functional->clone())
{ }

DotProcessorGenerator3D &DotProcessorGenerator3D::operator=(DotProcessorGenerator3D const &rhs)
{
    delete functional;
    functional = rhs.functional->clone();
    return *this;
}

BlockDomain::DomainT DotProcessorGenerator3D::appliesTo() const
{
    return functional->appliesTo();
}

void DotProcessorGenerator3D::rescale(double dxScale, double dtScale)
{
    functional->rescale(dxScale, dtScale);
}

void DotProcessorGenerator3D::setscale(int dxScale, int dtScale)
{
    functional->setscale(dxScale, dtScale);
}

void DotProcessorGenerator3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    functional->getModificationPattern(isWritten);
}

void DotProcessorGenerator3D::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    functional->getTypeOfModification(modified);
}

DataProcessor3D *DotProcessorGenerator3D::generate(std::vector<AtomicBlock3D *> atomicBlocks) const
{
    return new DotProcessor3D(functional->clone(), this->getDotList(), atomicBlocks);
}

DotProcessorGenerator3D *DotProcessorGenerator3D::clone() const
{
    return new DotProcessorGenerator3D(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D ************************* */

BoundedBoxProcessingFunctional3D::BoundedBoxProcessingFunctional3D() : dxScale(0), dtScale(0) { }

/** Operation is not applied to envelope by default. **/
BlockDomain::DomainT BoundedBoxProcessingFunctional3D::appliesTo() const
{
    return BlockDomain::bulk;
}

/** No rescaling is done by default. **/
void BoundedBoxProcessingFunctional3D::rescale(double dxScale, double dtScale) { }

void BoundedBoxProcessingFunctional3D::setscale(int dxScale_, int dtScale_)
{
    dxScale = dxScale_;
    dtScale = dtScale_;
}

int BoundedBoxProcessingFunctional3D::getDxScale() const
{
    return dxScale;
}

int BoundedBoxProcessingFunctional3D::getDtScale() const
{
    return dtScale;
}

void BoundedBoxProcessingFunctional3D::getModificationPattern(std::vector<bool> &isWritten) const
{
    std::vector<modif::ModifT> modified(isWritten.size());
    getTypeOfModification(modified);
    PLB_ASSERT(modified.size() == isWritten.size());
    for (pluint iBlock = 0; iBlock < isWritten.size(); ++iBlock) {
        isWritten[iBlock] = modified[iBlock] == modif::nothing ? false : true;
    }
}

BoxProcessingFunctional3D *BoundedBoxProcessingFunctional3D::getBulkProcessor() const
{
    return new BulkWrapperFunctional(this->clone());
}

BoxProcessingFunctional3D *BoundedBoxProcessingFunctional3D::getPlaneProcessor(
    int direction, int orientation) const
{
    return new PlaneWrapperFunctional(this->clone(), direction, orientation);
}

BoxProcessingFunctional3D *BoundedBoxProcessingFunctional3D::getEdgeProcessor(
    int plane, int normal1, int normal2) const
{
    return new EdgeWrapperFunctional(this->clone(), plane, normal1, normal2);
}

BoxProcessingFunctional3D *BoundedBoxProcessingFunctional3D::getCornerProcessor(
    int normalX, int normalY, int normalZ) const
{
    return new CornerWrapperFunctional(this->clone(), normalX, normalY, normalZ);
}

/* *************** Class BoundedBoxProcessingFunctional3D::BulkWrapperFunctional ** */

BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::BulkWrapperFunctional(
    BoundedBoxProcessingFunctional3D *boundedFunctional_) :
    boundedFunctional(boundedFunctional_)
{ }

BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::BulkWrapperFunctional(
    BulkWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone())
{ }

BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::~BulkWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional3D::BulkWrapperFunctional &
    BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::operator=(
        BulkWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    return *this;
}

void BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processBulkGeneric(domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional3D::BulkWrapperFunctional *
    BoundedBoxProcessingFunctional3D::BulkWrapperFunctional::clone() const
{
    return new BulkWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional ** */

BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::PlaneWrapperFunctional(
    BoundedBoxProcessingFunctional3D *boundedFunctional_, int direction_, int orientation_) :
    boundedFunctional(boundedFunctional_), direction(direction_), orientation(orientation_)
{ }

BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::PlaneWrapperFunctional(
    PlaneWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    direction(rhs.direction),
    orientation(rhs.orientation)
{ }

BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::~PlaneWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional &
    BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::operator=(
        PlaneWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    direction = rhs.direction;
    orientation = rhs.orientation;
    return *this;
}

void BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processPlaneGeneric(direction, orientation, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional *
    BoundedBoxProcessingFunctional3D::PlaneWrapperFunctional::clone() const
{
    return new PlaneWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional ** */

BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    BoundedBoxProcessingFunctional3D *boundedFunctional_, int plane_, int normal1_, int normal2_) :
    boundedFunctional(boundedFunctional_), plane(plane_), normal1(normal1_), normal2(normal2_)
{ }

BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::EdgeWrapperFunctional(
    EdgeWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    plane(rhs.plane),
    normal1(rhs.normal1),
    normal2(rhs.normal2)
{ }

BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::~EdgeWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional &
    BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::operator=(
        EdgeWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    plane = rhs.plane;
    normal1 = rhs.normal1;
    normal2 = rhs.normal2;
    return *this;
}

void BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processEdgeGeneric(plane, normal1, normal2, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional *
    BoundedBoxProcessingFunctional3D::EdgeWrapperFunctional::clone() const
{
    return new EdgeWrapperFunctional(*this);
}

/* *************** Class BoundedBoxProcessingFunctional3D::CornerWrapperFunctional ** */

BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::CornerWrapperFunctional(
    BoundedBoxProcessingFunctional3D *boundedFunctional_, int normalX_, int normalY_,
    int normalZ_) :
    boundedFunctional(boundedFunctional_), normalX(normalX_), normalY(normalY_), normalZ(normalZ_)
{ }

BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::CornerWrapperFunctional(
    CornerWrapperFunctional const &rhs) :
    boundedFunctional(rhs.boundedFunctional->clone()),
    normalX(rhs.normalX),
    normalY(rhs.normalY),
    normalZ(rhs.normalZ)
{ }

BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::~CornerWrapperFunctional()
{
    delete boundedFunctional;
}

BoundedBoxProcessingFunctional3D::CornerWrapperFunctional &
    BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::operator=(
        CornerWrapperFunctional const &rhs)
{
    delete boundedFunctional;
    boundedFunctional = rhs.boundedFunctional->clone();
    normalX = rhs.normalX;
    normalY = rhs.normalY;
    normalZ = rhs.normalZ;
    return *this;
}

void BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    boundedFunctional->processCornerGeneric(normalX, normalY, normalZ, domain, atomicBlocks);
}

BlockDomain::DomainT BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::appliesTo() const
{
    return boundedFunctional->appliesTo();
}

void BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::rescale(
    double dxScale, double dtScale)
{
    boundedFunctional->rescale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::setscale(int dxScale, int dtScale)
{
    boundedFunctional->setscale(dxScale, dtScale);
}

void BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    boundedFunctional->getModificationPattern(isWritten);
}

void BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    boundedFunctional->getTypeOfModification(modified);
}

BoundedBoxProcessingFunctional3D::CornerWrapperFunctional *
    BoundedBoxProcessingFunctional3D::CornerWrapperFunctional::clone() const
{
    return new CornerWrapperFunctional(*this);
}

void BoundedBoxProcessingFunctional3D::getGenerators(
    Box3D const &fullDomain, plint boundaryWidth,
    std::vector<BoxProcessorGenerator3D *> &generators)
{
    generators.resize(27);
    BlockSurface3D surf(fullDomain, boundaryWidth);
    generators[0] = new BoxProcessorGenerator3D(getBulkProcessor(), surf.bulk());

    generators[1] = new BoxProcessorGenerator3D(getPlaneProcessor(0, -1), surf.surface0N());
    generators[2] = new BoxProcessorGenerator3D(getPlaneProcessor(0, +1), surf.surface0P());
    generators[3] = new BoxProcessorGenerator3D(getPlaneProcessor(1, -1), surf.surface1N());
    generators[4] = new BoxProcessorGenerator3D(getPlaneProcessor(1, +1), surf.surface1P());
    generators[5] = new BoxProcessorGenerator3D(getPlaneProcessor(2, -1), surf.surface2N());
    generators[6] = new BoxProcessorGenerator3D(getPlaneProcessor(2, +1), surf.surface2P());

    generators[7] = new BoxProcessorGenerator3D(getEdgeProcessor(0, -1, -1), surf.edge0NN());
    generators[8] = new BoxProcessorGenerator3D(getEdgeProcessor(0, -1, 1), surf.edge0NP());
    generators[9] = new BoxProcessorGenerator3D(getEdgeProcessor(0, 1, -1), surf.edge0PN());
    generators[10] = new BoxProcessorGenerator3D(getEdgeProcessor(0, 1, 1), surf.edge0PP());
    generators[11] = new BoxProcessorGenerator3D(getEdgeProcessor(1, -1, -1), surf.edge1NN());
    generators[12] = new BoxProcessorGenerator3D(getEdgeProcessor(1, -1, 1), surf.edge1NP());
    generators[13] = new BoxProcessorGenerator3D(getEdgeProcessor(1, 1, -1), surf.edge1PN());
    generators[14] = new BoxProcessorGenerator3D(getEdgeProcessor(1, 1, 1), surf.edge1PP());
    generators[15] = new BoxProcessorGenerator3D(getEdgeProcessor(2, -1, -1), surf.edge2NN());
    generators[16] = new BoxProcessorGenerator3D(getEdgeProcessor(2, -1, 1), surf.edge2NP());
    generators[17] = new BoxProcessorGenerator3D(getEdgeProcessor(2, 1, -1), surf.edge2PN());
    generators[18] = new BoxProcessorGenerator3D(getEdgeProcessor(2, 1, 1), surf.edge2PP());

    generators[19] = new BoxProcessorGenerator3D(getCornerProcessor(-1, -1, -1), surf.cornerNNN());
    generators[20] = new BoxProcessorGenerator3D(getCornerProcessor(-1, -1, 1), surf.cornerNNP());
    generators[21] = new BoxProcessorGenerator3D(getCornerProcessor(-1, 1, -1), surf.cornerNPN());
    generators[22] = new BoxProcessorGenerator3D(getCornerProcessor(-1, 1, 1), surf.cornerNPP());
    generators[23] = new BoxProcessorGenerator3D(getCornerProcessor(1, -1, -1), surf.cornerPNN());
    generators[24] = new BoxProcessorGenerator3D(getCornerProcessor(1, -1, 1), surf.cornerPNP());
    generators[25] = new BoxProcessorGenerator3D(getCornerProcessor(1, 1, -1), surf.cornerPPN());
    generators[26] = new BoxProcessorGenerator3D(getCornerProcessor(1, 1, 1), surf.cornerPPP());
}

}  // namespace plb
