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
 * Interface for dataProcessing steps -- header file.
 */
#ifndef DATA_PROCESSOR_3D_H
#define DATA_PROCESSOR_3D_H

#include <algorithm>
#include <vector>

#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"

namespace plb {

// Forward declarations
class AtomicBlock3D;

/// DataProcessors are used to run extended operations on a lattice or data field
struct DataProcessor3D {
    virtual ~DataProcessor3D() { }
    /// Execute processing operation
    virtual void process() = 0;
    /// Clone Data Processor, on its dynamic type
    virtual DataProcessor3D *clone() const = 0;
    /// Extent of application area (0 for purely local operations)
    virtual plint extent() const;
    /// Extent of application area along a direction (0 or 1)
    virtual plint extent(int direction) const;
    /// Unique identifier for a given DataProcessor class. Produces the same ID as
    ///   the corresponding processor generator.
    virtual int getStaticId() const;
};

/// This is a factory class generating LatticeProcessors
/** The LatticeProcessorGenerator can be tailored (shifted/reduced) to
 *  a sublattice, after which the LatticeProcessor is generated. The
 *  LatticeProcessor itself is static, i.e. the coordinates of the
 *  sublattice to which it refers cannot be changed after construction.
 *  Instead, a new LatticeProcessor must be generated with the generator.
 */
struct DataProcessorGenerator3D {
    virtual ~DataProcessorGenerator3D();
    /// Shift the domain of application of this data processor.
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ) = 0;
    /// Multiply coordinates of the domain of application of this data processor.
    virtual void multiply(plint scale) = 0;
    /// Divide coordinates of the domain of application of this data processor.
    virtual void divide(plint scale) = 0;
    /// Extract a subdomain (in-place operation).
    /** \return True if original domain of application and subDomain intersect.
     */
    virtual bool extract(Box3D subDomain) = 0;
    /// Generate the data processor.
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks) const = 0;
    /// Clone DataProcessorGenerator, based on its dynamics type.
    virtual DataProcessorGenerator3D *clone() const = 0;
    /// Indicates whether data processor should be applied on envelope or not. Defaults to false.
    virtual BlockDomain::DomainT appliesTo() const;
    /// This function is obsolete, and has been replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    /// Specify the scale of the block on which the data processor is acting. Defaults to no
    /// rescaling.
    virtual void setscale(int dxScale_, int dtScale_);
    /// Tell which blocks are modified (written) when the processor is applied on them.
    /** This function is obsolete and should be replaced by getTypeOfModification().
     **/
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    /// Tell which blocks are modified and how by the processor. This method must
    /// be implemented in each data processor.
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    /// Unique identifier for a given DataProcessor class. Produces the same ID as
    ///   the corresponding data processor.
    virtual int getStaticId() const;
    /// Serialize content into a string-stream and overwrite data with result.
    virtual void serialize(Box3D &domain, std::string &data) const;
};

class BoxedDataProcessorGenerator3D : public DataProcessorGenerator3D {
public:
    BoxedDataProcessorGenerator3D(Box3D domain_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    Box3D getDomain() const;
    virtual void serialize(Box3D &domain, std::string &data) const;

private:
    Box3D domain;
};

class MultiBoxedDataProcessorGenerator3D : public DataProcessorGenerator3D {
public:
    MultiBoxedDataProcessorGenerator3D(std::vector<Box3D> const &domains_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    std::vector<Box3D> const &getDomains() const;
    virtual void serialize(Box3D &domain, std::string &data) const;

private:
    std::vector<Box3D> domains;
};

class DottedDataProcessorGenerator3D : public DataProcessorGenerator3D {
public:
    DottedDataProcessorGenerator3D(DotList3D const &dots_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    DotList3D const &getDotList() const;

private:
    DotList3D dots;
};

class ReductiveDataProcessorGenerator3D {
public:
    ReductiveDataProcessorGenerator3D();
    virtual ~ReductiveDataProcessorGenerator3D();
    /// Shift the domain of application of this data processor.
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ) = 0;
    /// Multiply coordinates of the domain of application of this data processor.
    virtual void multiply(plint scale) = 0;
    /// Divide coordinates of the domain of application of this data processor.
    virtual void divide(plint scale) = 0;
    /// Extract a subdomain (in-place operation).
    /** \return True if original domain of application and subDomain intersect.
     */
    virtual bool extract(Box3D subDomain) = 0;
    /// Generate the data processor.
    virtual DataProcessor3D *generate(std::vector<AtomicBlock3D *> atomicBlocks) = 0;
    /// Clone ReductiveDataProcessorGenerator, based on its dynamics type.
    virtual ReductiveDataProcessorGenerator3D *clone() const = 0;
    /// Const handle to statistics object.
    virtual BlockStatistics const &getStatistics() const = 0;
    /// Non-const handle to statistics object.
    virtual BlockStatistics &getStatistics() = 0;
    /// Indicates whether data processor should be applied on envelope or not. Defaults to false.
    virtual BlockDomain::DomainT appliesTo() const;
    /// This function is obsolete, and has been replaced by setscale.
    virtual void rescale(double dxScale, double dtScale);
    /// Tell which blocks are modified (written) when the processor is applied on them.
    /** This function is obsolete and should be replaced by getTypeOfModification().
     **/
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    /// Tell which blocks are modified and how by the processor. This method must be
    /// implemented in every data processor.
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const = 0;
    /// Specify the scale of the block on which the data processor is acting. Defaults to no
    /// rescaling.
    void setscale(int dxScale_, int dtScale_);
    /// Return the space scale of the subjacent block.
    int getDxScale() const;
    /// Return the time scale of the subjacent block.
    int getDtScale() const;
    /// Get the space dimensions of each reduced value.
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    /// Get the time dimensions of each reduced value.
    virtual void getDimensionsT(std::vector<int> &dimensions) const;
    virtual void serialize(Box3D &domain, std::string &data) const;

private:
    int dxScale, dtScale;
};

class BoxedReductiveDataProcessorGenerator3D : public ReductiveDataProcessorGenerator3D {
public:
    BoxedReductiveDataProcessorGenerator3D(Box3D domain_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    Box3D getDomain() const;
    virtual void serialize(Box3D &domain, std::string &data) const;

private:
    Box3D domain;
};

class MultiBoxedReductiveDataProcessorGenerator3D : public ReductiveDataProcessorGenerator3D {
public:
    MultiBoxedReductiveDataProcessorGenerator3D(std::vector<Box3D> const &domains_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    std::vector<Box3D> const &getDomains() const;
    virtual void serialize(Box3D &domain, std::string &data) const;

private:
    std::vector<Box3D> domains;
};

class DottedReductiveDataProcessorGenerator3D : public ReductiveDataProcessorGenerator3D {
public:
    DottedReductiveDataProcessorGenerator3D(DotList3D const &dots_);
    virtual void shift(plint deltaX, plint deltaY, plint deltaZ);
    virtual void multiply(plint scale);
    virtual void divide(plint scale);
    virtual bool extract(Box3D subDomain);
    DotList3D const &getDotList() const;

private:
    DotList3D dots;
};

}  // namespace plb

#endif  // DATA_PROCESSOR_3D_H
