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
 * 2D couplings -- header file.
 */
#ifndef COUPLING_2D_H
#define COUPLING_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

class Coupling2D {
public:
    Coupling2D(
        DataProcessorGenerator2D *generator_,
        std::vector<plint> const &multiBlocks_,  // Local ID (0, 1, 2, ...) inside Actions2D of
                                                 // multi-blocks to be coupled.
        std::vector<id_t> &allMultiBlocks);      // Global ID of all multi-blocks inside Actions2D.
    ~Coupling2D();
    Coupling2D(Coupling2D const &rhs);
    Coupling2D &operator=(Coupling2D const &rhs);
    void swap(Coupling2D &rhs);
    Coupling2D *clone() const;
    void execute();
    void generate(std::vector<id_t> &allMultiBlocks);

private:
    void clearDataProcessors();

private:
    DataProcessorGenerator2D *generator;
    std::vector<plint> multiBlocks;
    std::vector<DataProcessor2D *> dataProcessors;
};

struct Action2D {
    virtual ~Action2D() { }
    virtual void execute(std::vector<id_t> &allMultiBlocks) = 0;
    virtual void regenerate(std::vector<id_t> &allMultiBlocks) = 0;
    virtual Action2D *clone() const = 0;
};

class CouplingAction2D : public Action2D {
public:
    CouplingAction2D(Coupling2D *coupling_);
    CouplingAction2D(CouplingAction2D const &rhs);
    CouplingAction2D &operator=(CouplingAction2D const &rhs);
    void swap(CouplingAction2D &rhs);
    virtual ~CouplingAction2D();
    virtual CouplingAction2D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    Coupling2D *coupling;
};

template <typename T, template <typename U> class Descriptor>
class FullDomainCollideAndStreamAction2D : public Action2D {
public:
    FullDomainCollideAndStreamAction2D(plint blockId_);
    virtual FullDomainCollideAndStreamAction2D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

template <typename T, template <typename U> class Descriptor>
class CollideAndStreamAction2D : public Action2D {
public:
    CollideAndStreamAction2D(plint blockId_, Box2D domain_);
    virtual CollideAndStreamAction2D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    Box2D domain;
};

template <typename T, template <typename U> class Descriptor>
class FullDomainStreamAction2D : public Action2D {
public:
    FullDomainStreamAction2D(plint blockId_);
    virtual FullDomainStreamAction2D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

template <typename T, template <typename U> class Descriptor>
class StreamAction2D : public Action2D {
public:
    StreamAction2D(plint blockId_, Box2D domain_);
    virtual StreamAction2D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    Box2D domain;
};

template <typename T, template <typename U> class Descriptor>
class IncrementTimeAction2D : public Action2D {
public:
    IncrementTimeAction2D(plint blockId_);
    virtual IncrementTimeAction2D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

class CommunicateAction2D : public Action2D {
public:
    CommunicateAction2D(plint blockId_, modif::ModifT whichData_);
    virtual CommunicateAction2D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    modif::ModifT whichData;
};

class ExecuteInternalProcAction2D : public Action2D {
public:
    ExecuteInternalProcAction2D(plint blockId_, plint level_);
    virtual ExecuteInternalProcAction2D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId, level;
};

class EvaluateStatsAction2D : public Action2D {
public:
    EvaluateStatsAction2D(plint blockId_);
    virtual EvaluateStatsAction2D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    virtual void regenerate(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

class Actions2D {
public:
    Actions2D();
    ~Actions2D();
    Actions2D(Actions2D const &rhs);
    Actions2D &operator=(Actions2D const &rhs);
    void swap(Actions2D &rhs);
    Actions2D *clone() const;
    plint addBlock(MultiBlock2D &block);
    void replaceBlock(plint id, MultiBlock2D &block);
    plint addProcessor(
        BoxProcessingFunctional2D *functional, std::vector<plint> blockNums, Box2D domain);
    plint addProcessor(BoxProcessingFunctional2D *functional, plint blockNum1, Box2D domain);
    plint addProcessor(
        BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, Box2D domain);
    plint addProcessor(
        BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        Box2D domain);
    plint addProcessor(
        BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        plint blockNum4, Box2D domain);
    plint addProcessor(
        BoxProcessingFunctional2D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        plint blockNum4, plint blockNum5, Box2D domain);
    plint addInternalProcessors(plint blockNum, plint level);
    plint addCommunication(plint blockNum, modif::ModifT whichData);
    plint addEvaluateStats(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addCollideAndStream(plint blockNum, Box2D domain);
    template <typename T, template <typename U> class Descriptor>
    plint addCollideAndStream(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addStream(plint blockNum, Box2D domain);
    template <typename T, template <typename U> class Descriptor>
    plint addStream(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addIncrementTime(plint blockNum);
    void execute();
    void execute(plint actionId);
    void execute(plint actionFrom, plint actionTo);

private:
    std::vector<id_t> allMultiBlocks;
    std::vector<Action2D *> actions;
};

}  // namespace plb

#endif  // COUPLING_2D_H
