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
 * 3D couplings -- header file.
 */
#ifndef COUPLING_3D_H
#define COUPLING_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

class Actions3D;

class Coupling3D {
public:
    Coupling3D(
        DataProcessorGenerator3D *generator_,
        std::vector<plint> const &multiBlocks_,  // Local ID (0, 1, 2, ...) inside Actions3D of
                                                 // multi-blocks to be coupled.
        std::vector<id_t> &allMultiBlocks);      // Global ID of all multi-blocks inside Actions3D.
    ~Coupling3D();
    Coupling3D(Coupling3D const &rhs);
    Coupling3D &operator=(Coupling3D const &rhs);
    void swap(Coupling3D &rhs);
    Coupling3D *clone() const;
    void execute();
    void generate();

private:
    void clearDataProcessors();

private:
    DataProcessorGenerator3D *generator;
    std::vector<id_t> multiBlocks;
    std::vector<DataProcessor3D *> dataProcessors;
};

class ReductiveCoupling3D {
public:
    ReductiveCoupling3D(
        ReductiveDataProcessorGenerator3D *generator_,
        std::vector<plint> const &multiBlocks_,  // Local ID (0, 1, 2, ...) inside Actions3D of
                                                 // multi-blocks to be coupled.
        std::vector<id_t> &allMultiBlocks);      // Global ID of all multi-blocks inside Actions3D.
    ~ReductiveCoupling3D();
    ReductiveCoupling3D(ReductiveCoupling3D const &rhs);
    ReductiveCoupling3D &operator=(ReductiveCoupling3D const &rhs);
    void swap(ReductiveCoupling3D &rhs);
    ReductiveCoupling3D *clone() const;
    void execute();
    void generate();
    BlockStatistics const &statistics() const;

private:
    void clearData();

private:
    ReductiveDataProcessorGenerator3D *generator;
    std::vector<id_t> multiBlocks;
    CombinedStatistics *combinedStatistics;
    std::vector<ReductiveDataProcessorGenerator3D *> retainedGenerators;
    std::vector<BlockStatistics const *> individualStatistics;
    std::vector<DataProcessor3D *> dataProcessors;
};

struct Action3D {
    virtual ~Action3D() { }
    virtual void execute(std::vector<id_t> &allMultiBlocks) = 0;
    virtual Action3D *clone() const = 0;
};

class CouplingAction3D : public Action3D {
public:
    CouplingAction3D(Coupling3D *coupling_);
    CouplingAction3D(CouplingAction3D const &rhs);
    CouplingAction3D &operator=(CouplingAction3D const &rhs);
    void swap(CouplingAction3D &rhs);
    virtual ~CouplingAction3D();
    virtual CouplingAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    Coupling3D *coupling;
};

class ReductiveCouplingAction3D : public Action3D {
public:
    ReductiveCouplingAction3D(ReductiveCoupling3D *coupling_);
    ReductiveCouplingAction3D(ReductiveCouplingAction3D const &rhs);
    ReductiveCouplingAction3D &operator=(ReductiveCouplingAction3D const &rhs);
    void swap(ReductiveCouplingAction3D &rhs);
    virtual ~ReductiveCouplingAction3D();
    virtual ReductiveCouplingAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    BlockStatistics const &statistics() const;

private:
    ReductiveCoupling3D *coupling;
};

template <typename T, template <typename U> class Descriptor>
class FullDomainCollideAndStreamAction3D : public Action3D {
public:
    FullDomainCollideAndStreamAction3D(plint blockId_);
    virtual FullDomainCollideAndStreamAction3D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

template <typename T, template <typename U> class Descriptor>
class CollideAndStreamAction3D : public Action3D {
public:
    CollideAndStreamAction3D(plint blockId_, Box3D domain_);
    virtual CollideAndStreamAction3D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    Box3D domain;
};

template <typename T, template <typename U> class Descriptor>
class FullDomainStreamAction3D : public Action3D {
public:
    FullDomainStreamAction3D(plint blockId_);
    virtual FullDomainStreamAction3D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

template <typename T, template <typename U> class Descriptor>
class StreamAction3D : public Action3D {
public:
    StreamAction3D(plint blockId_, Box3D domain_);
    virtual StreamAction3D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    Box3D domain;
};

template <typename T, template <typename U> class Descriptor>
class IncrementTimeAction3D : public Action3D {
public:
    IncrementTimeAction3D(plint blockId_);
    virtual IncrementTimeAction3D<T, Descriptor> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

class CommunicateAction3D : public Action3D {
public:
    CommunicateAction3D(plint blockId_, modif::ModifT whichData_);
    virtual CommunicateAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
    modif::ModifT whichData;
};

class ExecuteInternalProcAction3D : public Action3D {
public:
    ExecuteInternalProcAction3D(plint blockId_, plint level_);
    virtual ExecuteInternalProcAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId, level;
};

class ExecuteAllInternalProcAction3D : public Action3D {
public:
    ExecuteAllInternalProcAction3D(plint blockId_);
    virtual ExecuteAllInternalProcAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

class EvaluateStatsAction3D : public Action3D {
public:
    EvaluateStatsAction3D(plint blockId_);
    virtual EvaluateStatsAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    plint blockId;
};

class SubActionsAction3D : public Action3D {
public:
    SubActionsAction3D(Actions3D *actions_);
    SubActionsAction3D(SubActionsAction3D const &rhs);
    SubActionsAction3D &operator=(SubActionsAction3D const &rhs);
    void swap(SubActionsAction3D &rhs);
    virtual ~SubActionsAction3D();
    virtual SubActionsAction3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
    Actions3D &get();

private:
    Actions3D *actions;
};

class Actions3D {
public:
    Actions3D();
    ~Actions3D();
    Actions3D(Actions3D const &rhs);
    Actions3D &operator=(Actions3D const &rhs);
    void swap(Actions3D &rhs);
    Actions3D *clone() const;
    plint addBlock(MultiBlock3D &block);
    plint addAction(Action3D *action);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, std::vector<plint> blockNums, Box3D domain);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, std::vector<plint> blockNums,
        std::vector<Box3D> const &domains);
    plint addProcessor(BoxProcessingFunctional3D *functional, plint blockNum1, Box3D domain);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, std::vector<Box3D> const &domains);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, Box3D domain);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2,
        std::vector<Box3D> const &domains);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        Box3D domain);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        std::vector<Box3D> const &domains);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        plint blockNum4, Box3D domain);
    plint addProcessor(
        BoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2, plint blockNum3,
        plint blockNum4, plint blockNum5, Box3D domain);
    plint addProcessor(
        ReductiveBoxProcessingFunctional3D *functional, std::vector<plint> blockNums, Box3D domain);
    plint addProcessor(
        ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, Box3D domain);
    plint addProcessor(
        ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2,
        Box3D domain);
    plint addProcessor(
        ReductiveBoxProcessingFunctional3D *functional, plint blockNum1, plint blockNum2,
        plint blockNum3, Box3D domain);
    plint addActions(Actions3D *subActions);
    Actions3D &getActions(plint id);
    plint addInternalProcessors(plint blockNum, plint level);  // Does not communicate.
    plint addAllInternalProcessorsAndComm(plint blockNum);     // Communicates.
    plint addCommunication(plint blockNum, modif::ModifT whichData);
    plint addEvaluateStats(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addCollideAndStream(plint blockNum, Box3D domain);
    template <typename T, template <typename U> class Descriptor>
    plint addCollideAndStream(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addStream(plint blockNum, Box3D domain);
    template <typename T, template <typename U> class Descriptor>
    plint addStream(plint blockNum);
    template <typename T, template <typename U> class Descriptor>
    plint addIncrementTime(plint blockNum);
    plint addSumReductionToVariable(plint actionID, double &result);
    template <class Function>
    plint addExecuteFunction(Function f);
    plint addNOP();
    void execute();
    void execute(plint actionId);
    void execute(plint actionFrom, plint actionTo);
    BlockStatistics const &statistics(plint actionId) const;

private:
    std::vector<id_t> allMultiBlocks;
    std::vector<Action3D *> actions;
};

class CopySumReductionToVariable3D : public Action3D {
public:
    CopySumReductionToVariable3D(Actions3D const &actions_, plint actionID_, double &result_);
    virtual CopySumReductionToVariable3D *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    Actions3D const &actions;
    plint actionID;
    double &result;
};

template <class Function>
class ExecuteFunctionAction : public Action3D {
public:
    ExecuteFunctionAction(Function f_);
    virtual ExecuteFunctionAction<Function> *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);

private:
    Function f;
};

class NOPaction : public Action3D {
public:
    virtual NOPaction *clone() const;
    virtual void execute(std::vector<id_t> &allMultiBlocks);
};

}  // namespace plb

#endif  // COUPLING_3D_H
