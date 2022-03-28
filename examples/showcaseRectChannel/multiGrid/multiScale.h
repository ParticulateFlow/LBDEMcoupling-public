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
 * Multi scale manager -- header file.
 */
#ifndef MULTI_SCALE_H
#define MULTI_SCALE_H

#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"

namespace plb {

template <typename T>
T from_reference(T value, int dimDx, int dimDt, int dxScale, int dtScale)
{
    int exponent = dxScale * dimDx + dtScale * dimDt;
    return value * util::twoToThePower(exponent);
}

template <typename T>
T to_reference(T value, int dimDx, int dimDt, int dxScale, int dtScale)
{
    int exponent = -dxScale * dimDx - dtScale * dimDt;
    return value * util::twoToThePower(exponent);
}

class MultiBlockManagement2D;
struct Box2D;

class MultiBlockManagement3D;
struct Box3D;

struct MultiScaleManager {
    virtual ~MultiScaleManager() { }
    virtual MultiScaleManager *clone() const = 0;

    virtual Box2D scaleBox(Box2D box, plint nLevel) const = 0;
    virtual Box3D scaleBox(Box3D box, plint nLevel) const = 0;

    virtual MultiBlockManagement2D scaleMultiBlockManagement(
        MultiBlockManagement2D const &multiBlockManagement, plint nLevel) const = 0;
    virtual MultiBlockManagement3D scaleMultiBlockManagement(
        MultiBlockManagement3D const &multiBlockManagement, plint nLevel) const = 0;

    virtual void scaleVelocity(Array<double, 2> &u, plint nLevel) const = 0;
    virtual void scaleVelocity(Array<double, 3> &u, plint nLevel) const = 0;

    virtual double scaleDeltaX(plint nLevel) const = 0;
    virtual double scaleDeltaT(plint nLevel) const = 0;
};

class PowerTwoMultiScaleManager : public MultiScaleManager {
public:
    virtual Box2D scaleBox(Box2D box, plint nLevel) const;
    virtual Box3D scaleBox(Box3D box, plint nLevel) const;

    virtual MultiBlockManagement2D scaleMultiBlockManagement(
        MultiBlockManagement2D const &multiBlockManagement, plint nLevel) const;
    virtual MultiBlockManagement3D scaleMultiBlockManagement(
        MultiBlockManagement3D const &multiBlockManagement, plint nLevel) const;

public:
    static plint twoToTheLevel(plint nLevel);
};

class ConvectiveMultiScaleManager : public PowerTwoMultiScaleManager {
public:
    virtual ConvectiveMultiScaleManager *clone() const;
    virtual void scaleVelocity(Array<double, 2> &u, plint nLevel) const;
    virtual void scaleVelocity(Array<double, 3> &u, plint nLevel) const;
    virtual double scaleDeltaX(plint nLevel) const;
    virtual double scaleDeltaT(plint nLevel) const;
};

namespace global {

class DefaultMultiScaleManager {
private:
    DefaultMultiScaleManager();
    ~DefaultMultiScaleManager();
    void set(MultiScaleManager *newManager);
    MultiScaleManager const &get() const;

private:
    MultiScaleManager *defaultManager;

    friend DefaultMultiScaleManager &accessDefaultMultiScaleManager();
    friend MultiScaleManager const &getDefaultMultiScaleManager();
    friend void setDefaultMultiScaleManager(MultiScaleManager *newMultiScaleManager);
};

DefaultMultiScaleManager &accessDefaultMultiScaleManager();

MultiScaleManager const &getDefaultMultiScaleManager();

void setDefaultMultiScaleManager(MultiScaleManager *newMultiScaleManager);

}  // namespace global

}  // namespace plb

#endif  // MULTI_SCALE_H
