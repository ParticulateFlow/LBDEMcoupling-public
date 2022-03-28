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

#ifndef BASIC_ALGORITHMS_H
#define BASIC_ALGORITHMS_H

#include <vector>

#include "core/globalDefs.h"

namespace plb {

namespace algorithm {

std::vector<plint> primeFactor(plint value);

std::vector<plint> evenRepartition(plint value, plint d);

}  // namespace algorithm

template <typename T>
class PIDController {
public:
    PIDController();
    // Kp: proportional coefficient.
    // Ki: integral coefficient.
    // Kd: derivative coefficient.
    T operator()(T target, T current, T Kp = 0.8, T Ki = 0.2, T Kd = 0.2);
    void saveState(std::string baseFileName, plint fileNamePadding, plint iIter);
    void loadState(std::string baseFileName, plint fileNamePadding, plint iIter);

private:
    T error;
    T sumErrors;
    T deltaError;
    T oldError;
};

template <typename T>
class Relaxation {
public:
    Relaxation();
    Relaxation(T omega_, T equilibrium_ = (T)0, T initialValue_ = (T)0);
    void setOmega(T omega_);
    void setEquilibrium(T equilibrium_);
    void setInitialValue(T initialValue_);
    T iterate();
    void saveState(std::string baseFileName, plint fileNamePadding, plint iIter);
    void loadState(std::string baseFileName, plint fileNamePadding, plint iIter);

private:
    T omega;
    T equilibrium;
    T previous;
    T next;
};

}  // namespace plb

#endif  // BASIC_ALGORITHMS_H
