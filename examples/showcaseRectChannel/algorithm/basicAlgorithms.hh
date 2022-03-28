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

#ifndef BASIC_ALGORITHMS_HH
#define BASIC_ALGORITHMS_HH

#include <cstdio>
#include <vector>

#include "algorithm/basicAlgorithms.h"
#include "core/globalDefs.h"
#include "io/imageWriter.h"
#include "parallelism/mpiManager.h"

namespace plb {

/* ************** class PIDController ********************* */

template <typename T>
PIDController<T>::PIDController() : error((T)0), sumErrors((T)0), deltaError((T)0), oldError((T)0)
{ }

template <typename T>
inline T PIDController<T>::operator()(T target, T current, T Kp, T Ki, T Kd)
{
    error = target - current;
    sumErrors += error;
    deltaError = error - oldError;
    oldError = error;
    return (Kp * error + Ki * sumErrors + Kd * deltaError);
}

template <typename T>
void PIDController<T>::saveState(std::string baseFileName, plint fileNamePadding, plint iIter)
{
    if (global::mpi().isMainProcessor()) {
        std::string fn = createFileName(baseFileName, iIter, fileNamePadding) + ".dat";
        FILE *fp = fopen(fn.c_str(), "wb");
        T data[4];
        data[0] = error;
        data[1] = sumErrors;
        data[2] = deltaError;
        data[3] = oldError;
        (void)fwrite(data, sizeof(T), 4, fp);
        fclose(fp);
    }
}

template <typename T>
void PIDController<T>::loadState(std::string baseFileName, plint fileNamePadding, plint iIter)
{
    T data[4];
    if (global::mpi().isMainProcessor()) {
        std::string fn = createFileName(baseFileName, iIter, fileNamePadding) + ".dat";
        FILE *fp = fopen(fn.c_str(), "rb");
        (void)fread(data, sizeof(T), 4, fp);
        fclose(fp);
    }
    global::mpi().bCast<char>((char *)&data[0], 4 * sizeof(T));
    error = data[0];
    sumErrors = data[1];
    deltaError = data[2];
    oldError = data[3];
}

/* ************** class Relaxation ********************* */

template <typename T>
Relaxation<T>::Relaxation() : omega((T)0), equilibrium((T)0), previous((T)0), next((T)0)
{ }

template <typename T>
Relaxation<T>::Relaxation(T omega_, T equilibrium_, T initialValue_) :
    omega(omega_), equilibrium(equilibrium_), previous(initialValue_), next((T)0)
{ }

template <typename T>
inline void Relaxation<T>::setOmega(T omega_)
{
    omega = omega_;
}

template <typename T>
inline void Relaxation<T>::setEquilibrium(T equilibrium_)
{
    equilibrium = equilibrium_;
}

template <typename T>
inline void Relaxation<T>::setInitialValue(T initialValue_)
{
    previous = initialValue_;
}

template <typename T>
inline T Relaxation<T>::iterate()
{
    next = ((T)1 - omega) * previous + omega * equilibrium;
    previous = next;
    return next;
}

template <typename T>
void Relaxation<T>::saveState(std::string baseFileName, plint fileNamePadding, plint iIter)
{
    if (global::mpi().isMainProcessor()) {
        std::string fn = createFileName(baseFileName, iIter, fileNamePadding) + ".dat";
        FILE *fp = fopen(fn.c_str(), "wb");
        T data[4];
        data[0] = omega;
        data[1] = equilibrium;
        data[2] = previous;
        data[3] = next;
        (void)fwrite(data, sizeof(T), 4, fp);
        fclose(fp);
    }
}

template <typename T>
void Relaxation<T>::loadState(std::string baseFileName, plint fileNamePadding, plint iIter)
{
    T data[4];
    if (global::mpi().isMainProcessor()) {
        std::string fn = createFileName(baseFileName, iIter, fileNamePadding) + ".dat";
        FILE *fp = fopen(fn.c_str(), "rb");
        (void)fread(data, sizeof(T), 4, fp);
        fclose(fp);
    }
    global::mpi().bCast<char>((char *)&data[0], 4 * sizeof(T));
    omega = data[0];
    equilibrium = data[1];
    previous = data[2];
    next = data[3];
}

}  // namespace plb

#endif  // BASIC_ALGORITHMS_HH
