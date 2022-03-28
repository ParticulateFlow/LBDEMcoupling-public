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
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_WRAPPER_2D_H
#define NTENSOR_ANALYSIS_WRAPPER_2D_H

#include <memory>

#include "atomicBlock/dataField2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisFunctional2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

template <typename T1, typename T2>
void copy(MultiNTensorField2D<T1> &field, MultiNTensorField2D<T2> &convertedField, Box2D domain);

template <typename T1, typename T2>
MultiNTensorField2D<T2> *copyConvert(MultiNTensorField2D<T1> &field, Box2D domain);

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_WRAPPER_2D_H
