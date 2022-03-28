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

#ifndef SERIALIZER_IO_HH
#define SERIALIZER_IO_HH

#include <iomanip>
#include <iosfwd>

#include "core/serializer.h"
#include "io/serializerIO.h"

namespace plb {

/* *************** Class AsciiWriter ******************************** */

template <typename T>
AsciiWriter<T>::AsciiWriter(std::ostream *ostr_, plint numDigits_) :
    ostr(ostr_), numDigits(numDigits_)
{ }

template <typename T>
AsciiWriter<T> *AsciiWriter<T>::clone() const
{
    return new AsciiWriter<T>(*this);
}

template <typename T>
void AsciiWriter<T>::writeHeader(pluint dataSize)
{ }

template <typename T>
void AsciiWriter<T>::writeData(char const *dataBuffer, pluint bufferSize)
{
    T tmp;
    for (pluint iData = 0; iData < bufferSize; iData += sizeof(T)) {
        for (pluint iVal = 0; iVal < sizeof(T); ++iVal) {
            *((char *)&tmp + iVal) = dataBuffer[iData + iVal];
        }
        if (numDigits == 0) {
            (*ostr) << tmp << " ";
        } else {
            (*ostr) << std::setprecision(numDigits) << tmp << " ";
        }
    }
}

/* *************** Class AsciiReader ******************************** */

template <typename T>
AsciiReader<T>::AsciiReader(std::istream *istr_) : istr(istr_)
{ }

template <typename T>
AsciiReader<T> *AsciiReader<T>::clone() const
{
    return new AsciiReader<T>(*this);
}

template <typename T>
void AsciiReader<T>::readHeader(pluint dataSize) const
{ }

template <typename T>
void AsciiReader<T>::readData(char *dataBuffer, pluint bufferSize) const
{
    T tmp;
    // Read the input data.
    for (pluint iData = 0; iData < bufferSize; iData += sizeof(T)) {
        (*istr) >> tmp;
        for (pluint iVal = 0; iVal < sizeof(T); ++iVal) {
            dataBuffer[iData + iVal] = *((char *)&tmp + iVal);
        }
    }
}

template <typename T>
void serializerToAsciiStream(
    DataSerializer const *serializer, std::ostream *ostr, plint numDigits, bool mainProcOnly)
{
    serializerToSink(serializer, new AsciiWriter<T>(ostr, numDigits), mainProcOnly);
}

template <typename T>
void asciiStreamToUnSerializer(
    std::istream *istr, DataUnSerializer *unSerializer, bool mainProcOnly)
{
    sourceToUnSerializer(new AsciiReader<T>(istr), unSerializer, mainProcOnly);
}

}  // namespace plb

#endif  // SERIALIZER_IO_HH
