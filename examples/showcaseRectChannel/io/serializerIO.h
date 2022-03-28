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

#ifndef SERIALIZER_IO_H
#define SERIALIZER_IO_H

#include <iomanip>
#include <iosfwd>

#include "core/globalDefs.h"
#include "core/serializer.h"

namespace plb {

/// Take a Serializer, convert into Base64 format (ASCII based binary representation), and stream
/// into output stream.
/** Ahead of the data, an integer value is encoded which stands for the total size of
 *  serialized data. For compatibility with the VTK file format (as of this writing),
 *  you can enforce that the type of this variable is converted to "unsigned int".
 *  Note that this may lead to errors on 64-bit platforms, if the total amount of
 *  data exceeds 2 GB.
 */
void serializerToBase64Stream(
    DataSerializer const *serializer, std::ostream *ostr, bool enforceUint = false,
    bool mainProcOnly = true);

void serializerToRawBinaryStream(
    DataSerializer const *serializer, std::ostream *ostr, bool mainProcOnly = true);

/// Take an input stream with Base64 encoded binary content, and stream into an unSerializer
/** If the integer value which indicates the amount of data to be unSerialized is of type
 *  "unsigned int", this fact can be enforced with the flag enforceUplint to ensure
 *  compatibility between 32-bit and 64-bit platforms.
 */
void base64StreamToUnSerializer(
    std::istream *istr, DataUnSerializer *unSerializer, bool enforceUint = false,
    bool mainProcOnly = true);

/// Take a Serializer, convert and stream into output in ASCII format.
/** Number of digits in the ASCII representation of numbers is given by the variable numDigits.
 */
template <typename T>
void serializerToAsciiStream(
    DataSerializer const *serializer, std::ostream *ostr, plint numDigits = 8,
    bool mainProcOnly = true);

/// Take an UnSerializer and fill it with data from an ASCII-format input stream.
template <typename T>
void asciiStreamToUnSerializer(
    std::istream *istr, DataUnSerializer *unSerializer, bool mainProcOnly = true);

/* *************** Class AsciiWriter ******************************** */

template <typename T>
class AsciiWriter : public SerializedWriter {
public:
    AsciiWriter(std::ostream *ostr_, plint numDigits_);
    virtual AsciiWriter<T> *clone() const;
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(char const *dataBuffer, pluint bufferSize);

private:
    std::ostream *ostr;
    plint numDigits;
};

/* *************** Class AsciiReader ******************************** */

template <typename T>
class AsciiReader : public SerializedReader {
public:
    AsciiReader(std::istream *istr_);
    virtual AsciiReader<T> *clone() const;
    virtual void readHeader(pluint dataSize) const;
    virtual void readData(char *dataBuffer, pluint bufferSize) const;

private:
    std::istream *istr;
};

}  // namespace plb

#endif  // SERIALIZER_IO_H
