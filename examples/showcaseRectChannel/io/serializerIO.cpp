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

#include "io/serializerIO.h"

#include <fstream>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/plbProfiler.h"
#include "io/base64.h"
#include "io/base64.hh"
#include "io/endianness.h"

namespace plb {

/* *************** Class Base64Writer ******************************** */

class Base64Writer : public SerializedWriter {
public:
    Base64Writer(std::ostream *ostr_, bool enforceUint_, bool switchEndianness_);
    Base64Writer(Base64Writer const &rhs);
    virtual Base64Writer *clone() const;
    ~Base64Writer();
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(char const *dataBuffer, pluint bufferSize);

private:
    std::ostream *ostr;
    bool enforceUint;
    bool switchEndianness;
    Base64Encoder<char> *dataEncoder;
};

Base64Writer::Base64Writer(std::ostream *ostr_, bool enforceUint_, bool switchEndianness_) :
    ostr(ostr_), enforceUint(enforceUint_), switchEndianness(switchEndianness_), dataEncoder(0)
{ }

Base64Writer::Base64Writer(Base64Writer const &rhs) :
    ostr(rhs.ostr), enforceUint(rhs.enforceUint), dataEncoder(0)
{
    if (rhs.dataEncoder) {
        dataEncoder = new Base64Encoder<char>(*rhs.dataEncoder);
    }
}

Base64Writer *Base64Writer::clone() const
{
    return new Base64Writer(*this);
}

Base64Writer::~Base64Writer()
{
    delete dataEncoder;
}

void Base64Writer::writeHeader(pluint dataSize)
{
    PLB_PRECONDITION(ostr && (bool)(*ostr));
    if (enforceUint) {
        Base64Encoder<unsigned int> sizeEncoder(*ostr, 1);
        PLB_PRECONDITION(dataSize <= std::numeric_limits<unsigned int>::max());
        unsigned int uintBinarySize = (unsigned int)dataSize;
        if (switchEndianness) {
            endianByteSwap(uintBinarySize);
        }
        sizeEncoder.encode(&uintBinarySize, 1);
    } else {
        Base64Encoder<pluint> sizeEncoder(*ostr, 1);
        if (switchEndianness) {
            endianByteSwap(dataSize);
        }
        sizeEncoder.encode(&dataSize, 1);
    }
    dataEncoder = new Base64Encoder<char>(*ostr, dataSize);
}

void Base64Writer::writeData(char const *dataBuffer, pluint bufferSize)
{
    global::profiler().start("io");
    dataEncoder->encode(dataBuffer, bufferSize);
    global::profiler().stop("io");
}

/* *************** Class RawBinaryWriter ******************************** */

class RawBinaryWriter : public SerializedWriter {
public:
    RawBinaryWriter(std::ostream *ostr_);
    RawBinaryWriter(RawBinaryWriter const &rhs);
    virtual RawBinaryWriter *clone() const;
    ~RawBinaryWriter();
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(char const *dataBuffer, pluint bufferSize);

private:
    std::ostream *ostr;
};

RawBinaryWriter::RawBinaryWriter(std::ostream *ostr_) : ostr(ostr_) { }

RawBinaryWriter::RawBinaryWriter(RawBinaryWriter const &rhs) : ostr(rhs.ostr) { }

RawBinaryWriter *RawBinaryWriter::clone() const
{
    return new RawBinaryWriter(*this);
}

RawBinaryWriter::~RawBinaryWriter() { }

void RawBinaryWriter::writeHeader(pluint dataSize)
{
    ostr->write((char const *)&dataSize, sizeof(dataSize));
}

void RawBinaryWriter::writeData(char const *dataBuffer, pluint bufferSize)
{
    global::profiler().start("io");
    ostr->write(dataBuffer, (int)bufferSize);
    global::profiler().stop("io");
}

/* *************** Class Base64Reader ******************************** */

class Base64Reader : public SerializedReader {
public:
    Base64Reader(std::istream *istr_, bool enforceUint_, bool switchEndianness_);
    Base64Reader(Base64Reader const &rhs);
    virtual Base64Reader *clone() const;
    ~Base64Reader();
    virtual void readHeader(pluint dataSize) const;
    virtual void readData(char *dataBuffer, pluint bufferSize) const;

private:
    std::istream *istr;
    bool enforceUint;
    bool switchEndianness;
    mutable Base64Decoder<char> *dataDecoder;
};

Base64Reader::Base64Reader(std::istream *istr_, bool enforceUint_, bool switchEndianness_) :
    istr(istr_), enforceUint(enforceUint_), switchEndianness(switchEndianness_), dataDecoder(0)
{ }

Base64Reader::Base64Reader(Base64Reader const &rhs) :
    istr(rhs.istr), enforceUint(rhs.enforceUint), dataDecoder(0)
{
    if (rhs.dataDecoder) {
        dataDecoder = new Base64Decoder<char>(*rhs.dataDecoder);
    }
}

Base64Reader *Base64Reader::clone() const
{
    return new Base64Reader(*this);
}

Base64Reader::~Base64Reader()
{
    delete dataDecoder;
}

void Base64Reader::readHeader(pluint dataSize) const
{
    PLB_PRECONDITION(istr && (bool)(*istr));
    pluint binarySize = 0;
    if (enforceUint) {
        unsigned int uintBinarySize;
        Base64Decoder<unsigned int> sizeDecoder(*istr, 1);
        sizeDecoder.decode(&uintBinarySize, 1);
        if (switchEndianness) {
            endianByteSwap(uintBinarySize);
        }
        binarySize = uintBinarySize;
    } else {
        Base64Decoder<pluint> sizeDecoder(*istr, 1);
        sizeDecoder.decode(&binarySize, 1);
        if (switchEndianness) {
            endianByteSwap(binarySize);
        }
    }
    PLB_PRECONDITION(binarySize == dataSize);

    dataDecoder = new Base64Decoder<char>(*istr, dataSize);
}

void Base64Reader::readData(char *dataBuffer, pluint bufferSize) const
{
    global::profiler().start("io");
    dataDecoder->decode(dataBuffer, bufferSize);
    global::profiler().stop("io");
}

/* *************** Free functions ************************************ */

void serializerToBase64Stream(
    DataSerializer const *serializer, std::ostream *ostr, bool enforceUint, bool mainProcOnly)
{
    serializerToSink(
        serializer,
        new Base64Writer(ostr, enforceUint, global::IOpolicy().getEndianSwitchOnBase64out()),
        mainProcOnly);
}

void serializerToRawBinaryStream(
    DataSerializer const *serializer, std::ostream *ostr, bool mainProcOnly)
{
    serializerToSink(serializer, new RawBinaryWriter(ostr), mainProcOnly);
}

void base64StreamToUnSerializer(
    std::istream *istr, DataUnSerializer *unSerializer, bool enforceUint, bool mainProcOnly)
{
    sourceToUnSerializer(
        new Base64Reader(istr, enforceUint, global::IOpolicy().getEndianSwitchOnBase64in()),
        unSerializer, mainProcOnly);
}

}  // namespace plb
