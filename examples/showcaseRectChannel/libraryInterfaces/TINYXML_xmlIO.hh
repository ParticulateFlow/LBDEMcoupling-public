/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
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
 * Input/Output in XML format -- generic code.
 */

#ifndef XML_IO_HH
#define XML_IO_HH

#include <cctype>
#include <iomanip>
#include <typeinfo>

#include "core/globalDefs.h"
#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "io/parallelIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "parallelism/mpiManager.h"

namespace plb {

template <typename T>
void XMLreaderProxy::read(T &value) const
{
    if (!reader)
        return;
    std::stringstream valueStr(reader->getText(id));
    T tmp = T();
    if (!(valueStr >> tmp)) {
        plbIOError(std::string("Cannot read value from XML element ") + reader->getName());
    }
    value = tmp;
}

template <>
inline void XMLreaderProxy::read<bool>(bool &value) const
{
    if (!reader)
        return;
    std::stringstream valueStr(reader->getText(id));
    std::string word;
    valueStr >> word;
    // Transform to lower-case, so that "true" and "false" are case-insensitive.
    word = util::tolower(word);
    if (word == "true") {
        value = true;
    } else if (word == "false") {
        value = false;
    } else {
        plbIOError(std::string("Cannot read boolean value from XML element ") + reader->getName());
    }
}

template <>
inline void XMLreaderProxy::read<std::string>(std::string &entry) const
{
    if (!reader)
        return;
    entry = reader->getText(id);
}

template <typename T>
bool XMLreaderProxy::readNoThrow(T &value) const
{
    if (!reader)
        return false;
    std::stringstream valueStr(reader->getText(id));
    T tmp = T();
    if (!(valueStr >> tmp)) {
        return false;
    }
    value = tmp;
    return true;
}

template <>
inline bool XMLreaderProxy::readNoThrow<bool>(bool &value) const
{
    if (!reader)
        return false;
    std::stringstream valueStr(reader->getText(id));
    std::string word;
    valueStr >> word;
    // Transform to lower-case, so that "true" and "false" are case-insensitive.
    word = util::tolower(word);
    if (word == "true") {
        value = true;
        return true;
    } else if (word == "false") {
        value = false;
        return true;
    }
    return false;
}

template <>
inline bool XMLreaderProxy::readNoThrow<std::string>(std::string &entry) const
{
    if (!reader)
        return false;
    entry = reader->getText(id);
    return true;
}

template <typename T>
void XMLreaderProxy::read(std::vector<T> &values) const
{
    if (!reader)
        return;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    std::vector<T> tmp(values);
    while (multiValueStr >> word) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            plbIOError(
                std::string("Cannot read value array from XML element ") + reader->getName());
        }
        tmp.push_back(value);
    }
    values.swap(tmp);
}

template <>
inline void XMLreaderProxy::read<bool>(std::vector<bool> &values) const
{
    if (!reader)
        return;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    std::vector<bool> tmp(values);
    while (multiValueStr >> word) {
        bool value = false;
        word = util::tolower(word);
        if (word == "true") {
            value = true;
        } else if (word == "false") {
            value = false;
        } else {
            plbIOError(
                std::string("Cannot read boolean value from XML element ") + reader->getName());
        }
        tmp.push_back(value);
    }
    values.swap(tmp);
}

template <typename T>
bool XMLreaderProxy::readNoThrow(std::vector<T> &values) const
{
    if (!reader)
        return false;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    std::vector<T> tmp(values);
    while (multiValueStr >> word) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            return false;
        }
        tmp.push_back(value);
    }
    values.swap(tmp);
    return true;
}

template <typename T, plint N>
void XMLreaderProxy::read(Array<T, N> &values) const
{
    if (!reader)
        return;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    values.resetToZero();
    plint i = 0;
    while (multiValueStr >> word && i < N) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            plbIOError(
                std::string("Cannot read value array from XML element ") + reader->getName());
        }
        values[i] = value;
        ++i;
    }
}

template <typename T, plint N>
bool XMLreaderProxy::readNoThrow(Array<T, N> &values) const
{
    if (!reader)
        return false;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    plint i = 0;
    while (multiValueStr >> word && i < N) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            return false;
        }
        values[i] = value;
        ++i;
    }
    return true;
}

template <typename T>
void XMLwriter::set(T const &value, plint precision)
{
    std::stringstream valuestr;
    if (precision >= 0) {
        valuestr << std::setprecision(precision);
    }
    valuestr << value;
    valuestr >> data_map[currentId].text;
}

template <typename T>
void XMLwriter::set(std::vector<T> const &values, plint precision)
{
    std::stringstream valuestr;
    if (precision >= 0) {
        valuestr << std::setprecision(precision);
    }
    for (pluint i = 0; i < values.size(); ++i) {
        if (i != 0) {
            valuestr << " ";
        }
        valuestr << values[i];
    }
    data_map[currentId].text = valuestr.str();
}

template <typename T, int N>
void XMLwriter::set(Array<T, N> const &values, plint precision)
{
    std::stringstream valuestr;
    if (precision >= 0) {
        valuestr << std::setprecision(precision);
    }
    for (pluint i = 0; i < N; ++i) {
        if (i != 0) {
            valuestr << " ";
        }
        valuestr << values[i];
    }
    data_map[currentId].text = valuestr.str();
}

template <>
inline void XMLwriter::set<bool>(bool const &value, plint precision)
{
    if (value) {
        data_map[currentId].text = "True";
    } else {
        data_map[currentId].text = "False";
    }
}

template <typename ostrT>
void XMLwriter::toOutputStream_parrallel(ostrT &ostr, int indent) const
{
    if (data_map.empty())
        return;
    if (isDocument) {
        ostr << "<?xml version=\"1.0\" ?>\n";
        std::vector<XMLwriter *> const &children = data_map.begin()->second.children;
        for (pluint iNode = 0; iNode < children.size(); ++iNode) {
            children[iNode]->toOutputStream_parrallel(ostr);
        }
    } else {
        std::map<plint, Data>::const_iterator it = data_map.begin();
        for (; it != data_map.end(); ++it) {
            std::vector<XMLwriter *> const &children = it->second.children;
            std::string const &text = it->second.text;
            std::string indentStr(indent, ' ');
            ostr << indentStr << "<" << name;
            if (data_map.size() > 1 || it->first != 0) {
                ostr << " id=\"" << it->first << "\"";
            }
            if (children.empty()) {
                ostr << ">";
            } else {
                ostr << ">\n";
            }
            if (!text.empty()) {
                if (children.empty()) {
                    ostr << " " << text << " ";
                } else {
                    ostr << indentStr << "    " << text << "\n";
                }
            }
            for (pluint iNode = 0; iNode < children.size(); ++iNode) {
                children[iNode]->toOutputStream_parrallel(ostr, indent + 4);
            }
            if (!children.empty()) {
                ostr << indentStr;
            }
            ostr << "</" << name << ">\n";
        }
    }
}

template <typename ostrT>
void XMLwriter::toOutputStream(ostrT &ostr, int indent) const
{
    if (!global::mpi().isMainProcessor())
        return;
    if (data_map.empty())
        return;

    if (isDocument) {
        ostr << "<?xml version=\"1.0\" ?>\n";
        std::vector<XMLwriter *> const &children = data_map.begin()->second.children;
        for (pluint iNode = 0; iNode < children.size(); ++iNode) {
            children[iNode]->toOutputStream(ostr);
        }
    } else {
        std::map<plint, Data>::const_iterator it = data_map.begin();
        for (; it != data_map.end(); ++it) {
            std::vector<XMLwriter *> const &children = it->second.children;
            std::string const &text = it->second.text;
            std::string indentStr(indent, ' ');
            ostr << indentStr << "<" << name;
            if (data_map.size() > 1 || it->first != 0) {
                ostr << " id=\"" << it->first << "\"";
            }
            if (children.empty()) {
                ostr << ">";
            } else {
                ostr << ">\n";
            }
            if (!text.empty()) {
                if (children.empty()) {
                    ostr << " " << text << " ";
                } else {
                    ostr << indentStr << "    " << text << "\n";
                }
            }
            for (pluint iNode = 0; iNode < children.size(); ++iNode) {
                children[iNode]->toOutputStream(ostr, indent + 4);
            }
            if (!children.empty()) {
                ostr << indentStr;
            }
            ostr << "</" << name << ">\n";
        }
    }
}

}  // namespace plb

#endif  // XML_IO_HH
