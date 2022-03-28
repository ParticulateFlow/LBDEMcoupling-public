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
 * Input/Output in XML format -- non-generic code.
 */

#include "libraryInterfaces/TINYXML_xmlIO.h"

#include <algorithm>
#include <cctype>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "io/parallelIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "parallelism/mpiManager.h"

namespace plb {

XMLreader XMLreader::notFound;

XMLreader::XMLreader(std::vector<TiXmlNode *> pParentVect)
{
    if (global::mpi().isMainProcessor()) {
        mainProcessorIni(pParentVect);
    } else {
        slaveProcessorIni();
    }
}

XMLreader::XMLreader(std::string fName)
{
    TiXmlDocument *doc = 0;
    bool loadOK = false;
    std::string errorMessage;
    if (global::mpi().isMainProcessor()) {
        doc = new TiXmlDocument(fName.c_str());
        loadOK = doc->LoadFile();
        if (!loadOK) {
            errorMessage = "Problem processing input XML file " + fName + ": " + doc->ErrorDesc();
        }
    }

    global::mpi().bCast(errorMessage);
    plbMainProcIOError(!loadOK, errorMessage);

    if (global::mpi().isMainProcessor()) {
        mainProcessorIni(doc);
        delete doc;
    } else {
        slaveProcessorIni();
    }
}

// Transform a string into a xml object
void XMLreader::XMLreader_parse_from_string(const char *raw_xml)
{
    TiXmlDocument *doc = new TiXmlDocument();
    doc->Parse(raw_xml, 0, TIXML_ENCODING_UTF8);
    if (global::mpi().isMainProcessor()) {
        mainProcessorIni(doc);
    } else {
        slaveProcessorIni();
    }
}

void XMLreader::mainProcessorIni(TiXmlNode *pParent)
{
    std::vector<TiXmlNode *> pParentVect;
    pParentVect.push_back(pParent);
    mainProcessorIni(pParentVect);
}

void XMLreader::mainProcessorIni(std::vector<TiXmlNode *> pParentVect)
{
    std::map<plint, TiXmlNode *> parents;
    for (pluint iParent = 0; iParent < pParentVect.size(); ++iParent) {
        PLB_PRECONDITION(
            pParentVect[iParent]->Type() == TiXmlNode::DOCUMENT
            || pParentVect[iParent]->Type() == TiXmlNode::ELEMENT);

        TiXmlElement *pParentElement = pParentVect[iParent]->ToElement();
        plint id = 0;
        if (pParentElement) {
            const char *attribute = pParentElement->Attribute("id");
            if (attribute) {
                std::stringstream attributestr(attribute);
                attributestr >> id;
            }
        }
        parents[id] = pParentVect[iParent];
    }

    plint numId = (plint)parents.size();
    global::mpi().bCast(&numId, 1);

    std::map<plint, TiXmlNode *>::iterator it = parents.begin();
    name = it->second->ValueStr();
    global::mpi().bCast(name);

    for (; it != parents.end(); ++it) {
        plint id = it->first;
        global::mpi().bCast(&id, 1);

        TiXmlNode *pParent = it->second;
        Data &data = data_map[id];
        data.text = "";

        typedef std::map<std::string, std::vector<TiXmlNode *> > ChildMap;
        ChildMap childMap;
        TiXmlNode *pChild;
        for (pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) {
            int type = pChild->Type();
            if (type == TiXmlNode::ELEMENT) {
                std::string name(pChild->Value());
                childMap[name].push_back(pChild);
            } else if (type == TiXmlNode::TEXT) {
                data.text = pChild->ToText()->ValueStr();
            }
        }
        global::mpi().bCast(data.text);
        plint numChildren = (plint)childMap.size();
        global::mpi().bCast(&numChildren, 1);

        for (ChildMap::iterator it = childMap.begin(); it != childMap.end(); ++it) {
            std::vector<TiXmlNode *> pChildVect = it->second;
            data.children.push_back(new XMLreader(pChildVect));
        }
    }
}

void XMLreader::slaveProcessorIni()
{
    plint numId = 0;
    global::mpi().bCast(&numId, 1);
    global::mpi().bCast(name);

    for (plint iId = 0; iId < numId; ++iId) {
        plint id = 0;
        global::mpi().bCast(&id, 1);
        Data &data = data_map[id];
        global::mpi().bCast(data.text);
        plint numChildren = 0;
        global::mpi().bCast(&numChildren, 1);

        for (plint iChild = 0; iChild < numChildren; ++iChild) {
            std::vector<TiXmlNode *> noParam;
            data.children.push_back(new XMLreader(noParam));
        }
    }
}

XMLreader::XMLreader()
{
    name = "XML node not found";
}

XMLreader::~XMLreader()
{
    std::map<plint, Data>::iterator it = data_map.begin();
    for (; it != data_map.end(); ++it) {
        std::vector<XMLreader *> &children = it->second.children;
        for (pluint iNode = 0; iNode < children.size(); ++iNode) {
            delete children[iNode];
        }
    }
}

void XMLreader::print(int indent) const
{
    std::string indentStr(indent, ' ');
    pcout << indentStr << "[" << name << "]" << std::endl;
    std::string text = getFirstText();
    if (!text.empty()) {
        pcout << indentStr << "  " << text << std::endl;
    }
    std::vector<XMLreader *> const &children = data_map.begin()->second.children;
    for (pluint iNode = 0; iNode < children.size(); ++iNode) {
        children[iNode]->print(indent + 2);
    }
}

XMLreaderProxy XMLreader::operator[](std::string name) const
{
    Data const &data = data_map.begin()->second;
    for (pluint iNode = 0; iNode < data.children.size(); ++iNode) {
        if (data.children[iNode]->name == name) {
            return XMLreaderProxy(data.children[iNode]);
        }
    }
    plbIOError(std::string("Element ") + name + std::string(" not found in XML file."));
    return XMLreaderProxy(0);
}

XMLreaderProxy XMLreader::getElement(std::string name, plint id) const
{
    std::map<plint, Data>::const_iterator it = data_map.find(id);
    if (it == data_map.end()) {
        std::stringstream idStr;
        idStr << id;
        plbIOError(std::string("Element with id ") + idStr.str() + std::string(" does not exist"));
    }
    std::vector<XMLreader *> const &children = it->second.children;
    for (pluint iNode = 0; iNode < children.size(); ++iNode) {
        if (children[iNode]->name == name) {
            return XMLreaderProxy(children[iNode]);
        }
    }
    plbIOError(std::string("Element ") + name + std::string(" not found in XML file."));
    return XMLreaderProxy(0);
}

std::string XMLreader::getName() const
{
    return name;
}

std::string XMLreader::getText() const
{
    return data_map.begin()->second.text;
}

std::string XMLreader::getText(plint id) const
{
    std::map<plint, Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        return it->second.text;
    } else {
        return "";
    }
}

plint XMLreader::getFirstId() const
{
    return data_map.begin()->first;
}

std::string XMLreader::getFirstText() const
{
    return data_map.begin()->second.text;
}

bool XMLreader::idExists(plint id) const
{
    std::map<plint, Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        return true;
    } else {
        return false;
    }
}

bool XMLreader::getNextId(plint &id) const
{
    std::map<plint, Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        ++it;
        if (it != data_map.end()) {
            id = it->first;
            return true;
        }
    }
    return false;
}

std::vector<XMLreader *> const &XMLreader::getChildren(plint id) const
{
    std::map<plint, Data>::const_iterator it = data_map.find(id);
    if (it == data_map.end()) {
        plbIOError(
            std::string("Cannot access id ") + util::val2str(id) + " in XML element " + name);
    }
    return it->second.children;
}

XMLreaderProxy::XMLreaderProxy(XMLreader const *reader_) : reader(reader_)
{
    if (reader) {
        id = reader->getFirstId();
    } else {
        id = 0;
    }
}

XMLreaderProxy::XMLreaderProxy(XMLreader const *reader_, plint id_) : reader(reader_), id(id_) { }

XMLreaderProxy XMLreaderProxy::operator[](std::string name) const
{
    if (!reader) {
        plbIOError(std::string("Cannot read value from XML element ") + name);
    }
    return reader->getElement(name, id);
}

XMLreaderProxy XMLreaderProxy::operator[](plint newId) const
{
    if (!reader) {
        plbIOError(std::string("Cannot read value from XML element"));
    }
    if (!reader->idExists(newId)) {
        std::stringstream newIdStr;
        newIdStr << newId;
        plbIOError(
            std::string("Id ") + newIdStr.str() + std::string(" does not exist in XML element"));
    }
    return XMLreaderProxy(reader, newId);
}

bool XMLreaderProxy::isValid() const
{
    return reader;
}

plint XMLreaderProxy::getId() const
{
    return id;
}

XMLreaderProxy XMLreaderProxy::iterId() const
{
    if (!reader) {
        plbIOError(std::string("Use of invalid XML element"));
    }
    plint newId = id;
    if (reader->getNextId(newId)) {
        return XMLreaderProxy(reader, newId);
    } else {
        return XMLreaderProxy(0);
    }
}

std::string XMLreaderProxy::getName() const
{
    if (!reader) {
        plbIOError(std::string("Cannot read value from XML element "));
    }
    return reader->getName();
}

std::vector<XMLreader *> const &XMLreaderProxy::getChildren() const
{
    if (!reader) {
        plbIOError(std::string("Cannot read value from XML element "));
    }
    return reader->getChildren(id);
}

XMLwriter::XMLwriter() : isDocument(true), currentId(0) { }

XMLwriter::~XMLwriter()
{
    std::map<plint, Data>::iterator it = data_map.begin();
    for (; it != data_map.end(); ++it) {
        std::vector<XMLwriter *> &children = it->second.children;
        for (pluint iNode = 0; iNode < children.size(); ++iNode) {
            delete children[iNode];
        }
    }
}

XMLwriter::XMLwriter(std::string name_) : isDocument(false), name(name_), currentId(0) { }

void XMLwriter::setString(std::string const &value)
{
    data_map[currentId].text = value;
}

XMLwriter &XMLwriter::operator[](std::string name)
{
    std::vector<XMLwriter *> &children = data_map[currentId].children;
    // If node already exists, simply return it.
    for (pluint iNode = 0; iNode < data_map[currentId].children.size(); ++iNode) {
        if (data_map[currentId].children[iNode]->name == name) {
            return *children[iNode];
        }
    }
    // Else, create and return it.
    children.push_back(new XMLwriter(name));
    return *children.back();
}

XMLwriter &XMLwriter::operator[](plint id)
{
    currentId = id;
    return *this;
}

void XMLwriter::print(std::string fName) const
{
    plb_ofstream ofile(fName.c_str());
    plbIOError(
        !ofile.is_open(),
        std::string("Could not open file ") + fName + std::string(" for write access"));
    toOutputStream(ofile);
}

std::string XMLwriter::sprint()
{
    std::stringstream a;
    toOutputStream_parrallel(a);
    return a.str();
}

}  // namespace plb
