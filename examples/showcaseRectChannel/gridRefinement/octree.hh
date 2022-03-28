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

#ifndef OCTREE_HH
#define OCTREE_HH

#include <cstring>
#include <vector>

#include "gridRefinement/octree.h"

namespace plb {

template <typename DataT>
OctreeNode<DataT>::OctreeNode(OctreeNode<DataT> *parent_, DataT *data_) :
    data(data_), parent(parent_)
{
    level =
        (parent != 0
             ? parent->level + 1
             : 0);  // The level must be >= 0 for normal octrees. Negative levels are reserved.
    isLeaf = true;
    memset(child, 0, sizeof child);
}

template <typename DataT>
int octreeChildType(OctreeNode<DataT> *node)
{
    int o = -1;
    if (node->parent != 0) {
        if (node->parent->level >= 0) {  // Normal octree.
            for (int i = 0; i < 8; i++) {
                if (node->parent->child[i] == node) {
                    o = i;
                    break;
                }
            }
        } else if (node->parent->level == -1)
        {  // The levels -1 and -2 are reserved for the octree periodic extensions.
            o = 7;
        } else if (node->parent->level == -2) {
            o = 0;
        }
    }
    return (o);
}

template <typename DataT>
void freeOctree(OctreeNode<DataT> *root)
{
    // Bottom-Up deallocation: Start deallocating at the leaves.
    for (int i = 0; i < 8; i++) {
        OctreeNode<DataT> *child = root->child[i];
        if (child != 0) {
            if (child->isLeaf) {
                delete child->data;
                delete child;
            } else {
                freeOctree(child);
            }
        }
    }
    delete root->data;
    delete root;
}

template <typename DataT, class AssignData>
void splitOctreeNode(OctreeNode<DataT> *node, AssignData &assignData)
{
    if (node->isLeaf) {
        for (int i = 0; i < 8; i++) {
            OctreeNode<DataT> *child =
                new OctreeNode<DataT>(node, 0);  // Do not know how to split the data at this point.
            node->child[i] = child;
            assignData(
                node->child[i]);  // Here the user-defined data is assigned to the child node.
        }
        node->isLeaf = false;
    }
}

template <typename DataT>
void getMinMaxOctreeLeafNodeLevels(OctreeNode<DataT> *root, int &minLevel, int &maxLevel)
{
    // minLevel and maxLevel must be properly initialized before calling this function.
    if (root == 0) {
        return;
    }

    if (root->isLeaf) {
        minLevel = std::min(minLevel, root->level);
        maxLevel = std::max(maxLevel, root->level);
    } else {
        for (int i = 0; i < 8; i++) {
            OctreeNode<DataT> *child = root->child[i];
            if (child != 0) {
                getMinMaxOctreeLeafNodeLevels(child, minLevel, maxLevel);
            }
        }
    }
}

template <typename DataT, class Process>
void processOctreePreOrder(OctreeNode<DataT> *root, Process &process)
{
    // Top-Down processing: Start processing at the root.
    process(root);
    for (int i = 0; i < 8; i++) {
        OctreeNode<DataT> *child = root->child[i];
        if (child != 0) {
            processOctreePreOrder(child, process);
        }
    }
}

template <typename DataT, class Process>
void processOctreePostOrder(OctreeNode<DataT> *root, Process &process)
{
    // Bottom-Up processing: Start processing at the leaves.
    for (int i = 0; i < 8; i++) {
        OctreeNode<DataT> *child = root->child[i];
        if (child != 0) {
            processOctreePostOrder(child, process);
        }
    }
    process(root);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualFaceNeighbor(OctreeNode<DataT> *node, int face)
{
    typedef OctreeTables OT;

    if (node == 0 || node->parent == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (OT::adj(face, octreeChildType(node))) {
        ancestor = getOctreeEqualFaceNeighbor(node->parent, face);
    } else {
        ancestor = node->parent;
    }

    return (ancestor != 0 ? ancestor->child[OT::reflect(face, octreeChildType(node))] : 0);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualEdgeNeighbor(OctreeNode<DataT> *node, int edge)
{
    typedef OctreeTables OT;

    if (node == 0 || node->parent == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (OT::adj(edge, octreeChildType(node))) {
        ancestor = getOctreeEqualEdgeNeighbor(node->parent, edge);
    } else if (OT::commonFace(edge, octreeChildType(node)) != OT::undef()) {
        ancestor =
            getOctreeEqualFaceNeighbor(node->parent, OT::commonFace(edge, octreeChildType(node)));
    } else {
        ancestor = node->parent;
    }

    return (ancestor != 0 ? ancestor->child[OT::reflect(edge, octreeChildType(node))] : 0);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualVertexNeighbor(OctreeNode<DataT> *node, int vertex)
{
    typedef OctreeTables OT;

    if (node == 0 || node->parent == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (OT::adj(vertex, octreeChildType(node))) {
        ancestor = getOctreeEqualVertexNeighbor(node->parent, vertex);
    } else if (OT::commonEdge(vertex, octreeChildType(node)) != OT::undef()) {
        ancestor =
            getOctreeEqualEdgeNeighbor(node->parent, OT::commonEdge(vertex, octreeChildType(node)));
    } else if (OT::commonFace(vertex, octreeChildType(node)) != OT::undef()) {
        ancestor =
            getOctreeEqualFaceNeighbor(node->parent, OT::commonFace(vertex, octreeChildType(node)));
    } else {
        ancestor = node->parent;
    }

    return (ancestor != 0 ? ancestor->child[OT::reflect(vertex, octreeChildType(node))] : 0);
}

template <typename DataT>
std::vector<OctreeNode<DataT> *> gatherOctreeEqualNeighbors(OctreeNode<DataT> *node)
{
    typedef OctreeTables OT;

    int numNeighbors = 26;
    std::vector<OctreeNode<DataT> *> neighbors(numNeighbors, (OctreeNode<DataT> *)0);

    if (node == 0 || node->parent == 0) {
        return (neighbors);
    }

    neighbors[OT::surface0N()] = getOctreeEqualFaceNeighbor(node, OT::surface0N());
    neighbors[OT::surface0P()] = getOctreeEqualFaceNeighbor(node, OT::surface0P());
    neighbors[OT::surface1N()] = getOctreeEqualFaceNeighbor(node, OT::surface1N());
    neighbors[OT::surface1P()] = getOctreeEqualFaceNeighbor(node, OT::surface1P());
    neighbors[OT::surface2N()] = getOctreeEqualFaceNeighbor(node, OT::surface2N());
    neighbors[OT::surface2P()] = getOctreeEqualFaceNeighbor(node, OT::surface2P());

    neighbors[OT::edge0NN()] = getOctreeEqualEdgeNeighbor(node, OT::edge0NN());
    neighbors[OT::edge0NP()] = getOctreeEqualEdgeNeighbor(node, OT::edge0NP());
    neighbors[OT::edge0PN()] = getOctreeEqualEdgeNeighbor(node, OT::edge0PN());
    neighbors[OT::edge0PP()] = getOctreeEqualEdgeNeighbor(node, OT::edge0PP());
    neighbors[OT::edge1NN()] = getOctreeEqualEdgeNeighbor(node, OT::edge1NN());
    neighbors[OT::edge1NP()] = getOctreeEqualEdgeNeighbor(node, OT::edge1NP());
    neighbors[OT::edge1PN()] = getOctreeEqualEdgeNeighbor(node, OT::edge1PN());
    neighbors[OT::edge1PP()] = getOctreeEqualEdgeNeighbor(node, OT::edge1PP());
    neighbors[OT::edge2NN()] = getOctreeEqualEdgeNeighbor(node, OT::edge2NN());
    neighbors[OT::edge2NP()] = getOctreeEqualEdgeNeighbor(node, OT::edge2NP());
    neighbors[OT::edge2PN()] = getOctreeEqualEdgeNeighbor(node, OT::edge2PN());
    neighbors[OT::edge2PP()] = getOctreeEqualEdgeNeighbor(node, OT::edge2PP());

    neighbors[OT::cornerNNN()] = getOctreeEqualVertexNeighbor(node, OT::cornerNNN());
    neighbors[OT::cornerNNP()] = getOctreeEqualVertexNeighbor(node, OT::cornerNNP());
    neighbors[OT::cornerNPN()] = getOctreeEqualVertexNeighbor(node, OT::cornerNPN());
    neighbors[OT::cornerNPP()] = getOctreeEqualVertexNeighbor(node, OT::cornerNPP());
    neighbors[OT::cornerPNN()] = getOctreeEqualVertexNeighbor(node, OT::cornerPNN());
    neighbors[OT::cornerPNP()] = getOctreeEqualVertexNeighbor(node, OT::cornerPNP());
    neighbors[OT::cornerPPN()] = getOctreeEqualVertexNeighbor(node, OT::cornerPPN());
    neighbors[OT::cornerPPP()] = getOctreeEqualVertexNeighbor(node, OT::cornerPPP());

    return (neighbors);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualFaceNeighbor(OctreeNode<DataT> *node, int face)
{
    typedef OctreeTables OT;

    if (node == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (node->parent != 0 && OT::adj(face, octreeChildType(node))) {
        ancestor = getOctreeGreaterEqualFaceNeighbor(node->parent, face);
    } else {
        ancestor = node->parent;
    }

    return (
        ancestor != 0 && !ancestor->isLeaf
            ? ancestor->child[OT::reflect(face, octreeChildType(node))]
            : ancestor);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualEdgeNeighbor(OctreeNode<DataT> *node, int edge)
{
    typedef OctreeTables OT;

    if (node == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (node->parent == 0) {
        ancestor = 0;
    } else if (OT::adj(edge, octreeChildType(node))) {
        ancestor = getOctreeGreaterEqualEdgeNeighbor(node->parent, edge);
    } else if (OT::commonFace(edge, octreeChildType(node)) != OT::undef()) {
        ancestor = getOctreeGreaterEqualFaceNeighbor(
            node->parent, OT::commonFace(edge, octreeChildType(node)));
    } else {
        ancestor = node->parent;
    }

    return (
        ancestor != 0 && !ancestor->isLeaf
            ? ancestor->child[OT::reflect(edge, octreeChildType(node))]
            : ancestor);
}

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualVertexNeighbor(OctreeNode<DataT> *node, int vertex)
{
    typedef OctreeTables OT;

    if (node == 0) {
        return (0);
    }

    OctreeNode<DataT> *ancestor = 0;
    if (node->parent == 0) {
        ancestor = 0;
    } else if (OT::adj(vertex, octreeChildType(node))) {
        ancestor = getOctreeGreaterEqualVertexNeighbor(node->parent, vertex);
    } else if (OT::commonEdge(vertex, octreeChildType(node)) != OT::undef()) {
        ancestor = getOctreeGreaterEqualEdgeNeighbor(
            node->parent, OT::commonEdge(vertex, octreeChildType(node)));
    } else if (OT::commonFace(vertex, octreeChildType(node)) != OT::undef()) {
        ancestor = getOctreeGreaterEqualFaceNeighbor(
            node->parent, OT::commonFace(vertex, octreeChildType(node)));
    } else {
        ancestor = node->parent;
    }

    return (
        ancestor != 0 && !ancestor->isLeaf
            ? ancestor->child[OT::reflect(vertex, octreeChildType(node))]
            : ancestor);
}

template <typename DataT>
std::vector<OctreeNode<DataT> *> gatherOctreeGreaterEqualNeighbors(OctreeNode<DataT> *node)
{
    typedef OctreeTables OT;

    int numNeighbors = 26;
    std::vector<OctreeNode<DataT> *> neighbors(numNeighbors, (OctreeNode<DataT> *)0);

    if (node == 0 || node->parent == 0) {
        return (neighbors);
    }

    neighbors[OT::surface0N()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface0N());
    neighbors[OT::surface0P()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface0P());
    neighbors[OT::surface1N()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface1N());
    neighbors[OT::surface1P()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface1P());
    neighbors[OT::surface2N()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface2N());
    neighbors[OT::surface2P()] = getOctreeGreaterEqualFaceNeighbor(node, OT::surface2P());

    neighbors[OT::edge0NN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge0NN());
    neighbors[OT::edge0NP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge0NP());
    neighbors[OT::edge0PN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge0PN());
    neighbors[OT::edge0PP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge0PP());
    neighbors[OT::edge1NN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge1NN());
    neighbors[OT::edge1NP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge1NP());
    neighbors[OT::edge1PN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge1PN());
    neighbors[OT::edge1PP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge1PP());
    neighbors[OT::edge2NN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge2NN());
    neighbors[OT::edge2NP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge2NP());
    neighbors[OT::edge2PN()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge2PN());
    neighbors[OT::edge2PP()] = getOctreeGreaterEqualEdgeNeighbor(node, OT::edge2PP());

    neighbors[OT::cornerNNN()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerNNN());
    neighbors[OT::cornerNNP()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerNNP());
    neighbors[OT::cornerNPN()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerNPN());
    neighbors[OT::cornerNPP()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerNPP());
    neighbors[OT::cornerPNN()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerPNN());
    neighbors[OT::cornerPNP()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerPNP());
    neighbors[OT::cornerPPN()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerPPN());
    neighbors[OT::cornerPPP()] = getOctreeGreaterEqualVertexNeighbor(node, OT::cornerPPP());

    return (neighbors);
}

template <typename DataT>
OctreePeriodicExtension<DataT>::OctreePeriodicExtension(
    OctreeNode<DataT> *originalRoot_, bool xPeriodic_, bool yPeriodic_, bool zPeriodic_) :
    originalRoot(originalRoot_),
    periodicExtensionRoot(0),
    parentOfOriginalRoot(0),
    levelOfOriginalRoot(0)
{
    if (originalRoot && (xPeriodic_ || yPeriodic_ || zPeriodic_)) {
        parentOfOriginalRoot = originalRoot->parent;
        levelOfOriginalRoot = originalRoot->level;
        recalibrateOriginalRootLevels(0);

        // Levels -1 and -2 are reserved for the octree periodic extensions only!
        // At level -2, only the child 0 is allocated, and at level -1 only the child 7 is
        // allocated. If any of these conventions is changed, the implementation of the
        // octreeChildType() function must also be updated.

        periodicExtensionRoot = new OctreeNode<DataT>(0, 0);
        periodicExtensionRoot->level = -2;

        periodicExtensionRoot->child[0] = new OctreeNode<DataT>(periodicExtensionRoot, 0);
        periodicExtensionRoot->isLeaf = false;

        originalRoot->parent = periodicExtensionRoot->child[0];
        periodicExtensionRoot->child[0]->child[7] = originalRoot;
        periodicExtensionRoot->child[0]->isLeaf = false;

        if (xPeriodic_) {
            periodicExtensionRoot->child[0]->child[3] = periodicExtensionRoot->child[0]->child[7];
            periodicExtensionRoot->child[4] = periodicExtensionRoot->child[0];
        }
        if (yPeriodic_) {
            periodicExtensionRoot->child[0]->child[5] = periodicExtensionRoot->child[0]->child[7];
            periodicExtensionRoot->child[0]->child[1] = periodicExtensionRoot->child[0]->child[3];
            periodicExtensionRoot->child[2] = periodicExtensionRoot->child[0];
            periodicExtensionRoot->child[6] = periodicExtensionRoot->child[4];
        }
        if (zPeriodic_) {
            periodicExtensionRoot->child[0]->child[6] = periodicExtensionRoot->child[0]->child[7];
            periodicExtensionRoot->child[0]->child[2] = periodicExtensionRoot->child[0]->child[3];
            periodicExtensionRoot->child[0]->child[4] = periodicExtensionRoot->child[0]->child[5];
            periodicExtensionRoot->child[0]->child[0] = periodicExtensionRoot->child[0]->child[1];
            periodicExtensionRoot->child[1] = periodicExtensionRoot->child[0];
            periodicExtensionRoot->child[5] = periodicExtensionRoot->child[4];
            periodicExtensionRoot->child[3] = periodicExtensionRoot->child[2];
            periodicExtensionRoot->child[7] = periodicExtensionRoot->child[6];
        }
    }
}

template <typename DataT>
OctreePeriodicExtension<DataT>::~OctreePeriodicExtension()
{
    restoreOriginalRootAndDestroyPeriodicExtension();
}

template <typename DataT>
OctreeNode<DataT> *OctreePeriodicExtension<DataT>::get() const
{
    return (originalRoot);
}

template <typename DataT>
OctreeNode<DataT> *OctreePeriodicExtension<DataT>::release()
{
    restoreOriginalRootAndDestroyPeriodicExtension();
    return (originalRoot);
}

template <typename DataT>
void OctreePeriodicExtension<DataT>::recalibrateOriginalRootLevels(int newRootLevel) const
{
    if (originalRoot->level == newRootLevel) {
        return;
    }

    originalRoot->parent = 0;

    struct RecalibrateLevels {
        RecalibrateLevels(int rootLevel_) : rootLevel(rootLevel_) { }
        void operator()(OctreeNode<DataT> *node)
        {
            if (node->parent) {
                node->level = node->parent->level + 1;
            } else {
                node->level = rootLevel;
            }
        }
        int rootLevel;
    };

    RecalibrateLevels recalibrateLevels(newRootLevel);
    processOctreePreOrder(originalRoot, recalibrateLevels);
}

template <typename DataT>
void OctreePeriodicExtension<DataT>::restoreOriginalRootAndDestroyPeriodicExtension()
{
    if (periodicExtensionRoot) {
        for (int i = 0; i < 8; i++) {
            periodicExtensionRoot->child[0]->child[i] = 0;
        }
        periodicExtensionRoot->child[0]->isLeaf = true;

        delete periodicExtensionRoot->child[0];
        delete periodicExtensionRoot;
        periodicExtensionRoot = 0;

        // Restore originalRoot to its state before the periodic extension.
        recalibrateOriginalRootLevels(levelOfOriginalRoot);
        originalRoot->parent = parentOfOriginalRoot;
    }
}

}  // namespace plb

#endif  // OCTREE_HH
