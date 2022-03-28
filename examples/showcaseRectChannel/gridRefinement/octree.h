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

#ifndef OCTREE_H
#define OCTREE_H

#include <vector>

namespace plb {

template <typename DataT>
struct OctreeNode {
    OctreeNode(OctreeNode<DataT> *parent_, DataT *data_);

    DataT *data;
    int level;
    bool isLeaf;
    OctreeNode<DataT> *parent;
    OctreeNode<DataT> *child[8];
};

struct OctreeTables {
    static int surface0N();
    static int surface0P();
    static int surface1N();
    static int surface1P();
    static int surface2N();
    static int surface2P();

    static int edge0NN();
    static int edge0NP();
    static int edge0PN();
    static int edge0PP();
    static int edge1NN();
    static int edge1NP();
    static int edge1PN();
    static int edge1PP();
    static int edge2NN();
    static int edge2NP();
    static int edge2PN();
    static int edge2PP();

    static int cornerNNN();
    static int cornerNNP();
    static int cornerNPN();
    static int cornerNPP();
    static int cornerPNN();
    static int cornerPNP();
    static int cornerPPN();
    static int cornerPPP();

    static int border();
    static int smaller();

    static int adj(int dir, int oct);
    static int reflect(int dir, int oct);
    static int commonFace(int dir, int oct);
    static int commonEdge(int dir, int oct);
    static int undef();

    static const int adjTab[26][8];
    static const int reflTab[26][8];
    static const int comFaceTab[26][8];
    static const int comEdgeTab[26][8];
};

template <typename DataT>
int octreeChildType(OctreeNode<DataT> *node);

template <typename DataT>
void freeOctree(OctreeNode<DataT> *root);

template <typename DataT, class AssignData>
void splitOctreeNode(OctreeNode<DataT> *node, AssignData &assignData);

template <typename DataT>
void getMinMaxOctreeLeafNodeLevels(OctreeNode<DataT> *root, int &minLevel, int &maxLevel);

template <typename DataT, class Process>
void processOctreePreOrder(OctreeNode<DataT> *root, Process &process);

template <typename DataT, class Process>
void processOctreePostOrder(OctreeNode<DataT> *root, Process &process);

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualFaceNeighbor(OctreeNode<DataT> *node, int face);

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualEdgeNeighbor(OctreeNode<DataT> *node, int edge);

template <typename DataT>
OctreeNode<DataT> *getOctreeEqualVertexNeighbor(OctreeNode<DataT> *node, int vertex);

template <typename DataT>
std::vector<OctreeNode<DataT> *> gatherOctreeEqualNeighbors(OctreeNode<DataT> *node);

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualFaceNeighbor(OctreeNode<DataT> *node, int face);

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualEdgeNeighbor(OctreeNode<DataT> *node, int edge);

template <typename DataT>
OctreeNode<DataT> *getOctreeGreaterEqualVertexNeighbor(OctreeNode<DataT> *node, int vertex);

template <typename DataT>
std::vector<OctreeNode<DataT> *> gatherOctreeGreaterEqualNeighbors(OctreeNode<DataT> *node);

template <typename DataT>
class OctreePeriodicExtension {
public:
    // The OctreePeriodicExtension takes full charge of the originalRoot_ and all its subtrees.
    // This means that the user must not use originalRoot_ or any of its subtrees until it is
    // released from the OctreePeriodicExtension. However, the OctreePeriodicExtension does not
    // take over the memory management of the originalRoot_ tree. The user needs to resume
    // its memory management after it is released from the OctreePeriodicExtension.
    OctreePeriodicExtension(
        OctreeNode<DataT> *originalRoot_, bool xPeriodic_, bool yPeriodic_, bool zPeriodic_);
    ~OctreePeriodicExtension();
    OctreeNode<DataT> *get() const;
    OctreeNode<DataT> *release();

private:
    void recalibrateOriginalRootLevels(int newRootLevel) const;
    void restoreOriginalRootAndDestroyPeriodicExtension();

private:
    OctreeNode<DataT> *originalRoot;
    OctreeNode<DataT> *periodicExtensionRoot;
    OctreeNode<DataT> *parentOfOriginalRoot;
    int levelOfOriginalRoot;
};

}  // namespace plb

#endif  // OCTREE_H
