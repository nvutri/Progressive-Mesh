#pragma once


#ifdef __APPLE__

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#include "MeshLib_Core/Mesh.h"
#include "MeshLib_Core/Iterators.h"

#else
#include <GL/glut.h>
#include "MeshLib_Core/Mesh.h"
#include "MeshLib_Core/Iterators.h"
#endif

#include "VSplitRecord.h"
#include <iostream>
#include <stdlib.h>
#include <pthread.h>
#include <set>

namespace XMeshLib {
    class PM { //Currently implemented for closed meshes only
    public:
        PM(Mesh *cMesh) {
            tMesh = cMesh;
            tmpMesh = NULL;
            currentMeshResolution = cMesh->numVertices();
        }

        void SetEdgePriority();

        void ProcessCoarsening(int targetVertSize = 4);

        void ProcessRefinement(int targetVertSize = -1); // by default, refine back to full resolution

        static void *FindAndCollapseEdge(void *);

        void GetNextCollapseEdges();

        Edge *GetNextCollapseEdge();

        bool CheckEdgeCollapseCondition(Edge *e);

        Vertex *EdgeCollapse(Edge *e, VSplitRecord &vsRec); //Currently for closed surfaces only
        int VertexSplit(VSplitRecord &vsRec);

        bool WriteVsplitRecord(const char filename[], std::vector<VSplitRecord> &vsRecList);

        bool SaveMesh(const char filename[]);

        // The temporal mesh stored in tMesh cannot be directly written or copied,
        // due to the unremoved isolated vertices (I keep them on purpose).
        // This following function converts it to a valid mesh and stores it in tmpMesh;
        void GetValidTmpMesh();

    public:
        Mesh *tMesh;
        Mesh *tmpMesh;
        int baseMeshResolution;
        int currentMeshResolution;
        std::vector<VSplitRecord> vsRecList;
        std::vector<int> tmpInd2OInd;
        // std::vector<int> oInd2TmpInd;
    protected:
        Vertex *CreateVertex();

        Face *CreateFace(int faceIndex);

        Edge *CreateEdge(int edgeIndex, Halfedge *he0, Halfedge *he1);

        void DeleteVertex(Vertex *v);

        void DeleteFace(Face *f);

        void DeleteEdge(Edge *e);

    protected: //for debugging
        void printHE(Halfedge *he);

        void printE(Edge *e);
    };
}
