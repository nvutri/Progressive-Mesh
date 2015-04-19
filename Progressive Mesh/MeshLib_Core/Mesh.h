#ifndef _MESH_H_
#define _MESH_H_

#include "Edge.h"
#include "Face.h"
#include "Halfedge.h"
#include "Vertex.h"
#include "Point.h"

#if defined __GNUC__ || defined __APPLE__

#include <ext/hash_map>

#else
#include <hash_map>
#endif

#include <list>
#include <vector>

class Mesh {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////								Methods										//////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:
    //(1) Constructor and Destructor
    Mesh();

    //Mesh::Mesh(std::vector<Face*>& faces);
    ~Mesh();

    //(2) I/O
    int numVertices() {
        return m_verts.size();
    }                            //number of vertices
    int numEdges() {
        return m_edges.size();
    }                            //number of edges
    int numFaces() {
        return m_faces.size();
    }                            //number of faces
    void copyTo(Mesh &targetMesh);                                        //copy current mesh to the target mesh
    bool readMFile(const char inFile[]);                                    //read a mesh from inFile
    bool writeMFile(const char outFile[]);                                    //write a mesh to outFile
    void clear();

    //(3) BASIC OPERATIONS
    //Check whether an element is on the boundary:
    bool isBoundary(Vertex *v) const {
        return v->boundary();
    }

    bool isBoundary(Edge *e) const {
        return e->boundary();
    }

    bool isBoundary(Halfedge *he) const {
        return !(he->twin());
    }

    // Indexing elements
    Vertex *indVertex(unsigned int ind) {
        return (ind >= m_verts.size() ? NULL : m_verts[ind]);
    }

    Face *indFace(unsigned int ind) {
        return (ind >= m_faces.size() ? NULL : m_faces[ind]);
    }

    Edge *indEdge(unsigned int ind) {
        return (ind >= m_edges.size() ? NULL : m_edges[ind]);
    }

    Edge *vertexEdge(Vertex *v0, Vertex *v1);            //To find an edge sharing v0 and v1
    Halfedge *vertexHalfedge(Vertex *srcV, Vertex *trgV);    //To find a half-edge from v0 to v1
    Edge *idEdge(int vid0, int vid1);

    Halfedge *idHalfedge(int srcVid, int trgVid);


public:
    //(4) Mesh Modification Operators
    Vertex *createVertex(int vertexId);

    Face *createFace(int faceId);

    Face *createFace(int faceId, Vertex *verts[]);

    Face *createFace(int faceId, int vIds[]);

    Edge *createEdge();

    Edge *createEdge(Halfedge *he0, Halfedge *he1);


    void LabelBoundaryVertices();

    ////(6) Some Basic Mesh Processing Operations
    ////	Flip the input edge, and return it.
    ////	if the flip works, return the resultant edge, otherwise return NULL
    //Edge *				EdgeFlip(Edge * e);

    //// Split the given edge(input) into two.
    ////	  return :  newly created vertex or NULL(when the split fails);
    //Vertex *			EdgeSplit(Edge * e);
    //Vertex *			FaceSplit(Face * f, double bary[3]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////								Variables										//////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public:
    //stdext::hash_map<Vertex *, int> vertIndex;
    std::vector<Edge *> m_edges;        // edge container
    std::vector<Vertex *> m_verts;        // vertex container
    std::vector<Face *> m_faces;        // face container

    //a temporary container to store halfedges surrounding some vertices
    std::vector<std::vector<Halfedge *>> v_adjInHEList;

protected:
    friend class MeshVertexIterator;

    friend class MeshEdgeIterator;

    friend class MeshFaceIterator;

    friend class MeshHalfedgeIterator;

    friend class MeshUtility;

    friend class MeshIO;
};


#endif //_XMESHLIB_CORE_MESH_H_ defined
