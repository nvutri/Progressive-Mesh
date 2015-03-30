#ifndef _EDGE_H_
#define _EDGE_H_

#include <string>

class Halfedge;

class Edge {
public:


    Edge() {
        m_halfedge[0] = NULL;
        m_halfedge[1] = NULL;
        m_propertyIndex = -1;
    }

    Edge(Halfedge *he0, Halfedge *he1) {
        m_halfedge[0] = he0;
        m_halfedge[1] = he1;
        m_propertyIndex = -1;
    }

    ~Edge() {
        ;
    }

    // Pointers for Halfedge Data Structure
    Halfedge *&he(int i) {
        return m_halfedge[i];
    }

    Halfedge *&twin(Halfedge *he) {
        return (he == m_halfedge[0]) ? (m_halfedge[1]) : (m_halfedge[0]);
    }

    // Computed by Halfedge Data Structure
    bool boundary() {
        return (!m_halfedge[0] || !m_halfedge[1]);
    }

    //optional
    int &index() {
        return m_propertyIndex;
    }

    std::string &PropertyStr() {
        return m_propertyStr;
    }

protected:
    //for Halfedge Data Structure
    Halfedge *m_halfedge[2];        // for boundary edge, m_halfedge[1]=NULL

    //optional
    std::string m_propertyStr;        // a string to store edge properties
    int m_propertyIndex;    // index to Property array
};

#endif
