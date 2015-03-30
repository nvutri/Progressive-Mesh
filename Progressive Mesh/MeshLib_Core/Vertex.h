#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "Halfedge.h"
#include "Point.h"

class Vertex {
public:
    Vertex() {
        m_halfedge = NULL;
        m_boundary = false;
        m_propertyIndex = -1;
    }

    ~Vertex() {
        ;
    }

    Point &point() {
        return m_point;
    }

    //Pointers for Halfedge Data Structure
    Halfedge *&he() {
        return m_halfedge;
    }

    //Computed by Halfedge Data Structure
    bool &boundary() {
        return m_boundary;
    }    //whether this is a boundary vertex

    //Rotation operations
    Halfedge *most_ccw_in_halfedge();

    Halfedge *most_ccw_out_halfedge();

    Halfedge *most_clw_in_halfedge();

    Halfedge *most_clw_out_halfedge();

    //optional
    int &index() {
        return m_propertyIndex;
    }

    std::string &PropertyStr() {
        return m_propertyStr;
    }


protected:
    //for Halfedge Data Structure
    Point m_point;
    Halfedge *m_halfedge;

    //optional
    bool m_boundary; // whether this is a boundary vertex
    std::string m_propertyStr;
    int m_propertyIndex; // index to Property array
};

#endif

