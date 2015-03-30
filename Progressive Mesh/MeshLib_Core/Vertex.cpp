#include "Vertex.h"

Halfedge *Vertex::most_ccw_in_halfedge() {
    if (!m_boundary)
        return m_halfedge; //for interior vertex, randomly pick one, any halfedge can be most ccw

    Halfedge *he = m_halfedge;
    Halfedge *nhe = he->ccw_rotate_about_target();
    Halfedge *startHe = he;
    while (nhe) {
        he = nhe;
        nhe = nhe->ccw_rotate_about_target();
        if (he == startHe) return startHe;  //This should not happen when the mesh is a valid manifold. This check avoids endless loop in handling nonmanifold
    }
    return he;
}

Halfedge *Vertex::most_clw_in_halfedge() {
    if (!m_boundary)
        return m_halfedge;
    Halfedge *he = m_halfedge;
    Halfedge *nhe = he->clw_rotate_about_target();
    Halfedge *startHe = he;
    while (nhe) {
        he = nhe;
        nhe = nhe->clw_rotate_about_target();
        if (he == startHe) return startHe;  //This should not happen when the mesh is a valid manifold. This check avoids endless loop in handling nonmanifold
    }
    return he;
}

Halfedge *Vertex::most_ccw_out_halfedge() {
    if (!m_boundary)
        return m_halfedge->twin();
    Halfedge *he = m_halfedge->twin();
    if (!he)
        he = m_halfedge->next();
    Halfedge *startHe = he;
    Halfedge *nhe = he->ccw_rotate_about_source();
    while (nhe) {
        he = nhe;
        nhe = nhe->ccw_rotate_about_source();
        if (he == startHe) return startHe;  //This should not happen when the mesh is a valid manifold. This check avoids endless loop in handling nonmanifold
    }
    return he;
}

Halfedge *Vertex::most_clw_out_halfedge() {
    if (!m_boundary)
        return m_halfedge->twin();
    Halfedge *he = m_halfedge->twin();
    if (!he)
        he = m_halfedge->next();
    Halfedge *startHe = he;
    Halfedge *nhe = he->clw_rotate_about_source();
    while (nhe) {
        he = nhe;
        nhe = nhe->clw_rotate_about_source();
        if (he == startHe) return startHe;  //This should not happen when the mesh is a valid manifold. This check avoids endless loop in handling nonmanifold
    }
    return he;
}