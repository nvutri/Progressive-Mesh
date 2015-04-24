#include "Mesh.h"
#include <iostream>
#include <sstream>

#pragma warning (disable : 4996)
#pragma warning (disable : 4018)

Mesh::Mesh() {
    ;
}

Mesh::~Mesh() {
    clear();
}

void Mesh::clear() {
    for (std::vector<Face *>::iterator fiter = m_faces.begin(); fiter != m_faces.end(); ++fiter) {
        Face *f = *fiter;
        delete f;
    }
    for (std::vector<Edge *>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++eiter) {
        Edge *e = *eiter;
        Halfedge *he1 = e->he(0);
        Halfedge *he2 = e->he(1);
        delete he1;
        if (he2) delete he2;
        delete e;
    }
    for (std::vector<Vertex *>::iterator viter = m_verts.begin(); viter != m_verts.end(); ++viter) {
        Vertex *v = *viter;
        delete v;
    }
    m_verts.clear();
    m_edges.clear();
    m_faces.clear();
}

Edge *Mesh::vertexEdge(Vertex *v0, Vertex *v1) {
    //First, check the right most side
    Halfedge *he0 = v0->most_clw_out_halfedge();
    if (he0->target() == v1)
        return he0->edge();

    Halfedge *he = he0->ccw_rotate_about_source();
    while (he != he0 && he) {
        if (he->target() == v1)
            return he->edge();
        he = he->ccw_rotate_about_source();
    }

    if (!v0->boundary())
        return NULL;  //for an interior vertex v0, the entire one-ring has been checked

    //Finally, check the left most side
    he = v0->most_ccw_in_halfedge();
    if (he->source() == v1)
        return he->edge();
    else
        return NULL;
}

Edge *Mesh::idEdge(int id0, int id1) {
    Vertex *v0 = indVertex(id0);
    Vertex *v1 = indVertex(id1);
    if (!v0 || !v1)
        return NULL;
    else
        return vertexEdge(v0, v1);
}

Halfedge *Mesh::vertexHalfedge(Vertex *v0, Vertex *v1) {
    Halfedge *he0 = v0->most_clw_out_halfedge();
    if (he0->target() == v1) return he0;
    Halfedge *he = he0->ccw_rotate_about_source();
    while (he != he0 && he) {
        if (he->target() == v1)
            return he;
        he = he->ccw_rotate_about_source();
    }
    return NULL;  //the entire one-ring of v0 has been checked
}


Halfedge *Mesh::idHalfedge(int srcVID, int trgVID) { //Check the surrounding outgoing half-edges from the source vertex
    Vertex *v0 = indVertex(srcVID);
    Vertex *v1 = indVertex(trgVID);
    if (!v0 || !v1)
        return NULL;
    else
        return vertexHalfedge(v0, v1);
}

// Create new geometric simplexes
Vertex *Mesh::createVertex(int vertexId) {
    Vertex *v = new Vertex;
    vertexId -= 1;
    v->index() = vertexId;
    if (vertexId >= m_verts.size()) {
        m_verts.resize(vertexId + 1);
    }
    m_verts[vertexId] = v;
    return v;
}

Face *Mesh::createFace(int faceId) {
    Face *f = new Face();
    faceId -= 1;
    f->index() = faceId;
    if (faceId >= m_faces.size()) {
        m_faces.resize(faceId + 1);
    }
    m_faces[faceId] = f;
    return f;
}

Edge *Mesh::createEdge() {
    Edge *e = new Edge();
    e->index() = m_edges.size();
    m_edges.push_back(e);
    return e;
}

Edge *Mesh::createEdge(Halfedge *he0, Halfedge *he1) {
    Edge *e = new Edge(he0, he1);
    e->index() = m_edges.size();
    m_edges.push_back(e);
    return e;
}


Face *Mesh::createFace(int faceId, int v[3]) {
    Vertex *verts[3];
    for (int i = 0; i < 3; i++) {
        verts[i] = indVertex(v[i]);
        if (!verts[i]) {
            std::cerr << "Error: invalid vertex id: " << v[i] << " provided when creating face "
                    << m_faces.size() << " !" << std::endl;
            return NULL;
        }
    }
    Face *f = createFace(faceId, verts);
    return f;
}

Face *Mesh::createFace(int faceId, Vertex *verts[3]) {
    int i;
    Face *f = createFace(faceId);
    //create Half-edges
    Halfedge *hes[3];
    for (i = 0; i < 3; i++) {
        hes[i] = new Halfedge();
        hes[i]->target() = verts[i];
        verts[i]->he() = hes[i];
        std::vector<Halfedge *> &adjInHEList = v_adjInHEList[verts[i]->index()];
        adjInHEList.push_back(hes[i]);
    }
    //linking to each other, and linking to the face
    for (i = 0; i < 3; i++) {
        hes[i]->next() = hes[(i + 1) % 3];
        hes[i]->prev() = hes[(i + 2) % 3];
        hes[i]->face() = f;
        f->he() = hes[i];
    }
    //Linking these halfedges with edges
    for (i = 0; i < 3; i++) {
        Vertex *ev0 = verts[i]; // target
        Vertex *ev1 = verts[(i + 2) % 3]; // source
        Edge *e = NULL;
        // The new halfedge is [ev1, ev0]:
        // should we generate a new edge with these two vertices?
        // not necessary: if [ev0, ev1] was an existing halfedge
        // So, we shall check whether halfedge [ev0,ev1] was created and added into the mesh before,
        // We had the container of one-ring incoming half-edges of ev1
        std::vector<Halfedge *> &adjInHEList = v_adjInHEList[ev1->index()];
        for (int j = 0; j < adjInHEList.size(); ++j) {
            Halfedge *he = adjInHEList[j];
            if (he->source() == ev0) {
                e = he->edge();
                break;
            }
        }
        if (e) {
            // [ev0, ev1] exists hence the edge exists
            if (!(e->he(1)))
                e->he(1) = hes[i];
            else {// Because when we create an edge the first time, the halfedge becomes its he[0], and he[1] should=NULL;
                std::cerr << "Error: an edge appears more than twice. Non-manifold surfaces!" << std::endl;
                exit(0);
            }
        }
        else // [ev0, ev1] was not found, create this new edge
            e = createEdge(hes[i], NULL);
        hes[i]->edge() = e;
    }
    return f;
}

void Mesh::LabelBoundaryVertices() {// we should do this once, after the half-edge data structure has been created
    for (std::vector<Edge *>::iterator eiter = m_edges.begin(); eiter != m_edges.end(); ++eiter) {
        Edge *edge = *eiter;
        Halfedge *he[2];
        he[0] = edge->he(0);
        he[1] = edge->he(1);
        if (!he[1]) {
            he[0]->target()->boundary() = true;
            he[0]->source()->boundary() = true;
        }
    }
}

bool Mesh::readMFile(const char inputFile[]) {
    std::cout << "Reading mesh " << inputFile << " ...";
    FILE *fp = fopen(inputFile, "r");
    if (!fp) {
        std::cerr << "Can't open file " << inputFile << "!" << std::endl;
        return false;
    }

    char line[1024];
    clear();

    while (!feof(fp)) {
        fgets(line, 1024, fp);
        if (!strlen(line)) continue;
        char *str = strtok(line, " \r\n");
        if (!str) continue;

        if (!strcmp(str, "Vertex")) { //parsing a line of vertex element
            str = strtok(NULL, " \r\n{");
            int vertexId = atoi(str);
            Vertex *v = createVertex(vertexId); // Note: the index in M File starts with 1, while our default index starts with 0
            // Parsing (x,y,z)
            for (int i = 0; i < 3; i++) {
                str = strtok(NULL, " \r\n{");
                v->point()[i] = atof(str);
            }

            // Storing the neighboring halfedges, for efficient edge searching in the initialization stage
            std::vector<Halfedge *> heList;
            if (v->index() >= v_adjInHEList.size()) {
                v_adjInHEList.resize(v->index() + 1);
            }
            v_adjInHEList[v->index()] = heList;

            // Parsing the property string
            str = strtok(NULL, "\r\n");
            if (!str || strlen(str) == 0) continue;
            std::string s(str);
            int sp = s.find("{");
            int ep = s.find("}");
            if (sp >= 0 && ep >= 0)
                v->PropertyStr() = s.substr(sp + 1, ep - sp - 1);
            continue;
        }
        else if (!strcmp(str, "Face")) { //parsing a line of face element

            str = strtok(NULL, " \r\n");
            if (!str || strlen(str) == 0) continue;
            int faceId = atoi(str);
            int vids[3];
            for (int i = 0; i < 3; i++) {
                str = strtok(NULL, " \r\n{");
                vids[i] = atoi(str) - 1; //Note: the index in M File starts with 1, while our default index starts with 0
            }
            Face *f = createFace(faceId, vids);
            if (!str) continue;
            str = strtok(NULL, "\r\n");
            if (!str || strlen(str) == 0) continue;
            std::string s(str);
            int sp = s.find("{");
            int ep = s.find("}");
            if (sp >= 0 && ep >= 0)
                f->PropertyStr() = s.substr(sp + 1, ep - sp - 1);
            continue;
        }
    }

    LabelBoundaryVertices();

    fclose(fp);

    int heInd = 0;
    for (std::vector<Edge *>::iterator eit = m_edges.begin(); eit != m_edges.end(); ++eit) {
        Edge *e = *eit;
        Halfedge *he0 = e->he(0);
        he0->index() = heInd++;
        Halfedge *he1 = e->he(1);
        if (he1)
            he1->index() = heInd++;
    }

    // After initialization, the neighboring halfedges list is not needed anymore; remove them to save space
    v_adjInHEList.resize(m_verts.size());
    for (int i = 0; i < v_adjInHEList.size(); ++i)
        v_adjInHEList[i].clear();
    v_adjInHEList.clear();

    printf("Done!\n");
    return true;
}


bool Mesh::writeMFile(const char outputFile[]) {
    FILE *fp = fopen(outputFile, "w");
    if (!fp) {
        std::cerr << "Cannot open file " << outputFile << "to write!" << std::endl;
        return false;
    }

    std::cout << "Writing mesh " << outputFile << " ...";
    std::vector<Vertex *>::iterator vit;
    for (vit = m_verts.begin(); vit != m_verts.end(); ++vit) {
        Vertex *ver = *vit;
        std::ostringstream oss;
        oss.precision(6);
        oss.setf(std::ios::fixed, std::ios::floatfield);  //setting the output precision: now 6
        oss << "Vertex " << ver->index() + 1 << " " << ver->point()[0] << " " << ver->point()[1] << " " << ver->point()[2];
        fprintf(fp, "%s ", oss.str().c_str());
        if (ver->PropertyStr().size() > 0)
            fprintf(fp, "{%s}", ver->PropertyStr().c_str());
        fprintf(fp, "\n");
    }

    std::vector<Face *>::iterator fit;
    for (fit = m_faces.begin(); fit != m_faces.end(); ++fit) {
        Face *face = *fit;
        Halfedge *the0 = face->he();
        Halfedge *the1 = the0->next();
        int v0 = the0->source()->index() + 1;
        int v1 = the0->target()->index() + 1;
        int v2 = the1->target()->index() + 1;
        fprintf(fp, "Face %d %d %d %d ", face->index() + 1, v0, v1, v2);
        if (face->PropertyStr().size() > 0)
            fprintf(fp, "{%s}", face->PropertyStr().c_str());
        fprintf(fp, "\n");
    }

    fclose(fp);
    std::cout << "Done!" << std::endl;
    return true;
}

void Mesh::copyTo(Mesh &tMesh) {
    std::cout << "Copying the mesh...";

    for (std::vector<Vertex *>::iterator viter = m_verts.begin();
         viter != m_verts.end(); ++viter) {
        Vertex *v = *viter;
        Vertex *nv = tMesh.createVertex(v->index());
        nv->point() = v->point();
        nv->PropertyStr() = v->PropertyStr();
        nv->boundary() = v->boundary();
    }

    std::vector<Face *>::iterator fiter = m_faces.begin();
    for (; fiter != m_faces.end(); ++fiter) {
        Face *f = *fiter;
        Face *nf = tMesh.createFace(f->index());
        Halfedge *he[3];
        Halfedge *nhe[3];
        he[0] = f->he();
        he[1] = he[0]->next();
        he[2] = he[1]->next();
        for (int j = 0; j < 3; ++j) {
            nhe[j] = new Halfedge();
            Vertex *v1 = he[j]->target();
            Vertex *nv1 = tMesh.indVertex(v1->index());
            nv1->he() = nhe[j];
            nhe[j]->target() = nv1;
            nhe[j]->face() = nf;
        }
        for (int j = 0; j < 3; ++j) {
            nhe[j]->next() = nhe[(j + 1) % 3];
            nhe[j]->prev() = nhe[(j + 2) % 3];
        }
        nf->he() = nhe[0];
        nf->PropertyStr() = f->PropertyStr();
    }

    for (fiter = m_faces.begin(); fiter != m_faces.end(); ++fiter) {
        Face *f = *fiter;
        Face *nf = tMesh.indFace(f->index());
        Halfedge *he[3];
        Halfedge *nhe[3];
        he[0] = f->he();
        he[1] = he[0]->next();
        he[2] = he[1]->next();
        nhe[0] = nf->he();
        nhe[1] = nhe[0]->next();
        nhe[2] = nhe[1]->next();
        for (int i = 0; i < 3; ++i) {
            Edge *e = he[i]->edge();
            if (he[i] == e->he(0)) {
                if (e->boundary()) {
                    Edge *ne = tMesh.createEdge(nhe[i], NULL);
                    nhe[i]->edge() = ne;
                    ne->PropertyStr() = e->PropertyStr();
                    continue;
                }
                //get its twin: twin_he
                Halfedge *twin_he, *twin_nhe;
                Face *tf = he[i]->twin()->face();
                Face *tnf = tMesh.indFace(tf->index());
                twin_he = tf->he();
                twin_nhe = tnf->he();
                for (int j = 0; j < 3; ++j) {
                    if (twin_he->edge() == e) break;
                    twin_he = twin_he->next();
                    twin_nhe = twin_nhe->next();
                }
                Edge *ne = tMesh.createEdge(nhe[i], twin_nhe);
                nhe[i]->edge() = ne;
                twin_nhe->edge() = ne;
                ne->PropertyStr() = e->PropertyStr();
            }
        }
    }
    std::cout << "Done!" << std::endl;
}
