#include <thrust/execution_policy.h>
#include <thrust/remove.h>
#include <thrust/set_operations.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "PM.h"
#include <cassert>
#include <fstream>

#include "XLibCommon.h"

#define NUM_THREADS 16

pthread_mutex_t mutexLock;

using namespace XMeshLib;


thrust::host_vector<int> collapseEdges;
int *visibleEdges;
int *visibleEdgesEnd;
struct thread_data {
    int thread_id;
    PM *pm_ptr;
    Edge *edge;
    XMeshLib::VSplitRecord vsRec;
};

void PM::SetEdgePriority() {
}

bool PM::CheckEdgeCollapseCondition(Edge *e) {
    Halfedge *he = e->he(0);
    Vertex *v1 = he->source();
    Vertex *v2 = he->target();
    if (!(v1->visible && v2->visible))
      return false;
    std::vector<int> set1;
    std::vector<int> set2;
    for (VertexVertexIterator vvit(v1); !vvit.end(); ++vvit) {
        Vertex *vv1 = *vvit;
        if (vv1->visible) {
          set1.push_back(vv1->index());
        }
    }
    for (VertexVertexIterator vvit(v2); !vvit.end(); ++vvit) {
        Vertex *vv2 = *vvit;
        if (vv2->visible) {
          set2.push_back(vv2->index());
        }
    }
    std::sort(set1.begin(), set1.end());
    std::sort(set2.begin(), set2.end());
    std::vector<int> intSet(set1.size() + set2.size());
    std::vector<int>::iterator it = std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), intSet.begin());
    return (it - intSet.begin()) == 2;
}

Vertex *PM::EdgeCollapse(Edge *e, VSplitRecord &vsRec) {
    Halfedge *phe[6];
    phe[0] = e->he(0);
    phe[1] = e->he(1);
    assert(phe[1]);
    Face *f1 = phe[0]->face();
    Face *f2 = phe[1]->face();
    vsRec.face1 = f1->index();
    vsRec.face2 = f2->index();
    phe[2] = phe[0]->next();
    phe[4] = phe[0]->prev();
    phe[3] = phe[1]->prev();
    phe[5] = phe[1]->next();
    Halfedge *dhe2 = phe[2]->twin();
    Halfedge *dhe3 = phe[3]->twin();
    Halfedge *dhe4 = phe[4]->twin();
    Halfedge *dhe5 = phe[5]->twin();
    Vertex *vs = phe[0]->target();
    Vertex *vt = phe[1]->target();
    Vertex *vl = phe[2]->target();
    Vertex *vr = phe[5]->target();
    // merge vt to vs
    // (1) link all vt's neighboring halfedge to vs
    std::vector<Halfedge *> tmpHe;
    for (VertexInHalfedgeIterator hit(vt); !hit.end(); ++hit) {
        Halfedge *he = *hit;
        if (he == phe[1])
            continue;
        tmpHe.push_back(he);
    }
    for (unsigned int i = 0; i < tmpHe.size(); ++i)
        tmpHe[i]->target() = vs;

    // (2) set the new twin halfedges of dhe2 and dhe3 to dhe4 and dhe5, respectively
    Edge *e1 = dhe2->edge();
    Edge *e2 = dhe3->edge();
    Edge *de1 = dhe4->edge();
    Edge *de2 = dhe5->edge();
    e1->he(0) = dhe2;
    e1->he(1) = dhe4;
    dhe4->edge() = e1;
    e2->he(0) = dhe3;
    e2->he(1) = dhe5;
    dhe5->edge() = e2;

    // (3) set related vertices' halfedges
    vs->he() = dhe2;
    vl->he() = dhe4;
    vr->he() = dhe3;

    // (4) delete he[0~6], delete the two faces, and delete three edges
    DeleteFace(f1);
    DeleteFace(f2);

    vsRec.vs = vs;
    vsRec.vl = vl;
    vsRec.vr = vr;
    vsRec.vt = vt;
    vsRec.old_vs_pt = vs->point();
    vs->point() = (vs->point() + vt->point()) / 2;
    DeleteVertex(vt);
    vsRec.eIndex = e->index();
    vsRec.de1Index = de1->index();
    vsRec.de2Index = de2->index();
    DeleteEdge(e);
    DeleteEdge(de1);
    DeleteEdge(de2);
    return NULL;
}


int PM::VertexSplit(VSplitRecord &vsRec) {
    // (0) Obtain prestored primitives from vsplitRecord,
    // and undeleted primitives from halfedge data structure
    Vertex *vt = vsRec.vt;
    Vertex *vs = vsRec.vs;    // tMesh->indVertex(vsRec.vs_ind);
    Vertex *vl = vsRec.vl;    // tMesh->indVertex(vsRec.vl_ind);
    Vertex *vr = vsRec.vr;    // tMesh->indVertex(vsRec.vr_ind);
    Halfedge *dhe4 = tMesh->vertexHalfedge(vs, vl);
    Halfedge *dhe3 = tMesh->vertexHalfedge(vs, vr);
    if (!dhe3 || !dhe4) {
        return -1;
    }
    Halfedge *dhe2 = dhe4->twin();
    Halfedge *dhe5 = dhe3->twin();
    vs->point() = vsRec.old_vs_pt;

    // (2) Create new primitives
    assert(!tMesh->m_verts[vt->index()]->visible);
    tMesh->m_verts[vt->index()]->visible = true;
    Halfedge *phe[6];
    for (int i = 0; i < 6; ++i) {
        phe[i] = new Halfedge;
    }
    Edge *e_sl = dhe2->edge();
    Edge *e_sr = dhe3->edge();
    Edge *e_tl = CreateEdge(vsRec.eIndex, dhe4, phe[4]);
    Edge *e_tr = CreateEdge(vsRec.de1Index, dhe5, phe[5]);
    Edge *e_st = CreateEdge(vsRec.de2Index, phe[0], phe[1]);
    Face *f0 = CreateFace(vsRec.face1);
    Face *f1 = CreateFace(vsRec.face2);

    // (2) Collect following ccw order, about vt,
    // all in halfedges from dhe4->prev() to dhe5.
    // These he's target should be linked to vt
    Halfedge *che = dhe4->prev();
    std::vector<Halfedge *> cheList;
    while (che != dhe5) {
        cheList.push_back(che);
        che = che->ccw_rotate_about_target();
    };
    cheList.push_back(dhe5);

    // (3) Link primitives
    // a. between halfedges
    phe[0]->next() = phe[2];
    phe[0]->prev() = phe[4];
    phe[2]->next() = phe[4];
    phe[2]->prev() = phe[0];
    phe[4]->next() = phe[0];
    phe[4]->prev() = phe[2];
    phe[1]->next() = phe[5];
    phe[1]->prev() = phe[3];
    phe[3]->next() = phe[1];
    phe[3]->prev() = phe[5];
    phe[5]->next() = phe[3];
    phe[5]->prev() = phe[1];
    // b. between halfedges and edges
    e_sl->he(0) = dhe2;
    e_sl->he(1) = phe[2];
    e_sr->he(0) = dhe3;
    e_sr->he(1) = phe[3];
    phe[2]->edge() = e_sl;
    phe[3]->edge() = e_sr;
    phe[0]->edge() = e_st;
    phe[1]->edge() = e_st;
    phe[4]->edge() = e_tl;
    phe[5]->edge() = e_tr;
    dhe4->edge() = e_tl;
    dhe5->edge() = e_tr;

    // c. between halfedges and vertices
    vs->he() = dhe2;
    vt->he() = dhe5;
    vl->he() = dhe4;
    vr->he() = dhe3;
    phe[0]->target() = vs;
    phe[2]->target() = vl;
    phe[4]->target() = vt;
    phe[1]->target() = vt;
    phe[3]->target() = vs;
    phe[5]->target() = vr;
    for (unsigned int i = 0; i < cheList.size(); ++i)
        cheList[i]->target() = vt;
    // d. between halfedges and faces
    phe[0]->face() = phe[2]->face() = phe[4]->face() = f0;
    phe[1]->face() = phe[3]->face() = phe[5]->face() = f1;
    f0->he() = phe[0];
    f1->he() = phe[1];
    return 0;
}


Face *PM::CreateFace(int faceIndex) {
    tMesh->m_faces[faceIndex]->visible = true;
    return tMesh->m_faces[faceIndex];
}

Vertex *PM::CreateVertex() {
    Vertex *v = new Vertex();
    int vNum = tMesh->m_verts.size();
    if (!vNum)
        v->index() = 0;
    else
        v->index() = tMesh->m_verts[vNum - 1]->index() + 1;
    tMesh->m_verts.push_back(v);
    return v;
}

Edge *PM::CreateEdge(int edgeIndex, Halfedge *he0, Halfedge *he1) {

    Edge *e = tMesh->m_edges[edgeIndex];
    e->visible = true;
    e->setHalfEdge(he0, he1);
    return e;
}

void PM::DeleteFace(Face *f) {
    f->visible = false;
}

void PM::DeleteVertex(Vertex *v) {
    v->visible = false;
}

void PM::DeleteEdge(Edge *e) {
    e->visible = false;
}

void PM::printE(Edge *e) {
    printHE(e->he(0));
    printHE(e->he(1));
}

void PM::printHE(Halfedge *he) {
    if (he == NULL)
        std::cout << "NULL \n";
    else
        std::cout << "src: " << he->source()->index() << " trg: " << he->target()->index() << "\n";
}

bool PM::WriteVsplitRecord(const char filename[], std::vector<VSplitRecord> &vsRecList) {
    std::ofstream output(filename);
    if (!output.good()) {
        std::cerr << "Can't open file " << filename << "!" << std::endl;
        return false;
    }
    output << vsRecList.size() << "\n";
    for (unsigned int i = 0; i < vsRecList.size(); ++i) {
        VSplitRecord &vRec = vsRecList[i];
        output << vRec.vs->index() + 1 << " " << vRec.vt->index() + 1 << " " << vRec.vl->index() + 1 << " " << vRec.vr->index() + 1 << "\n";
    }
    output.close();
    return true;
}

bool PM::SaveMesh(const char filename[]) {
    std::cout << "Writing mesh " << filename << " ...";
    std::ofstream output(filename);
    if (!output.good()) {
        std::cerr << "Can't open file " << filename << "!" << std::endl;
        return false;
    }
    int vSize = tMesh->numVertices();
    int fSize = tMesh->numFaces();
    int eSize = tMesh->numEdges();
    std::cout << "\nvSize: " << vSize << "\nfSize: " << fSize << "\neSize: " << eSize << std::endl;
    for (int i = 0; i < vSize; ++i) {
        Vertex *v = tMesh->m_verts[i];
        if (!v) continue;
        output << "Vertex " << v->index() + 1 << " "
                << v->point().v[0] << " " << v->point().v[1] << " " << v->point().v[2] << "\n";
        if (!v->PropertyStr().empty()) {
            output << " {" << v->PropertyStr() << "}";
            output << "\n";
        }
    }
    for (int i = 0; i < fSize; ++i) {
        Face *f = tMesh->m_faces[i];
        Halfedge *the0 = f->he();
        Halfedge *the1 = the0->next();
        int vid0 = the0->source()->index() + 1;
        int vid1 = the0->target()->index() + 1;
        int vid2 = the1->target()->index() + 1;
        output << "Face " << f->index() + 1 << " " << vid0 << " " << vid1 << " " << vid2 << "\n";
        if (!f->PropertyStr().empty()) {
            output << " {" << f->PropertyStr() << "}";
            output << "\n";
        }
    }
    for (int i = 0; i < eSize; ++i) {
        Edge *e = tMesh->m_edges[i];
        int vid0 = e->he(0)->source()->index() + 1;
        int vid1 = e->he(0)->target()->index() + 1;
        if (!e->PropertyStr().empty()) {
            output << "Edge " << vid0 << " " << vid1 << " {" << e->PropertyStr() << "}\n";
        }
    }
    output.close();
    return true;
}

Edge *PM::GetNextCollapseEdge() {
    std::pair<Edge *, double> minE;
    minE.first = NULL;
    minE.second = 1e10;
    for (MeshEdgeIterator eit(tMesh); !eit.end(); ++eit) {
        Edge *e = *eit;
        double clen = minE.second;
        if (e->he(0)) {
            Point &p0 = e->he(0)->source()->point();
            Point &p1 = e->he(0)->target()->point();
            clen = (p0 - p1).norm();
        } else {
            Point &p0 = e->he(1)->target()->point();
            Point &p1 = e->he(1)->source()->point();
            clen = (p0 - p1).norm();
        }
        if (clen < minE.second) {
            bool canCollapse = CheckEdgeCollapseCondition(e);
            if (canCollapse) {
                minE.second = clen;
                minE.first = e;
            }
        }
    }
    if (!minE.first)
        return NULL;
    else
        return minE.first;
}

void PM::GetNextCollapseEdges(int numEdges) {
    double LIMIT = 1e10;
    std::set<Edge *> chosenEdges;
    collapseEdges.clear();

    for (int *visibleEdgesIt = visibleEdges; visibleEdgesIt != visibleEdgesEnd; ++visibleEdgesIt) {
        Edge *e = tMesh->indEdge(*visibleEdgesIt);
        if (e->visible) {
            std::set<Edge *>::iterator ePos = std::find(chosenEdges.begin(), chosenEdges.end(), e);
            if (ePos == chosenEdges.end()) {
                double clen;
                Point &p0 = e->he(0)->target()->point();
                Point &p1 = e->he(1)->target()->point();
                clen = (p0 - p1).norm();
                bool canCollapse = CheckEdgeCollapseCondition(e);
                if (clen <LIMIT && canCollapse) {
                    collapseEdges.push_back(e->index());
                    chosenEdges.insert(e);
                    Halfedge *phe[6];
                    phe[0] = e->he(0);
                    phe[1] = e->he(1);
                    phe[2] = phe[0]->next();
                    phe[4] = phe[0]->prev();
                    phe[3] = phe[1]->prev();
                    phe[5] = phe[1]->next();
                    std::vector<Vertex *> chosenVertices;

                    chosenVertices.push_back(phe[0]->target());
                    chosenVertices.push_back(phe[1]->target());
                    chosenVertices.push_back(phe[2]->target());
                    chosenVertices.push_back(phe[5]->target());
                    for (int vIndex = 0; vIndex < chosenVertices.size(); ++vIndex) {
                        for (VertexFaceIterator faceIt(chosenVertices.at(vIndex)); !faceIt.end(); ++faceIt) {
                            for (FaceEdgeIterator edgeIt(*faceIt); !edgeIt.end(); ++edgeIt) {
                                chosenEdges.insert(*edgeIt);
                            }
                        }
                    }
//                        std::cout << collapseEdges.size() << " " << e->index() << std::endl;
                    if (collapseEdges.size() == numEdges) {
                        return;
                    }
                }
            }
        }
    }
}

void *PM::FindAndCollapseEdge(void *threadarg) {
    struct thread_data *my_data = (struct thread_data *) threadarg;
    Edge *cE = my_data->edge;
    if (!cE)
        return NULL;
    my_data->pm_ptr->EdgeCollapse(cE, my_data->vsRec);
    return NULL;
}

void PM::ProcessCoarsening(int targetDisplacement) {
    pthread_t threads[NUM_THREADS];
    struct thread_data td[NUM_THREADS];

    pthread_mutex_init(&mutexLock, NULL);
    int targetVertSize = currentMeshResolution - targetDisplacement;
    for (int i = 0; i < NUM_THREADS; ++i) {
        td[i].thread_id = i;
        td[i].pm_ptr = this;
    }
    visibleEdges = new int [tMesh->numEdges()];
    int counter = 0;
    for (MeshEdgeIterator eit(tMesh); !eit.end(); ++eit) {
        visibleEdges[counter++] = (*eit)->index();
    }
    visibleEdgesEnd = thrust::set_difference(
      thrust::host,
      visibleEdges, visibleEdges + tMesh->numEdges(),
      collapseEdges.begin(), collapseEdges.end(),
      visibleEdges);
    int iteration = 0;
    for (; currentMeshResolution > targetVertSize;) {
        GetNextCollapseEdges(NUM_THREADS);
        visibleEdgesEnd = thrust::set_difference(
            thrust::host,
            visibleEdges, visibleEdgesEnd,
            collapseEdges.begin(), collapseEdges.end(),
            visibleEdges);
        int edgeCollapseSize = NUM_THREADS;
        for (int i = 0; i < edgeCollapseSize; i++) {
            // std::cout << "Thread " << i << " " << collapseEdges[i]->index() << std::endl;
            td[i].edge = tMesh->indEdge(collapseEdges[i]);
            int rc = pthread_create(&threads[i], NULL, &PM::FindAndCollapseEdge, (void *) &td[i]);
            if (rc) {
                std::cout << "Error:unable to create thread," << rc << std::endl;
                exit(-1);
            }
        }
        for (int i = 0; i < edgeCollapseSize; ++i) {
            int rc = pthread_join(threads[i], NULL);
        }
        for (int i = 0; i < edgeCollapseSize; ++i) {
            vsRecList.push_back(td[i].vsRec);
        }
        currentMeshResolution -= edgeCollapseSize;
    }
    // std::cout << "Done Coarsing" << std::endl;
}

void PM::ProcessRefinement(int targetDisplacement) {
    std::vector <XMeshLib::VSplitRecord> unrefined;
    for (int iter = 0; iter < targetDisplacement && !vsRecList.empty(); ++iter) {
        XMeshLib::VSplitRecord &vsRec = vsRecList.back();
        int err_code = VertexSplit(vsRec);
        if (err_code < 0) {
            unrefined.push_back(vsRec);
        }
        vsRecList.pop_back();
    }
    vsRecList.insert(vsRecList.end(), unrefined.begin(), unrefined.end());
    currentMeshResolution += targetDisplacement;
}


void PM::GetValidTmpMesh() {
    if (tmpMesh)
        delete tmpMesh;
    tmpMesh = new Mesh;
    tmpInd2OInd.clear();

    for (std::vector<Vertex *>::iterator viter = tMesh->m_verts.begin(); viter != tMesh->m_verts.end(); ++viter) {
        Vertex *v = *viter;
        Vertex *nv = tmpMesh->createVertex(v->index());
        if (!v) {
            nv->PropertyStr() = "invalid";
            continue;
        }
        nv->point() = v->point();
        nv->PropertyStr() = v->PropertyStr();
        nv->boundary() = v->boundary();
    }

    std::vector<Face *>::iterator fiter = tMesh->m_faces.begin();
    int cfind = 0;
    for (; fiter != tMesh->m_faces.end(); ++fiter) {
        Face *f = *fiter;
        f->index() = cfind++;
        Face *nf = tmpMesh->createFace(f->index());
        Halfedge *he[3];
        Halfedge *nhe[3];
        he[0] = f->he();
        he[1] = he[0]->next();
        he[2] = he[1]->next();
        for (int j = 0; j < 3; ++j) {
            nhe[j] = new Halfedge();
            Vertex *v1 = he[j]->target();
            Vertex *nv1 = tmpMesh->indVertex(v1->index());
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

    for (fiter = tMesh->m_faces.begin(); fiter != tMesh->m_faces.end(); ++fiter) {
        Face *f = *fiter;
        Face *nf = tmpMesh->indFace(f->index());
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
                    Edge *ne = tmpMesh->createEdge(nhe[i], NULL);
                    nhe[i]->edge() = ne;
                    ne->PropertyStr() = e->PropertyStr();
                    continue;
                }
                //get its twin: twin_he
                Halfedge *twin_he, *twin_nhe;
                Face *tf = he[i]->twin()->face();
                Face *tnf = tmpMesh->indFace(tf->index());
                twin_he = tf->he();
                twin_nhe = tnf->he();
                for (int j = 0; j < 3; ++j) {
                    if (twin_he->edge() == e) break;
                    twin_he = twin_he->next();
                    twin_nhe = twin_nhe->next();
                }
                Edge *ne = tmpMesh->createEdge(nhe[i], twin_nhe);
                nhe[i]->edge() = ne;
                twin_nhe->edge() = ne;
                ne->PropertyStr() = e->PropertyStr();
            }
        }
    }
    for (std::vector<Vertex *>::iterator vit = tmpMesh->m_verts.begin(); vit != tmpMesh->m_verts.end();) {
        Vertex *v = *vit;
        if (v->PropertyStr().compare("invalid") == 0) {
            vit = tmpMesh->m_verts.erase(vit);
        }
        else {
            tmpInd2OInd.push_back(v->index());
            v->index() = tmpInd2OInd.size() - 1;
            ++vit;
        }
    }
}
