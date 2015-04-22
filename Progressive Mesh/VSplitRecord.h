#pragma once

#ifdef __APPLE__

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#include "MeshLib_Core/Mesh.h"
#include "MeshLib_Core/Iterators.h"

#else
//#include <GL\glut.h>
#include "MeshLib_Core\Mesh.h"
#include "MeshLib_Core\Iterators.h"
#endif

namespace XMeshLib {
    class VSplitRecord {
    public:
        VSplitRecord() {
            ;
        }

        ~VSplitRecord() {
            ;
        }

        Vertex *vs;
        Vertex *vl;
        Vertex *vr;
        Vertex *vt;
        Point old_vs_pt;

        Point &old_vt_pt() {
            return vt->point();
        }
    };
}
