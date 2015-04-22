#ifdef __APPLE__

#include "MeshLib_Core/Mesh.h"

#else
#include <GL/glut.h>
#include "MeshLib_Core\Mesh.h"
#include "MeshLib_Core\Iterators.h"
#endif

#include "PM.h"
#include "Render.h"
using namespace XMeshLib;

int main(int argc, char **argv) {
    Mesh *cmesh = new Mesh;
    if (argv[1]) {
        cmesh->readMFile(argv[1]);
        PM cpm(cmesh);
        cpm.ProcessCoarsening(1000);
        cpm.ProcessRefinement();
        delete cmesh;
    }
    Mesh *baseMesh = new Mesh;
    baseMesh->readMFile("baseMesh.m");
    std::cout << baseMesh->numVertices() << " " << baseMesh->numFaces() << std::endl;
    Render::begin(argc, argv, baseMesh);
    return 0;
}
