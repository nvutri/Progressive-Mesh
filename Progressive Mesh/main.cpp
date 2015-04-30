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
        int baseMeshResolution = atoi(argv[2]);
        int displaceMentResolution = cpm.currentMeshResolution - baseMeshResolution;
        std::cout << displaceMentResolution << std::endl;
        cpm.ProcessCoarsening(displaceMentResolution);
        Render::begin(argc, argv, &cpm);
        delete cmesh;
    }
    return 0;
}
