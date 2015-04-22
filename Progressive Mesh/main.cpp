#ifdef __APPLE__

#include "MeshLib_Core/Mesh.h"

#else
#include "MeshLib_Core\Mesh.h"
#include "MeshLib_Core\Iterators.h"
#endif

#include "PM.h"
using namespace XMeshLib;

int main(int argc, char **argv) {
    Mesh *cmesh = new Mesh;
    cmesh->readMFile(argv[1]);
    PM cpm(cmesh);
    cpm.ProcessCoarsening(100);
    cpm.ProcessRefinement();
    delete cmesh;
	system("pause");
    return 0;
}
