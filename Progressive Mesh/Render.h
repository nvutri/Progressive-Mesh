#ifdef __APPLE__

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#include "MeshLib_Core/Mesh.h"
#include "MeshLib_Core/Iterators.h"

#else
#include "MeshLib_Core\Mesh.h"
#include "MeshLib_Core\Iterators.h"
#endif

#include <iostream>

typedef struct {
    Point v1, v2, v3;
} FaceNormal;

class Render {
public:
    Render(Mesh *cMesh) {
    }

    ~Render() {
    }

public:
    static void MyInit();

    static void SetCamera(void);

    static void Render_BoxAndAxes();

    static void Render_Mesh();

    static void ComputeBoundingBox();

    static void ComputeNormal();

    static Point ComputeVertexNormal(Vertex *);

    static void BFSQueueComponent(Face *);

    static int CalculateComponents();

    static void BFSQueueBoundary(Edge *edge);

    static int CalculateBoundaries();

    static void GaussianColorCode(Vertex *);

    static void VertexNormalDisplay(Point pt);

    static void begin(int argc, char **argv, Mesh *renderMesh);

    static void display();

    static void handleKeypress(unsigned char key, int x, int y);

    static void reshape(int w, int h);

    static void mouseClick(int button, int state, int x, int y);

    static void mouseMove(int x, int y);
};