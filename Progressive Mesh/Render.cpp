/**
* Student: Tri Nguyen
* LSU ID: 89-6534216
* PAWS: tngu418
*/
#include "Render.h"
#include "PM.h"
#include <fstream>
#include <sstream>
#include <queue>
#include <set>

float whratio;
int win_height, win_width;
GLUquadricObj *obj;

/* Some variables to measure mouse movement */
int mousePositionX0 = 0, mousePositionY0 = 0;
int mouseButton = 0;

/* Some variables to describe the object dimension and camera placement */
float cameraPosition[3] = {0, 0, 2};    // Default camera position
float objCenter[3] = {0, 0, 0};         // Default object center
float boxMin[3] = {0, 0, 0};
float boxMax[3] = {0, 0, 0};
float axislen = 1.414;

/* Some variables to control interactive transformations */
float my_Rotation_x = 0, my_Rotation_y = 0;
float my_Translation[3] = {0, 0, 0};

std::vector<FaceNormal> normalVertexes;
std::vector<double> vertexesGaussian;

bool pressedVertexNormalize = false;
bool pressedGaussianCurvature = false;
Mesh *pmesh;
XMeshLib::PM *pCPM;

void Render::display() {
    SetCamera();
    Render_BoxAndAxes();
    Render_Mesh();
    glutSwapBuffers();
}

// Called when a key is pressed
void Render::handleKeypress(unsigned char key, int x, int y) {
    if (key == 102) {
        // f key.
        pressedVertexNormalize = false;
        display();
    } else if (key == 118) {
        // v key.
        pressedVertexNormalize = true;
        display();
    } else if (key == 107) {
        // k key.
        pressedGaussianCurvature = !pressedGaussianCurvature;
        display();
    } else if (key == 117) {
        pCPM->ProcessRefinement(100);
        ComputeNormal();
        display();
    } else if (key == 100) {
        pCPM->ProcessCoarsening(100);
        ComputeNormal();
        display();
    }
}

void Render::reshape(int w, int h) {
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    win_height = h;
    win_width = w;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    whratio = float((double) w / (double) h);    //A Commonly Suggested Setting: set ratio in gluPerspective to the aspect ratio of the associated viewport
    gluPerspective(60, whratio, axislen * 0.01, axislen * 5);
    glMatrixMode(GL_MODELVIEW);    //change back to modelview
    glutPostRedisplay();
}

void Render::begin(int argc, char **argv, XMeshLib::PM *cpm) {
    pCPM = cpm;
    pmesh = cpm->tMesh;
    ComputeBoundingBox();
    ComputeNormal();

//    CalculateComponents();
//    CalculateBoundaries();

    // OpenGL Routines.
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(700, 700);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Homework 1 submitted by XXX");
    MyInit();

    // Register Callback Functions.
    glutDisplayFunc(display);
    glutKeyboardFunc(handleKeypress);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);

    glutMainLoop();
}


void Render::SetCamera(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(cameraPosition[0], cameraPosition[1], cameraPosition[2], objCenter[0], objCenter[1], objCenter[2], 0, 1, 0);

    glTranslatef(my_Translation[0], my_Translation[1], my_Translation[2]);

    glTranslatef(objCenter[0], objCenter[1], objCenter[2]);    //before doing rotation to the object, move the object center to the origin

    glRotatef(my_Rotation_y, 0.0, 1.0, 0.0);
    glRotatef(my_Rotation_x, 1.0, 0.0, 0.0);

    glTranslatef(-objCenter[0], -objCenter[1], -objCenter[2]);
}


void Render::mouseMove(int x, int y) {
    double movingScale = axislen / win_height;  // just a scaling factor to make the mouse moving not too sensitive
    /* Rotation */
    if (mouseButton == GLUT_LEFT_BUTTON) {
        ////////////do something////////////////
        my_Rotation_y += x - mousePositionX0;
        my_Rotation_x += y - mousePositionY0;
    }

    /* xy translation */
    if (mouseButton == GLUT_MIDDLE_BUTTON) {
        ////////////do something ////////////////
        my_Translation[0] += movingScale * (x - mousePositionX0);
        my_Translation[1] -= movingScale * (y - mousePositionY0);
    }

    /* zoom in and out */
    if (mouseButton == GLUT_RIGHT_BUTTON) { // suppose we want to make moving up as zooming out
        my_Translation[2] += movingScale * (y - mousePositionY0);
    }
    mousePositionX0 = x;
    mousePositionY0 = y;
    glutPostRedisplay();
}

void Render::mouseClick(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_LEFT_BUTTON;
    else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_MIDDLE_BUTTON;
    else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
        mouseButton = GLUT_RIGHT_BUTTON;

    mousePositionX0 = x;
    mousePositionY0 = y;
    return;
}

/*
Rendering box and axes
 */
void Render::Render_BoxAndAxes() {

    float axiswidth = axislen / 100;

    glMatrixMode(GL_MODELVIEW);

    glColor3f(1, 1, 1);

    glPushMatrix();
    //bounding box
    glTranslatef(objCenter[0], objCenter[1], objCenter[2]);
    glutWireCube(axislen);
    glTranslatef(-axislen / 2, -axislen / 2, -axislen / 2);
    glutSolidSphere(axiswidth * 1.5, 10, 10);

    //x-axis
    glColor3f(1, 0, 0);
    glPushMatrix();
    glRotatef(90, 0, 1, 0);
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(axislen, 0, 0);
    glRotatef(90, 0, 1, 0);
    glutWireCone(axiswidth * 1.5, axiswidth * 3, 10, 10);
    glPopMatrix();


    //y-axis
    glColor3f(0, 1, 0);
    glPushMatrix();
    glTranslatef(0, axislen, 0);
    glRotatef(-90, 1, 0, 0);
    glutWireCone(axiswidth * 1.5, axiswidth * 3, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glRotatef(-90, 1, 0, 0);
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();

    //z-axis
    glColor3f(0, 0, 1);
    glPushMatrix();
    glTranslatef(0, 0, axislen);
    glRotatef(-90, 0, 0, 1);
    glutWireCone(axiswidth * 1.5, axiswidth * 3, 10, 10);
    glPopMatrix();

    glPushMatrix();
    gluCylinder(obj, axiswidth, axiswidth, axislen, 5, 5);
    glPopMatrix();

    glPopMatrix();
}

void Render::GaussianColorCode(Vertex *v) {
    if (pressedGaussianCurvature) {
        double gc = vertexesGaussian[v->index()];
        if (gc < -0.05) {
            glColor3f(0.0f, 0.5f, 0.0f);
        } else if (gc > 0.05) {
            glColor3f(0.5f, 0.0f, 0.0f);
        } else {
            glColor3f(0.7f, 0.7f, 0.7f);
        }
    }
}

void Render::VertexNormalDisplay(Point pt) {
    if (pressedVertexNormalize) {
        glNormal3dv(pt.v);
    }
}

/*
Render triangle meshes.
 */
void Render::Render_Mesh() {
    glColor3f(0.7f, 0.7f, 0.7f);
    //traverse all the face and draw them
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit) {
        Face *f = *fit;
        if (f) {
            Halfedge *he = f->he();
            Vertex *v1 = he->source();
            Vertex *v2 = he->target();
            Vertex *v3 = he->next()->target();
            glBegin(GL_TRIANGLES);
            glNormal3d(0, 0, 0);
            FaceNormal normFace = normalVertexes[f->index()];
            VertexNormalDisplay(normFace.v1);
            GaussianColorCode(v1);
            glVertex3dv(v1->point().v);
            VertexNormalDisplay(normFace.v2);
            GaussianColorCode(v2);
            glVertex3dv(v2->point().v);
            VertexNormalDisplay(normFace.v3);
            GaussianColorCode(v3);
            glVertex3dv(v3->point().v);
            glEnd();
        }
    }
}


/*
Computer object bounding box for correct view.
*/
void Render::ComputeBoundingBox() {
    objCenter[0] = objCenter[1] = objCenter[2] = 0;
    boxMin[0] = boxMin[1] = boxMin[2] = 1e5;
    boxMax[0] = boxMax[1] = boxMax[2] = -1e5;
    int vNum = int(pmesh->numVertices());

    for (MeshVertexIterator vit(pmesh); !vit.end(); ++vit) {
        Vertex *v = *vit;
        if (v) {
            Point &p = v->point();
            for (int j = 0; j < 3; ++j) {
                float value = float(p.v[j]);
                objCenter[j] += value;
                if (boxMax[j] < value)
                    boxMax[j] = value;
                if (boxMin[j] > value)
                    boxMin[j] = value;
            }
        }

    }
    axislen = float(sqrt((boxMax[2] - boxMin[2]) * (boxMax[2] - boxMin[2]) +
            (boxMax[1] - boxMin[1]) * (boxMax[1] - boxMin[1]) +
            (boxMax[0] - boxMin[0]) * (boxMax[0] - boxMin[0])));

    objCenter[0] /= vNum;
    objCenter[1] /= vNum;
    objCenter[2] /= vNum;

    cameraPosition[0] = objCenter[0];
    cameraPosition[1] = objCenter[1];
    cameraPosition[2] = objCenter[2] + float(axislen * 1.5);

    std::cout << objCenter[0]
            << " " << objCenter[1]
            << " " << objCenter[2]
            << " " << cameraPosition[0]
            << " " << cameraPosition[1]
            << " " << cameraPosition[2] << "\n";
}


void Render::MyInit() {
    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0, 0, 0, 0);
    glShadeModel(GL_SMOOTH);
    glPolygonMode(GL_FRONT, GL_FILL);

    obj = gluNewQuadric();    //only for drawing spheres and cones

    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHT0);

    // Create light components
    GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
    GLfloat diffuseLight[] = {0.8f, 0.8f, 0.8, 1.0f};
    GLfloat specularLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat position[] = {
            cameraPosition[0],
            cameraPosition[1],
            cameraPosition[2],
            1.0f
    }; // the light is on the camera position

    // Assign created components to GL_LIGHT0
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
}


/*
Computer normal for all faces.
*/
void Render::ComputeNormal() {
    normalVertexes.resize((unsigned long) pmesh->numFaces() * 2);
    vertexesGaussian.resize((unsigned long) pmesh->numVertices() * 2);
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit) {
        Face *f = *fit;
        if (f) {
            Halfedge *he = f->he();
            // Get the vertex normal for every vertex.
            FaceNormal fn;
            fn.v1 = ComputeVertexNormal(he->source());
            fn.v2 = ComputeVertexNormal(he->target());
            fn.v3 = ComputeVertexNormal(he->next()->target());

            // Save the normal vector to a similar structure
            // for using glNormal later.
            if (f->index() >= normalVertexes.size()) {
                normalVertexes.resize(f->index() * 2);
            }
            normalVertexes[f->index()] = fn;
        }
    }
}

/**
* Calculate vertex normal for a vertex.
*/
Point Render::ComputeVertexNormal(Vertex *v) {
    Point total_weighted_n_face;
    float total_a = 0;
    for (VertexFaceIterator vfit(v); !vfit.end(); ++vfit) {
        Face *f = *vfit;
        Halfedge *he = f->he();
        int counter = 0;
        while (he->source()->index() != v->index() && counter < 3) {
            he = he->next();
            ++counter;
        }
        Point &p1 = he->source()->point();
        Point &p2 = he->target()->point();
        Point &p3 = he->next()->target()->point();
        Point p1p2 = p2 - p1;
        Point p1p3 = p3 - p1;
        Point n_face = (p1p2 ^ p1p3) / ((p1p2 ^ p1p3).norm());
        double a_p1p2 = acos((p1p2 * p1p3) / (p1p2.norm() * p1p3.norm()));
        n_face = n_face * a_p1p2;
        total_weighted_n_face = total_weighted_n_face + n_face;
        total_a += a_p1p2;
    }
    vertexesGaussian[v->index()] = 2 * M_PI - total_a;
    return total_weighted_n_face / total_a;
}


/**
* Calculate number of components using BFS.
*/
const std::string TRUE_MARK("1");
const std::string FALSE_MARK;

int Render::CalculateComponents() {
    int counter = 0;
    for (MeshFaceIterator fit(pmesh); !fit.end(); ++fit) {
        Face *f = *fit;
        if (f) {
            std::string propertyStr = f->PropertyStr();
            if (propertyStr.compare(FALSE_MARK) == 0) {
                ++counter;
                f->PropertyStr() = TRUE_MARK;
                BFSQueueComponent(f);
            }
        }

    }
    std::cout << "Components: " << counter << std::endl;
    return counter;
}

void Render::BFSQueueComponent(Face *f) {
    std::queue<Face *> faceQueue;
    faceQueue.push(f);
    int numFaces = 0;
    std::set<Edge *> edgeSet;
    std::set<Vertex *> vertexSet;
    while (!faceQueue.empty()) {
        Face *fEntry = faceQueue.front();
        ++numFaces;
        vertexSet.insert(fEntry->he()->source());
        vertexSet.insert(fEntry->he()->target());
        vertexSet.insert(fEntry->he()->next()->target());
        edgeSet.insert(fEntry->he()->edge());
        edgeSet.insert(fEntry->he()->next()->edge());
        edgeSet.insert(fEntry->he()->prev()->edge());
        faceQueue.pop();
        Halfedge *he = fEntry->he();
        std::vector<Halfedge *> neighborHalfEdges;
        neighborHalfEdges.push_back(he->ccw_rotate_about_source());
        neighborHalfEdges.push_back(he->clw_rotate_about_source());
        neighborHalfEdges.push_back(he->ccw_rotate_about_target());
        neighborHalfEdges.push_back(he->clw_rotate_about_target());
        for (int i = 0; i < neighborHalfEdges.size(); ++i)
            if (neighborHalfEdges[i] && neighborHalfEdges[i]->face()) {
                Face *neighborFace = neighborHalfEdges[i]->face();
                if (neighborFace->PropertyStr().compare(FALSE_MARK) == 0) {
                    neighborFace->PropertyStr() = TRUE_MARK;
                    faceQueue.push(neighborFace);
                }
            }
    }
    int genus = 1 - (numFaces - (int) edgeSet.size() + (int) vertexSet.size()) / 2;
    std::cout << "Genus: " << genus << std::endl;
}


/*
Calculate number of boundaries using BFS.
 */

int Render::CalculateBoundaries() {
    int counter = 0;
    for (MeshEdgeIterator fit(pmesh); !fit.end(); ++fit) {
        Edge *edge = *fit;
        std::string propertyStr = edge->PropertyStr();
        if (edge->boundary() && propertyStr.compare(FALSE_MARK) == 0) {
            ++counter;
            edge->PropertyStr() = TRUE_MARK;
            BFSQueueBoundary(edge);
        }
    }
    std::cout << "Boundaries: " << counter << std::endl;
    return counter;
}

void Render::BFSQueueBoundary(Edge *edge) {
    std::queue<Edge *> edgeQueue;
    edgeQueue.push(edge);
    while (!edgeQueue.empty()) {
        Edge *edgeEntry = edgeQueue.front();
        edgeQueue.pop();
        Halfedge *he = edgeEntry->he(0);
        if (!he)
            he = edgeEntry->he(1);
        Vertex *source = he->source();
        std::vector<Halfedge *> neighborHalfEdges;
        neighborHalfEdges.push_back(source->most_ccw_in_halfedge());
        neighborHalfEdges.push_back(source->most_ccw_out_halfedge());
        neighborHalfEdges.push_back(source->most_clw_in_halfedge());
        neighborHalfEdges.push_back(source->most_clw_out_halfedge());
        for (int i = 0; i < neighborHalfEdges.size(); ++i)
            if (neighborHalfEdges[i] && neighborHalfEdges[i]->face()) {
                Edge *neighborEdge = neighborHalfEdges[i]->edge();
                if (neighborEdge->PropertyStr().compare(FALSE_MARK) == 0) {
                    neighborEdge->PropertyStr() = TRUE_MARK;
                    edgeQueue.push(neighborEdge);
                }
            }
    }
}
