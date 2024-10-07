#include <vector>
#include "raylib.h"

struct WFModel {
    std::vector<Vector3> verts;
    std::vector<std::vector<int>> faces;
    std::vector<Vector3> vert_normals;

    WFModel(const char* filename);
    int nverts();
    int nfaces();
    Vector3 vert(int idx);
    Vector3& refvert(int idx);
    std::vector<int> face(int idx);

    void _compute_vertex_normals();
};