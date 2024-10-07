#include "model.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <raymath.h>

WFModel::WFModel(const char* filename) : verts(){
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;

    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vector3 v;
            iss >> v.x >> v.y >> v.z;

            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;

            int itrash;
            int idx;

            iss >> trash;

            while (iss >> idx >> trash >> itrash >> trash >> itrash) {
                idx--;
                f.push_back(idx);
            }
            this->faces.push_back(f);
        }
    }

    this->_compute_vertex_normals();

    std::cerr << "# v# " << verts.size() << "f# " << faces.size() << std::endl;
}

int WFModel::nfaces() {
    return (int)this->faces.size();
}

int WFModel::nverts() {
    return (int)this->verts.size();
}

Vector3 WFModel::vert(int idx) {
    return this->verts[idx];
}

Vector3& WFModel::refvert(int idx) {
    return this->verts[idx];
}

std::vector<int> WFModel::face(int idx) {
    return this->faces[idx];
}

void WFModel::_compute_vertex_normals() {
    this->vert_normals.resize(nverts(), Vector3Zero());

    for (int i = 0; i < this->nfaces(); i++) {
        std::vector<int> face = this->face(i);

        Vector3 v0 = vert(face[0]);
        Vector3 v1 = vert(face[1]);
        Vector3 v2 = vert(face[2]);

        Vector3 normal = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(v2, v0), Vector3Subtract(v1, v0)));

        vert_normals[face[0]] = Vector3Add(vert_normals[face[0]], normal);
        vert_normals[face[1]] = Vector3Add(vert_normals[face[1]], normal);
        vert_normals[face[2]] = Vector3Add(vert_normals[face[2]], normal);
    }

    for (auto &n : vert_normals) {
        n = Vector3Normalize(n);
    }
}





