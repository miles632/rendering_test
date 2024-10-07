#include <iostream>
#include <cassert>
#include <random>

#include "camera.hpp"
#include "model.hpp"
#include "raylib.h"
#include "raymath.h"

#define S_WIDTH 800
#define S_HEIGHT 800

const Vector3 light_dir(1,1,1);
const Vector3 eye(1,1,3);
const Vector3 center(0,0,0);
const Vector3 up(0,1,0);

float angle = 0.0f;
Vector3 axis = {0, 0, 1};
float rotationSpeed = 0.01f * DEG2RAD;

template <typename N, typename V>
N vector_cross(V t0, V t1) {
    return t0.x * t1.y - t0.y * t1.x;
}
template <typename N, typename V>
N edge_cross(V t0a, V t1b, V t2p) {
    V ab = {.x = (t1b.x - t0a.x),.y = (t1b.y - t0a.y)};
    V ap = {.x = (t2p.x - t0a.x),.y = (t2p.y - t0a.y)};

    return vector_cross<N, V>(ab, ap);
}


void line(int x0, int y0, int x1, int y1, Image &image, Color color) {
    bool steep = false;
    int x_delta = std::abs(x1 - x0);
    int y_delta = std::abs(y1 - y0);

    if (y_delta > x_delta) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        std::swap(x_delta, y_delta);
        steep = true;
    }

    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int d = 2 * y_delta - x_delta;
    int y_step = (y1 > y0) ? 1 : -1;
    int y = y0;

    for (int x = x0; x <= x1; x++) {
        if (steep) {
            ImageDrawPixel(&image, y, x, color);
        } else {
            ImageDrawPixel(&image, x, y, color);
        }

        if (d > 0) {
            y += y_step;
            d -= 2 * x_delta;
        }

        d += 2 * y_delta;
    }
}

Vector3 barycentric(Vector3 A, Vector3 B, Vector3 C, Vector3 P) {
    Vector3 v0 = Vector3Subtract(B, A);
    Vector3 v1 = Vector3Subtract(C, A);
    Vector3 v2 = Vector3Subtract(P, A);

    float d00 = Vector3DotProduct(v0, v0);
    float d01 = Vector3DotProduct(v0, v1);
    float d11 = Vector3DotProduct(v1, v1);
    float d20 = Vector3DotProduct(v2, v0);
    float d21 = Vector3DotProduct(v2, v1);

    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return (Vector3){ u, v, w};
}

Vector3 world_to_screen(Vector3 v) {
    // return Vector3{int(1-(v.x+1.)*S_WIDTH/2.+.5), int((v.y+1.)*S_HEIGHT/2.+.5), v.z};
    return Vector3{
        int((1.0 - (v.x + 1.0) / 2.0) * S_WIDTH + 0.5),  // Horizontal flip
        int((1.0 - (v.y + 1.0) / 2.0) * S_HEIGHT + 0.5), // Vertical flip
        v.z
    };
}

Color triangle_fragment(float *intensity, Vector3 bary_coords) {
    float pixel_intensity = bary_coords.x * intensity[0] + bary_coords.y * intensity[1] + bary_coords.z * intensity[2];
    pixel_intensity = std::max(0.0f, std::min(1.0f, pixel_intensity));

    Color color = (Color){
        static_cast<unsigned char>(pixel_intensity * 255),
        static_cast<unsigned char>(pixel_intensity * 255),
        static_cast<unsigned char>(pixel_intensity * 255),
        255  // Full alpha
    };
    return color;
}

void triangle(Vector3 *pts, float *intensity, float *zbuffer, Image &image, Color color) {
    Vector2 bboxmin{ std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()};
    Vector2 bboxmax{-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
    Vector2 clamp{static_cast<float>(image.width-1), static_cast<float>(image.height-1)};

    bboxmin.x = std::max(0.f, std::min(bboxmin.x, std::min(pts[0].x, std::min(pts[1].x, pts[2].x))));
    bboxmin.y = std::max(0.f, std::min(bboxmin.y, std::min(pts[0].y, std::min(pts[1].y, pts[2].y))));

    bboxmax.x = std::min(clamp.x, std::max(bboxmax.x, std::max(pts[0].x, std::max(pts[1].x, pts[2].x))));
    bboxmax.y = std::min(clamp.y, std::max(bboxmax.y, std::max(pts[0].y, std::max(pts[1].y, pts[2].y))));

    Vector3 P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vector3 bary = barycentric(pts[0], pts[1], pts[2], P);
            if (bary.x<0 || bary.y<0 || bary.z<0) continue;

            // Color color = triangle_fragment(intensity, bary)
            Color color = WHITE;

            P.z = 0;
            P.z += pts[0].z * bary.x;
            P.z += pts[1].z * bary.y;
            P.z += pts[2].z * bary.z;

            int idx = static_cast<int>(P.x + P.y * S_WIDTH);
            if (zbuffer[idx] < P.z) {
                zbuffer[idx] = P.z;
                ImageDrawPixel(&image, P.x, P.y, color);
            }
        }
    }
}



int main() {
    InitWindow(S_WIDTH, S_HEIGHT, "a");

    // imgbuf init
    Image imgbuf = GenImageColor(S_WIDTH, S_HEIGHT, BLACK);
    Texture2D text = LoadTextureFromImage(imgbuf);
    SetTargetFPS(144);

    // model init
    WFModel* model = new WFModel("obj/cube.obj");
    const Vector3 light_dir {0,1, -1};

    //zbuf init
    auto *zbuffer = new float[S_WIDTH * S_HEIGHT];
    for (int i = S_WIDTH*S_HEIGHT - 1; i >= 0; i--) {
        zbuffer[i] = -std::numeric_limits<float>::max();
    }

    // look_at(eye, center, up);
    // viewport_transformation_m(S_WIDTH/8, S_HEIGHT/8, S_WIDTH*3/4, S_HEIGHT*3/4, 1.0);

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        Image imgbuf = GenImageColor(S_WIDTH, S_HEIGHT, BLACK);
        for (int i = S_WIDTH*S_HEIGHT - 1; i >= 0; i--) {
            zbuffer[i] = -std::numeric_limits<float>::max();
        }

        angle+=rotationSpeed;
        if (angle >= 360.0f) {
            angle = fmod(angle, 360.0f);
        }

        Matrix rotationMatrix = MatrixRotate(axis, angle);

        // std::vector<int> face = model->face(0);
        // Vector3 vert1 = model->vert(face[0]);
        // std::cout << vert1.x << vert1.y << vert1.z << std::endl;

        for (int i = 0; i < model->nfaces(); i++) {
            std::vector<int> face = model->face(i);

            // Rotate vertices
            model->refvert(face[0]) = Vector3Transform(model->vert(face[0]), rotationMatrix);
            model->refvert(face[1]) = Vector3Transform(model->vert(face[1]), rotationMatrix);
            model->refvert(face[2]) = Vector3Transform(model->vert(face[2]), rotationMatrix);

            // Rotate normals too
            // model->vert_normals[face[0]] = Vector3Transform(model->vert_normals[face[0]], rotationMatrix);
            // model->vert_normals[face[1]] = Vector3Transform(model->vert_normals[face[1]], rotationMatrix);
            // model->vert_normals[face[2]] = Vector3Transform(model->vert_normals[face[2]], rotationMatrix);

            Vector3 world_coords0 = model->vert(face[0]);
            Vector3 world_coords1 = model->vert(face[1]);
            Vector3 world_coords2 = model->vert(face[2]);

            Vector3 screen_coords0 = world_to_screen(world_coords0);
            Vector3 screen_coords1 = world_to_screen(world_coords1);
            Vector3 screen_coords2 = world_to_screen(world_coords2);

            // Vector3 normal0 = model->vert_normals[face[0]];
            // Vector3 normal1 = model->vert_normals[face[1]];
            // Vector3 normal2 = model->vert_normals[face[2]];
            //
            // float intensity0 = Vector3DotProduct(Vector3Normalize(light_dir), normal0);
            // float intensity1 = Vector3DotProduct(Vector3Normalize(light_dir), normal1);
            // float intensity2 = Vector3DotProduct(Vector3Normalize(light_dir), normal2);
            // float intens[3] = { intensity0, intensity1, intensity2 };
            //
            Vector3 screen_coords[3] = { screen_coords0, screen_coords1, screen_coords2 };

            // if (intensity0 > 0 || intensity1 > 0 || intensity2 > 0) {
            //     // Use interpolation to smoothly shade the triangle based on vertex intensities
            //     triangle(screen_coords, intens, zbuffer, imgbuf, WHITE);
            // }
            float a = 1.0f;
            triangle(screen_coords, &a, zbuffer, imgbuf, WHITE);

        }


        // for (int i = 0; i < model->nfaces(); i++) {
        //     std::vector<int> face = model->face(i);
        //
        //     model->vert(face[0]) = Vector3Transform(model->vert(face[0]), rotationMatrix);
        //     model->vert(face[1]) = Vector3Transform(model->vert(face[1]), rotationMatrix);
        //     model->vert(face[2]) = Vector3Transform(model->vert(face[2]), rotationMatrix);
        //
        //     Vector3 world_coords0 = model->vert(face[0]);
        //     Vector3 world_coords1 = model->vert(face[1]);
        //     Vector3 world_coords2 = model->vert(face[2]);
        //
        //     Vector3 screen_coords0 = world_to_screen(world_coords0);
        //     Vector3 screen_coords1 = world_to_screen(world_coords1);
        //     Vector3 screen_coords2 = world_to_screen(world_coords2);
        //
        //     Vector3 normal0 = model->vert_normals[face[0]];
        //     Vector3 normal1 = model->vert_normals[face[1]];
        //     Vector3 normal2 = model->vert_normals[face[2]];
        //
        //     float intensity0 = Vector3DotProduct(Vector3Normalize(light_dir), normal0);
        //     float intensity1 = Vector3DotProduct(Vector3Normalize(light_dir), normal1);
        //     float intensity2 = Vector3DotProduct(Vector3Normalize(light_dir), normal2);
        //     float intens[3];
        //     intens[0] = intensity0; intens[1] = intensity1; intens[2] = intensity2;
        //
        //     Vector3 screen_coords[3];
        //     screen_coords[0] = screen_coords0; screen_coords[1] = screen_coords1; screen_coords[2] = screen_coords2;
        //
        //     if (intensity0 > 0 || intensity1 > 0 || intensity2 > 0) {
        //         // Use interpolation to smoothly shade the triangle based on vertex intensities
        //         triangle(screen_coords, intens, zbuffer, imgbuf, WHITE);
        //     }
        // }
        UpdateTexture(text, imgbuf.data);
        DrawTexture(text, 0, 0, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
