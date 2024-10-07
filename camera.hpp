
#include <raylib.h>
#include <raymath.h>

Matrix look_at(Vector3 eye, Vector3 center, Vector3 up) {
    Vector3 z = Vector3Normalize(Vector3Subtract(eye, center));
    Vector3 x = Vector3Normalize(Vector3CrossProduct(up, z)); // (1,0,0)
    Vector3 y = Vector3Normalize(Vector3CrossProduct(z, x)); //  (0,1,0)
    Matrix res = MatrixIdentity();

    res.m0 = x.x;   res.m4 = x.y;   res.m8 = x.z;
    res.m1 = y.x;   res.m5 = y.x;   res.m9 = y.z;
    res.m2 = z.x;   res.m6 = z.y;   res.m10 = z.z;

    res.m12 = -Vector3DotProduct(x, eye);
    res.m13 = -Vector3DotProduct(y, eye);
    res.m14 = -Vector3DotProduct(z, eye);

    return res;
}

Matrix viewport_transformation_m(int x, int y, int w, int h, float depth) {
    Matrix m = MatrixIdentity();

    m.m12 = x + w / 2.0f;
    m.m13 = y + h / 2.0f;
    m.m14 = depth / 2.0f;

    m.m0 = w / 2.0f; // x
    m.m5 = h / 2.0f; // y
    m.m10 = depth / 2.0f; // z

    return m;
}


