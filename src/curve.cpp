#include "curve.h"
#include "vertexrecorder.h"
#include <bits/stdc++.h>
using namespace std;

const float c_pi = 3.14159265358979323846f;
Curve evalBezier(const vector<Vector3f> &P, unsigned steps) {
    // Check
    if (P.size() < 4 || P.size() % 3 != 1) {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit(0);
    }

    cerr << "\t>>> evalBezier has been called with the following input:"
         << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (int i = 0; i < (int)P.size(); ++i) {
        cerr << "\t>>> " << P[i] << endl;
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;

    size_t pieces = (P.size() - 1) / 3;
    Curve R(steps * pieces + 1);
    // clang-format off
    const Matrix4f Mbezier{
         1, -3, 3, -1,
         0, 3, -6, 3, 
         0, 0, 3, -3, 
         0, 0, 0, 1,
    };
    // clang-format on

    for (size_t j = 0; j < pieces; ++j) {
        const Matrix4f Gbezier(Vector4f(P[3 * j], 0),
                               Vector4f(P[3 * j + 1], 0),
                               Vector4f(P[3 * j + 2], 0),
                               Vector4f(P[3 * j + 3], 0));

        for (unsigned i = 0, k = j * steps; i < steps; ++i, ++k) {
            float t = i / (float)steps;
            Vector4f Vt(1, t, t * t, t * t * t);
            Vector4f Tt(0, 1, 2 * t, 3 * t * t);

            R[k].V = (Gbezier * (Mbezier * Vt)).xyz();
            R[k].T = (Gbezier * (Mbezier * Tt)).xyz();
        }
    }

    const Matrix4f Gbezier(Vector4f(P[P.size() - 4], 0),
                           Vector4f(P[P.size() - 3], 0),
                           Vector4f(P[P.size() - 2], 0),
                           Vector4f(P[P.size() - 1], 0));
    Vector4f Vt(1, 1, 1, 1);
    Vector4f Tt(0, 1, 2, 3);
    R.back().V = (Gbezier * (Mbezier * Vt)).xyz();
    R.back().T = (Gbezier * (Mbezier * Tt)).xyz();

    if (!R.empty()) {
        R[0].T.normalize();
        Vector3f arbitrary_B(0, 0, 1);
        if (abs(Vector3f::dot(R[0].T, arbitrary_B)) > 0.99f) {
            arbitrary_B = Vector3f(1, 0, 0);
        }

        R[0].N = Vector3f::cross(arbitrary_B, R[0].T).normalized();
        R[0].B = Vector3f::cross(R[0].T, R[0].N).normalized();

        for (size_t i = 1; i < R.size(); ++i) {
            R[i].T.normalize();
            R[i].N = Vector3f::cross(R[i - 1].B, R[i].T).normalized();
            R[i].B = Vector3f::cross(R[i].T, R[i].N).normalized();
        }
    }

    return R;
}

Curve evalBspline(const vector<Vector3f> &P, unsigned steps) {
    // Check
    if (P.size() < 4) {
        cerr << "evalBspline must be called with 4 or more control points."
             << endl;
        exit(0);
    }

    cerr << "\t>>> evalBSpline has been called with the following input:"
         << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): " << endl;
    for (int i = 0; i < (int)P.size(); ++i) {
        cerr << "\t>>> " << P[i].x() << " " << P[i].y() << " " << P[i].z()
             << endl;
    }
    cerr << "\t>>> Steps (type steps): " << steps << endl;

    vector<Vector3f> bezierP;
    bezierP.reserve(3 * (P.size() - 3) + 1);

    // clang-format off
    Matrix4f Mbspline{
        1, -3, 3, -1, 
        4, 0, -6, 3, 
        1, 3, 3, -3, 
        0, 0, 0, 1,
    };
    Mbspline /= 6.0f;

    const Matrix4f Mbezier{
         1, -3, 3, -1,
         0, 3, -6, 3, 
         0, 0, 3, -3, 
         0, 0, 0, 1,
    };
    // clang-format on
    const Matrix4f MbsplineToBezier = Mbspline * Mbezier.inverse();

    for (size_t j = 0; j + 3 < P.size(); ++j) {
        const Matrix4f Gbspline(Vector4f(P[j], 0),
                                Vector4f(P[j + 1], 0),
                                Vector4f(P[j + 2], 0),
                                Vector4f(P[j + 3], 0));
        const Matrix4f Gbezier = Gbspline * MbsplineToBezier;

        if (j == 0) {
            bezierP.push_back(Gbezier.getCol(0).xyz());
        }
        bezierP.push_back(Gbezier.getCol(1).xyz());
        bezierP.push_back(Gbezier.getCol(2).xyz());
        bezierP.push_back(Gbezier.getCol(3).xyz());
    }

    return evalBezier(bezierP, steps);
}

Curve evalCircle(float radius, unsigned steps) {
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).

    // Preallocate a curve with steps+1 CurvePoints
    Curve R(steps + 1);

    // Fill it in counterclockwise
    for (unsigned i = 0; i <= steps; ++i) {
        // step from 0 to 2pi
        float t = 2.0f * c_pi * float(i) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f(cos(t), sin(t), 0);

        // Tangent vector is first derivative
        R[i].T = Vector3f(-sin(t), cos(t), 0);

        // Normal vector is second derivative
        R[i].N = Vector3f(-cos(t), -sin(t), 0);

        // Finally, binormal is facing up.
        R[i].B = Vector3f(0, 0, 1);
    }

    return R;
}

void recordCurve(const Curve &curve, VertexRecorder *recorder) {
    const Vector3f WHITE(1, 1, 1);
    for (int i = 0; i < (int)curve.size() - 1; ++i) {
        recorder->record_poscolor(curve[i].V, WHITE);
        recorder->record_poscolor(curve[i + 1].V, WHITE);
    }
}
void recordCurveFrames(const Curve &curve, VertexRecorder *recorder,
                       float framesize) {
    Matrix4f T;
    const Vector3f RED(1, 0, 0);
    const Vector3f GREEN(0, 1, 0);
    const Vector3f BLUE(0, 0, 1);

    const Vector4f ORGN(0, 0, 0, 1);
    const Vector4f AXISX(framesize, 0, 0, 1);
    const Vector4f AXISY(0, framesize, 0, 1);
    const Vector4f AXISZ(0, 0, framesize, 1);

    for (int i = 0; i < (int)curve.size(); ++i) {
        T.setCol(0, Vector4f(curve[i].N, 0));
        T.setCol(1, Vector4f(curve[i].B, 0));
        T.setCol(2, Vector4f(curve[i].T, 0));
        T.setCol(3, Vector4f(curve[i].V, 1));

        // Transform orthogonal frames into model space
        Vector4f MORGN = T * ORGN;
        Vector4f MAXISX = T * AXISX;
        Vector4f MAXISY = T * AXISY;
        Vector4f MAXISZ = T * AXISZ;

        // Record in model space
        recorder->record_poscolor(MORGN.xyz(), RED);
        recorder->record_poscolor(MAXISX.xyz(), RED);

        recorder->record_poscolor(MORGN.xyz(), GREEN);
        recorder->record_poscolor(MAXISY.xyz(), GREEN);

        recorder->record_poscolor(MORGN.xyz(), BLUE);
        recorder->record_poscolor(MAXISZ.xyz(), BLUE);
    }
}
