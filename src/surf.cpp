#include "surf.h"
#include "vertexrecorder.h"
#include <bits/stdc++.h>
using namespace std;

namespace {

// We're only implenting swept surfaces where the profile curve is
// flat on the xy-plane.  This is a check function.
static bool checkFlat(const Curve &profile) {
    for (unsigned i = 0; i < profile.size(); i++)
        if (profile[i].V[2] != 0.0 || profile[i].T[2] != 0.0 ||
            profile[i].N[2] != 0.0)
            return false;

    return true;
}

// use these to solve sign problems
static void addOrientedTriangle(Surface &surface, unsigned a, unsigned b,
                                unsigned c) {
    Vector3f faceNormal = Vector3f::cross(surface.VV[b] - surface.VV[a],
                                          surface.VV[c] - surface.VV[a]);
    Vector3f vertexNormal = surface.VN[a] + surface.VN[b] + surface.VN[c];

    if (Vector3f::dot(faceNormal, vertexNormal) < 0) {
        surface.VF.push_back(Tup3u(a, c, b));
    } else {
        surface.VF.push_back(Tup3u(a, b, c));
    }
}

static void addOrientedQuad(Surface &surface, unsigned prevLower,
                            unsigned currLower, unsigned currUpper,
                            unsigned prevUpper) {
    addOrientedTriangle(surface, prevLower, currLower, currUpper);
    addOrientedTriangle(surface, prevLower, currUpper, prevUpper);
}
} // namespace

// DEBUG HELPER
Surface quad() {
    Surface ret;
    ret.VV.push_back(Vector3f(-1, -1, 0));
    ret.VV.push_back(Vector3f(+1, -1, 0));
    ret.VV.push_back(Vector3f(+1, +1, 0));
    ret.VV.push_back(Vector3f(-1, +1, 0));

    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));
    ret.VN.push_back(Vector3f(0, 0, 1));

    ret.VF.push_back(Tup3u(0, 1, 2));
    ret.VF.push_back(Tup3u(0, 2, 3));
    return ret;
}

mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());

Vector3f rndVector3f() {
    uniform_real_distribution<> rng(-1, 1);
    Vector3f e(rng(rnd), rng(rnd), rng(rnd));
    e.normalize();
    return e;
}

Surface makeSurfRev(const Curve &profile, unsigned steps) {
    if (!checkFlat(profile)) {
        cerr << "surfRev profile curve must be flat on xy plane." << endl;
        exit(0);
    }

    cerr << "\t>>> makeSurfRev called\n";
    cerr << "\t>>> Steps (type steps): " << steps << endl;

    // return quad();

    Surface surface;

    static constexpr float PI = 3.14159265358979323846f;

    for (unsigned i = 0, k = 0; i <= steps; ++i) {
        float t = 2.0f * PI * float(i) / steps;
        float sinT = sin(t), cosT = cos(t);
        Matrix4f M{
            Vector4f(cosT, 0.0f, sinT, 0.0f),
            Vector4f(0.0f, 1.0f, 0.0f, 0.0f),
            Vector4f(-sinT, 0.0f, cosT, 0.0f),
            Vector4f(0.0f, 0.0f, 0.0f, 1.0f),
        };

        for (size_t j = 0; j < profile.size(); j++, k++) {
            auto &[V, T, N, B] = profile[j];

            // P' = M * P
            Vector4f P(V.x(), V.y(), V.z(), 1.0f);
            Vector4f P_prime = M * P;

            surface.VV.push_back(P_prime.xyz());

            // N' = normalize((M^-1)^T * N)
            Vector4f N_vec(-N.x(), -N.y(), -N.z(), 0.0f);
            auto M_inv_T = M.inverse().transposed();
            auto N_prime = (M_inv_T * N_vec).xyz().normalized();
            surface.VN.push_back(N_prime);

            if (i && j) {
                addOrientedQuad(surface, k - (unsigned)profile.size() - 1,
                                k - 1, k, k - (unsigned)profile.size());
            }
        }
    }

    cerr << "\t>>> " << surface.VV.size() << "\n";

    return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep) {
    Surface surface;

    if (!checkFlat(profile)) {
        cerr << "genCyl profile curve must be flat on xy plane." << endl;
        exit(0);
    }
    bool isClosed = false;
    float totalAngle = 0.0;
    if (sweep.size() > 1 &&
        (sweep.front().V - sweep.back().V).absSquared() < 1e-4) {
        isClosed = true;
        float x = Vector3f::dot(sweep.back().N, sweep.front().N);
        float y = Vector3f::dot(sweep.back().N, sweep.front().B);
        totalAngle = atan2(y, x);
    }

    for (size_t i = 0, k = 0; i < sweep.size(); ++i) {
        const auto &s = sweep[i];

        Vector3f N_corr = s.N;
        Vector3f B_corr = s.B;

        if (isClosed) {
            float theta = -totalAngle * float(i) / float(sweep.size() - 1);
            float cosTheta = cos(theta);
            float sinTheta = sin(theta);

            N_corr = cosTheta * s.N + sinTheta * s.B;
            B_corr = -sinTheta * s.N + cosTheta * s.B;
            cerr << "Corrected (" << s.N.x() << ", " << s.N.y() << ", "
                 << s.N.z() << ") to (" << N_corr.x() << ", " << N_corr.y()
                 << ", " << N_corr.z() << ")" << endl;
        }

        for (size_t j = 0; j < profile.size(); ++j, ++k) {
            const auto &p = profile[j];

            Matrix4f M{
                Vector4f(N_corr, 0),
                Vector4f(B_corr, 0),
                Vector4f(s.T, 0),
                Vector4f(s.V, 1),
            };

            Vector4f P(p.V.x(), p.V.y(), p.V.z(), 1.0f);
            Vector4f P_prime = M * P;

            surface.VV.push_back(P_prime.xyz());

            Vector4f N(-p.N.x(), -p.N.y(), -p.N.z(), 0.0f);
            Matrix4f M_inv_T = M.inverse().transposed();
            Vector4f N_prime = M_inv_T * N;
            surface.VN.push_back(N_prime.xyz().normalized());

            if (i && j) {
                addOrientedQuad(surface, k - (unsigned)profile.size() - 1,
                                k - 1, k, k - (unsigned)profile.size());
            }
        }
    }
    return surface;
}

void recordSurface(const Surface &surface, VertexRecorder *recorder) {
    const Vector3f WIRECOLOR(0.4f, 0.4f, 0.4f);
    for (int i = 0; i < (int)surface.VF.size(); i++) {
        recorder->record(surface.VV[surface.VF[i][0]],
                         surface.VN[surface.VF[i][0]], WIRECOLOR);
        recorder->record(surface.VV[surface.VF[i][1]],
                         surface.VN[surface.VF[i][1]], WIRECOLOR);
        recorder->record(surface.VV[surface.VF[i][2]],
                         surface.VN[surface.VF[i][2]], WIRECOLOR);
    }
}

void recordNormals(const Surface &surface, VertexRecorder *recorder,
                   float len) {
    const Vector3f NORMALCOLOR(0, 1, 1);
    for (int i = 0; i < (int)surface.VV.size(); i++) {
        recorder->record_poscolor(surface.VV[i], NORMALCOLOR);
        recorder->record_poscolor(surface.VV[i] + surface.VN[i] * len,
                                  NORMALCOLOR);
    }
}

void outputObjFile(ostream &out, const Surface &surface) {

    for (int i = 0; i < (int)surface.VV.size(); i++)
        out << "v  " << surface.VV[i][0] << " " << surface.VV[i][1] << " "
            << surface.VV[i][2] << endl;

    for (int i = 0; i < (int)surface.VN.size(); i++)
        out << "vn " << surface.VN[i][0] << " " << surface.VN[i][1] << " "
            << surface.VN[i][2] << endl;

    out << "vt  0 0 0" << endl;

    for (int i = 0; i < (int)surface.VF.size(); i++) {
        out << "f  ";
        for (unsigned j = 0; j < 3; j++) {
            unsigned a = surface.VF[i][j] + 1;
            out << a << "/" << "1" << "/" << a << " ";
        }
        out << endl;
    }
}
