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
		// step from 0 to 2pi
		float t = 2.0f * PI * float(i) / steps;
		float sinT = sin(t), cosT = cos(t);

		// Initialize position
		// We're pivoting counterclockwise around the y-axis
		
		for (size_t j = 0; j < profile.size(); ++j, ++k) {
			auto &[V, T, N, B] = profile[j];
			surface.VV.push_back(Vector3f(V.x() * cosT, V.y(), V.x() * sinT));
			// surface.VN.push_back(Vector3f(T.x() * cosT, T.y(), T.x() * sinT));
			surface.VN.push_back(Vector3f(-N.x() * cosT, -N.y(), -N.x() * sinT));
			// surface.VN.push_back(rndVector3f());
			if (i && j) {
				surface.VF.push_back(Tup3u(k - (unsigned)profile.size() - 1, k - 1, k));
				surface.VF.push_back(Tup3u(k - (unsigned)profile.size() - 1, k, k - (unsigned)profile.size()));
				/* surface.VF.push_back(Tup3u(k - (unsigned)profile.size() - 1, k - (unsigned)profile.size(), k - 1));
				surface.VF.push_back(Tup3u(k - (unsigned)profile.size(), k - 1, k)); */
			}
		}
	}

	cerr << "\t>>> " << surface.VV.size() << "\n";

	return surface;
}

Surface makeGenCyl(const Curve &profile, const Curve &sweep) {
	Surface surface;
	surface = quad();

	if (!checkFlat(profile)) {
		cerr << "genCyl profile curve must be flat on xy plane." << endl;
		exit(0);
	}

	// TODO: Here you should build the surface.  See surf.h for details.

	cerr << "\t>>> makeGenCyl called (but not implemented).\n\t>>> Returning "
			"empty surface."
		 << endl;

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
