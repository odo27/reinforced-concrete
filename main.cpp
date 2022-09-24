//main.cpp
#include <iostream>
#include "beam.h"
using namespace std;

int main() {
	beam a(21, 400, 300, 550, 500, 50, 1500);
	cout << "Beam A" << endl;
	point a_pcr = a.crack();
	cout << "M_cr : " << a_pcr.moment_ << "[kN-m]";
	cout << "  Phi_cr: " << a_pcr.curvature_ << "[1/m]" << endl;
	point a_pye = a.yield_elastic(1e-4);
	cout << "M_y1 : " << a_pye.moment_ << "[kN-m]";
	cout << "  Phi_y1:  " << a_pye.curvature_ << "[1/m]" << endl;
	point a_pyp = a.yield_plastic(1e-4);
	cout << "M_y2 : " << a_pyp.moment_ << "[kN-m]";
	cout << "  Phi_y2:  " << a_pyp.curvature_ << "[1/m]" << endl;
	point a_pr = a.rupture(1e-4);
	cout << "M_u  : " << a_pr.moment_ << "[kN-m]";
	cout << "  Phi_u :   " << a_pr.curvature_ << "[1/m]" << endl << endl;
	cout << "Beam B" << endl;
	beam b(21, 400, 300, 550, 500, 50, 500);
	point b_pcr = b.crack();
	cout << "M_cr : " << b_pcr.moment_ << "[kN-m]";
	cout << "  Phi_cr: " << b_pcr.curvature_ << "[1/m]" << endl;
	point b_pye = b.yield_elastic(1e-4);
	cout << "M_y1 : " << b_pye.moment_ << "[kN-m]";
	cout << "  Phi_y1:  " << b_pye.curvature_ << "[1/m]" << endl;
	point b_pyp = b.yield_plastic(1e-4);
	cout << "M_y2 : " << b_pyp.moment_ << "[kN-m]";
	cout << "  Phi_y2:  " << b_pyp.curvature_ << "[1/m]" << endl;
	point b_pr = b.rupture(1e-4);
	cout << "M_u  : " << b_pr.moment_ << "[kN-m]";
	cout << "  Phi_u :   " <<b_pr.curvature_ << "[1/m]" << endl;
}