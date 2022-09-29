//main.cpp
#include <iostream>
#include "beam.h"
#include "dbeam.h"
using namespace std;

int main() {
	cout << "Beam A" << endl;
	beam a(21, 400, 300, 550, 500, 1500);
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
	beam b(21, 400, 300, 550, 500, 500);
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
	cout << "  Phi_u :   " << b_pr.curvature_ << "[1/m]" << endl << endl;
	cout << "Beam C" << endl;
	dbeam c(24, 400, 400, 600, 540, 3600, 60, 1000);
	point c_cr = c.crack();
	cout << "M_cr : " << c_cr.moment_ << "[kN-m]";
	cout << "  Phi_cr: " << c_cr.curvature_ << "[1/m]" << endl;
	point c_pye = c.yield_elastic(1e-4);
	cout << "M_y1 : " << c_pye.moment_ << "[kN-m]";
	cout << "  Phi_y1:  " << c_pye.curvature_ << "[1/m]" << endl;
	point c_pyp = c.yield_plastic(1e-4);
	cout << "M_y2 : " << c_pyp.moment_ << "[kN-m]";
	cout << "  Phi_y2:  " << c_pyp.curvature_ << "[1/m]" << endl;
	point c_u = c.rupture(1e-4);
	cout << "M_u  :  " << c_u.moment_ << "[kN-m]";
	cout << "  Phi_u :   " << c_u.curvature_ << "[1/m]" << endl;
}