//beam.h
#pragma once
#include <math.h>

struct point {
	double moment_;
	double curvature_;
	point(double moment, double curvature) :moment_(moment), curvature_(curvature) {}
};

class beam {
private:
	int f_ck_, f_y_; //[MPa]
	int b_, h_, d_, d_p_; //[mm]
	int A_s_; //[mm^2]
	double E_c_, E_s_; //[MPa]
	double eps_y_, eps_cu_; //[unitless]
public:
	beam(int f_ck, int f_y, int b, int h, int d, int d_p, int A_s) {
		f_ck_ = f_ck;
		f_y_ = f_y;
		b_ = b;
		h_ = h;
		d_ = d;
		d_p_ = d_p;
		A_s_ = A_s;
		if (f_ck < 40)
			E_c_ = 8500 * pow(f_ck + 4, 1.0 / 3);
		else if (f_ck < 60)
			E_c_ = 8500 * pow(f_ck + 0.1 * f_ck, 1.0 / 3);
		else
			E_c_ = 8500 * pow(f_ck + 6, 1.0 / 3);
		E_s_ = 2e5;
		eps_y_ = f_y_ / E_s_;
		eps_cu_ = 0.003;
	}
	point crack() {
		double f_r = 0.63 * sqrt(f_ck_);
		double I_g = b_ * pow(h_, 3) / 12;
		double y_t = h_ / 2.0;
		double eps_r = f_r / E_c_;
		double M_cr = f_r * I_g / y_t / 1e6; //[kN-m]
		double phi_cr = eps_r / y_t * 1e3; //[1/m]
		point point_cr(M_cr, phi_cr);
		return point_cr;
	}
	point yield_elastic(double accuracy) {
		int T_s = A_s_ * f_y_;
		double key = d_ / 2.0;
		double c = key;
		while (1) {
			double phi_y = eps_y_ / (d_ - c);
			double eps_cm = phi_y * c;
			double sig_cm = E_c_ * eps_cm;
			double C_c = sig_cm * b_ * c / 2;
			if (abs(C_c - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c < T_s)
				c += key;
			else
				c -= key;
		}
		double M_y = T_s * (d_ - c / 3) / 1e6; //[kN-m]
		double phi_y = eps_y_ / (d_ - c) * 1e3; //[1/m]
		point point_y(M_y, phi_y);
		return point_y;
	}
	point yield_plastic(double accuracy) {
		int T_s = A_s_ * f_y_;
		double key = d_ / 2.0;
		double c = key;
		int n = 1e4;
		while (1) {
			double phi_y = eps_y_ / (d_ - c);
			double C_c = 0;
			for (int k = 1; k <= n; k++) {
				double y_k = (k - 0.5) * c / n;
				double eps_c = phi_y * y_k;
				double sig_c = f_ck_ * (2 * eps_c / 0.002 - pow(eps_c / 0.002, 2));
				C_c += sig_c * b_ * c / n;
			}
			if (abs(C_c - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c < T_s)
				c += key;
			else
				c -= key;
		}
		double M_y = T_s * (d_ - h_ / 2.0);
		double phi_y = eps_y_ / (d_ - c);
		for (int k = 1; k <= n; k++) {
			double y_k = (k - 0.5) * c / n;
			double eps_c = phi_y * y_k;
			double sig_c = f_ck_ * (2 * eps_c / 0.002 - pow(eps_c / 0.002, 2));
			M_y += sig_c * b_ * c / n * (y_k + h_ / 2.0 - c);
		}
		M_y /= 1e6; //[kN-m]
		phi_y *= 1e3; //[1/m]
		point point_y(M_y, phi_y);
		return point_y;
	}
	point rupture(double accuracy) {
		int T_s = A_s_ * f_y_;
		double beta;
		double sig_c = 0.85 * f_ck_;
		if (f_ck_ < 28)
			beta = 0.85;
		else if (f_ck_ < 56)
			beta = 0.85 - (f_ck_ - 28.0) / 140;
		else
			beta = 0.65;
		double key = d_ / 2.0;
		double c = key;
		while (1) {
			double phi_u = eps_cu_ / c;
			double eps_s = phi_u * (d_ - c);
			double C_c = sig_c * beta * c * b_;
			if (abs(C_c - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c < T_s)
				c += key;
			else
				c -= key;
		}
		double M_u = T_s * (d_ - beta * c / 2) / 1e6; //[kN-m]
		double phi_u = eps_cu_ / c * 1e3; //[1/m]
		point point_u(M_u, phi_u);
		return point_u;
	}
};