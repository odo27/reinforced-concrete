//dbeam.h
#pragma once
#include "beam.h"

class dbeam : public beam {
private:
	int d_p_;
	int A_s_p_;
public:
	dbeam(int f_ck, int f_y, int b, int h, int d, int A_s, int d_p, int A_s_p) : beam(f_ck, f_y, b, h, d, A_s) {
		d_p_ = d_p;
		A_s_p_ = A_s_p;
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
			double eps_s_p = phi_y * (c - d_p_);
			if (eps_s_p > eps_y_)
				eps_s_p = eps_y_;
			double sigma_s_p = E_s_ * eps_s_p;
			double C_s = sigma_s_p * A_s_p_;
			if (abs(C_c + C_s - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c + C_s < T_s)
				c += key;
			else
				c -= key;
		}
		double phi_y = eps_y_ / (d_ - c);
		double eps_cm = phi_y * c;
		double sig_cm = E_c_ * eps_cm;
		double C_c = sig_cm * b_ * c / 2;
		double eps_s_p = phi_y * (c - d_p_);
		if (eps_s_p > eps_y_)
			eps_s_p = eps_y_;
		double sigma_s_p = E_s_ * eps_s_p;
		double C_s = sigma_s_p * A_s_p_;
		double M_y = C_c * (d_ - c / 3) + C_s * (d_ - d_p_);
		M_y /= 1e6; //[kN-m]
		phi_y *= 1e3; //[1/m]
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
			double eps_s_p = phi_y * (c - d_p_);
			if (eps_s_p > eps_y_)
				eps_s_p = eps_y_;
			double sigma_s_p = E_s_ * eps_s_p;
			double C_s = sigma_s_p * A_s_p_;
			if (abs(C_c + C_s - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c + C_s < T_s)
				c += key;
			else
				c -= key;
		}
		double phi_y = eps_y_ / (d_ - c);
		double eps_s_p = phi_y * (c - d_p_);
		if (eps_s_p > eps_y_)
			eps_s_p = eps_y_;
		double sigma_s_p = E_s_ * eps_s_p;
		double C_s = sigma_s_p * A_s_p_;
		double M_y = C_s * (h_ / 2.0 - d_p_) + T_s * (d_ - h_ / 2.0);
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
			double C_c = sig_c * beta * c * b_;
			double eps_s_p = phi_u * (c - d_p_);
			if (eps_s_p > eps_y_)
				eps_s_p = eps_y_;
			double sigma_s_p = E_s_ * eps_s_p;
			double C_s = sigma_s_p * A_s_p_;
			if (abs(C_c + C_s - T_s) < accuracy)
				break;
			key /= 2;
			if (C_c + C_s < T_s)
				c += key;
			else
				c -= key;
		}
		double phi_u = eps_cu_ / c;
		double C_c = sig_c * beta * c * b_;
		double eps_s_p = phi_u * (c - d_p_);
		if (eps_s_p > eps_y_)
			eps_s_p = eps_y_;
		double sigma_s_p = E_s_ * eps_s_p;
		double C_s = sigma_s_p * A_s_p_;
		double M_u = C_c * (d_ - beta * c / 2) + C_s * (d_ - d_p_);
		M_u /= 1e6; //[kN-m]
		phi_u *= 1e3; //[1/m]
		point point_u(M_u, phi_u);
		return point_u;
	}
};