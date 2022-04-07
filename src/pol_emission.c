/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "functions.h"
#include "gsl/gsl_sf_hyperg.h"
#include "parameters.h"
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <stdio.h>
///////////////////////////////////////////////////////////////////////rho

double f_m(double X) {
    return 2.011 * exp(-pow(X, 1.035) / 4.7) -
           cos(X * 0.5) * exp(-pow(X, 1.2) / 2.73) - 0.011 * exp(-X / 47.2) +
           (0.011 * exp(-X / 47.2) - pow(2., -1. / 3.) / pow(3., 23. / 6.) *
                                         10000. * M_PI * pow(X, -8. / 3.)) *
               0.5 * (1. + tanh(10. * log(X / 120.)));
}

// Bessel function approximations:
double bessel_appr(int n, double x) {
    if (n == 0)
        return -log(x / 2.) - 0.5772;

    if (n == 1)
        return 1. / x;

    if (n == 2)
        return 2. / x / x;

    exit(0);
}

// B.15:
double DeltaJ_5(double X) {
    return 0.43793091 * log(1. + 0.00185777 * pow(X, 1.50316886));
}


double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B){
	
	double rho_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		
		double F_factor_K3_5=(1-exp(-pow(X_kappa,0.84)/30.)-sin(X_kappa/10.)*exp(-3*pow(X_kappa,0.471)/2.));
		double F_factor_K4=(1-exp(-pow(X_kappa,0.84)/18.)-sin(X_kappa/6.)*exp(-7*pow(X_kappa,0.5)/4.));
		double F_factor_K4_5=(1-exp(-pow(X_kappa,0.84)/12.)-sin(X_kappa/4.)*exp(-2*pow(X_kappa,0.525)));
		double F_factor_K5=(1-exp(-pow(X_kappa,0.84)/8.)-sin(3*X_kappa/8.)*exp(-9*pow(X_kappa,0.541)/4.));
			
		double Q_factor_K3_5=(17*w-3*pow(w,1./2.)+7*pow(w,1./2.)*exp(-5*w));
		double Q_factor_K4=((46./3.)*w-(5./3.)*pow(w,1./2.)+(17./3.)*pow(w,1./2.)*exp(-5*w));
		double Q_factor_K4_5=(14*w-(13./8.)*pow(w,1./2.)+(9./2.)*pow(w,1./2.)*exp(-5*w));
		double Q_factor_K5=((25./2.)*w-pow(w,1./2.)+5*pow(w,1./2.)*exp(-5*w));

		if (kappa==3.5)
			factor_F=F_factor_K3_5;;
			factor_Q=Q_factor_K3_5;
		elif (kappa==4)
			factor_F=F_factor_K4;
			factor_Q=Q_factor_K4;
		elif (kappa==4.5)
			factor_F=F_factor_K4_5;
			factor_Q=Q_factor_K4_5;
		elif (kappa==5)
			factor_F=F_factor_K5;
			factor_Q=Q_factor_K5;
		else
			print('error');
		
    		return -n_e*pow(ELECTRON_CHARGE,2)*pow(nuc,2)*pow(sth,2)/(m*SPEED_OF_LIGHT*pow(nu,3))*factor_F*factor_Q;
		}
	
	double rho_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		double wp2 = 4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    		double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    		double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    		double Thetaer = 1. / theta_e;

    		return 2. * M_PI * nu / 2. / SPEED_OF_LIGHT * wp2 * omega0 * omega0 /
           		pow(2. * M_PI * nu, 4.) * f_m(Xe) *
           		(bessel_appr(1, Thetaer) / bessel_appr(2, Thetaer) +
            		6. * theta_e) *
           		sin(theta_B) * sin(theta_B);
		}
	if(DF==kappa)
		return rho_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return rho_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}

double rho_V(double theta_e, double n_e, double nu, double B, double theta_B){
	
	double rho_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		
		double G_factor_K3_5=(1-0.17*log(1+0.447*pow(X_kappa,-1./2.)));
		double G_factor_K4=(1-0.17*log(1+0.391*pow(X_kappa,-1./2.)));
		double G_factor_K4_5=(1-0.17*log(1+0.348*pow(X_kappa,-1./2.)));
		double G_factor_K5=(1-0.17*log(1+0.313*pow(X_kappa,-1./2.)));
		
		double V_factor_K3_5=((pow(w,2)+2*w+1)/((25./8.)*pow(w,2)+4*w+1));
		double V_factor_K4=((pow(w,2)+54*w+50)/((30./11.)*pow(w,2)+134*w+50));
		double V_factor_K4_5=((pow(w,2)+43*w+38)/((7./3.)*pow(w,2)+(185./2.)*w+38));
		double V_factor_K5=((w+(13./14.))/(2*w+(13./14.)));

		if (kappa==3.5)
			factor_G=G_factor_K3_5;;
			factor_V=V_factor_K3_5;
		elif (kappa==4)
			factor_G=G_factor_K4;
			factor_V=V_factor_K4;
		elif (kappa==4.5)
			factor_G=G_factor_K4_5;
			factor_V=V_factor_K4_5;
		elif (kappa==5)
			factor_G=G_factor_K5;
			factor_V=V_factor_K5;
		else
			print('error');
		
    		return 2*n_e*pow(ELECTRON_CHARGE,2)*nuc*cos(theta)/(m*SPEED_OF_LIGHT*pow(nu,2))*bessel_appr(0, double (pow(w,-1)))/(bessel_appr(2, double (pow(w,-1))))*factor_G*factor_V;
		}
	
	double rho_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		double wp2 = 4. * M_PI * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE / ELECTRON_MASS;
    		double omega0 = ELECTRON_CHARGE * B / ELECTRON_MASS / SPEED_OF_LIGHT;
    		double Xe = theta_e * sqrt(sqrt(2.) * sin(theta_B) *
                               (1.e3 * omega0 / 2. / M_PI / nu));
    		double Thetaer = 1. / theta_e;

    		return 2.0 * M_PI * nu / SPEED_OF_LIGHT * wp2 * omega0 /
           		pow(2. * M_PI * nu, 3) * (bessel_appr(0, Thetaer) - DeltaJ_5(Xe)) /
           		bessel_appr(2, Thetaer) * cos(theta_B);
		}
	if(DF==kappa)
		return rho_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return rho_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}

///////////J for thermal and kappa
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double J_I(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double J_S_I_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double J_low_factor_I=1;
		double J_high_factor_I=1;
		double J_S_factor_I=1;
		double J_x_I=3*pow(kappa,-3./2.);
		double J_I_low_kappa= pow(X_kappa,1./3.)*sth*4*M_PI* tgamma(kappa - 4. / 3.) /
            		(pow(3, 7. / 3.) * tgamma(kappa - 2.)) * J_low_factor_I;
		double J_I_high_kappa= pow(X_kappa, -(kappa - 2) / 2.) * sth * pow(3, (kappa - 1) / 2.) *
             		(kappa - 2.) * (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
             		tgamma(kappa / 4. + 4. / 3.)*J_high_factor_I;
		
		return pow((pow(J_I_low_kappa, -J_x_I) + pow(J_I_high_kappa, -J_x_I)), -1. / J_x_I)*J_S_factor_I;
		}
	double j_I_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
		double I_I(double x) {
    			return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) + 0.9977 * pow(x, -2. / 3.)) *
           			exp(-1.8899 * pow(x, 1. / 3.));
				}
    		double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e;

    		double x = nu / nu_c;

    		return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           		SPEED_OF_LIGHT / theta_e / theta_e * I_I(x);
		}
	if(DF==kappa)
		return J_S_I_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return j_I_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}
		
double J_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double J_S_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double J_low_factor_Q=0.5;
		double J_high_factor_Q=(pow((4/5),2.)+kappa/50);
		double J_S_factor_Q=(-1);
		double J_x_Q=(37/10)*pow(kappa,-8./5.);
		double J_Q_low_kappa= pow(X_kappa,1./3.)*sth*4*M_PI* tgamma(kappa - 4. / 3.) /
            				(pow(3, 7. / 3.) * tgamma(kappa - 2.)) * J_low_factor_Q;
		double J_Q_high_kappa= pow(X_kappa, -(kappa - 2) / 2.) * sth * pow(3, (kappa - 1) / 2.) *
             				(kappa - 2.) * (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
             				tgamma(kappa / 4. + 4. / 3.)*J_high_factor_Q;
		
		return pow((pow(J_Q_low_kappa, -J_x_Q) + pow(J_Q_high_kappa, -J_x_Q)), -1. / J_x_Q)*J_S_factor_Q;
		}
	double j_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		double I_Q(double x) {
    			return 2.5651 *
           			(1 + 0.93193 * pow(x, -1. / 3.) + 0.499873 * pow(x, -2. / 3.)) *
           			exp(-1.8899 * pow(x, 1. / 3.));
				}
		double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e;

    		double x = nu / nu_c;

    		return n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / 2. / sqrt(3.) /
           		SPEED_OF_LIGHT / theta_e / theta_e * I_Q(x);
		}
	if(DF==kappa)
		return J_S_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return j_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}

double J_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double J_S_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double J_low_factor_V=(pow(0.75,2.)*pow((pow(sth,-12./25.)-1),12./25.)*pow(kappa,-66./125.)*pow(X_kappa,-7./20.)/w);
		double J_high_factor_V=(pow((7/8),2.)*pow((pow(sth,-5./2.)-1),11./25.)*pow(kappa,-11./25.)*pow(X_kappa,-1./2.)/w);
		double J_S_factor_V=(sgn(cos(theta)));
		double J_x_V=3*pow(kappa,-3./2.);
		double J_V_low_kappa= pow(X_kappa,1./3.)*sth*4*M_PI* tgamma(kappa - 4. / 3.) /
            				(pow(3, 7. / 3.) * tgamma(kappa - 2.)) * J_low_factor_V;
		double J_V_high_kappa= pow(X_kappa, -(kappa - 2) / 2.) * sth * pow(3, (kappa - 1) / 2.) *
             				(kappa - 2.) * (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
             				tgamma(kappa / 4. + 4. / 3.)*J_high_factor_V;
		
		return pow((pow(J_V_low_kappa, -J_x_V) + pow(J_V_high_kappa, -J_x_V)), -1. / J_x_V)*J_S_factor_V;
		}
	double j_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		double I_V(double x) {
    			return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
            			0.0292545 * pow(x, -0.5) + 2.03773 * pow(x, -1. / 3.)) *
           			exp(-1.8899 * pow(x, 1. / 3.));
				}
		double nu_c = 3.0 * ELECTRON_CHARGE * B * sin(theta_B) /
                      (4.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT) * theta_e *
                      theta_e;

    		double x = nu / nu_c;

    		return 2. * n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nu / tan(theta_B) /
           		3. / sqrt(3.) / SPEED_OF_LIGHT / theta_e / theta_e / theta_e *
           		I_V(x);
			}

	if(DF==kappa)
		return J_S_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return j_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}
///////////A for thermal and kappa
double a = kappa - 1. / 3.;
double b = kappa + 1.;
double c = kappa + 2. / 3.;
double z = -kappa * w;
double hyp2F1 = pow(1. - z, -a) * tgamma(c) * tgamma(b - a) /
                        (tgamma(b) * tgamma(c - a)) *
                        gsl_sf_hyperg_2F1(a, c - b, a - b + 1., 1. / (1. - z)) +
                    pow(1. - z, -b) * tgamma(c) * tgamma(a - b) /
                        (tgamma(a) * tgamma(c - b)) *
                        gsl_sf_hyperg_2F1(b, c - a, b - a + 1., 1. / (1. - z));
double A_I(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double A_S_I_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double A_low_factor_I=1;
		double A_high_factor_I=(pow((3/kappa),19./4.)+(3/5));
		double A_S_factor_I=1;
		double A_x_I=pow((-7./4. + kappa*8./5.),-43./50.);
		double A_I_low_kappa= pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) * (10. / 41.) *
            		2*M_PI / pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
            		(kappa - 1.) * kappa / (3. * kappa - 1.) * tgamma(5. / 3.) * hyp2F1*A_low_factor_I;
		double A_I_high_kappa= pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
             		((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
             		(2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) *A_high_Factor_I;

		
		return pow((pow(A_I_low_kappa, -A_x_I) + pow(A_I_high_kappa, -A_x_I)), -1. / A_x_I)*A_S_factor_I;
		}
	double A_I_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
		?
		}
	if(DF==kappa)
		return A_S_I_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return A_I_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}
		
double A_Q(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double A_S_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double A_low_factor_Q=(25/48);
		double A_high_factor_Q=(pow(21,2.)*pow(kappa,-pow(12./5.,2.))+(11/20));
		double A_S_factor_Q=(-1);
		double A_x_Q=(7./5.)*pow(kappa,-23./20.);
		double A_Q_low_kappa= pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) * (10. / 41.) *
            		2*M_PI / pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
            		(kappa - 1.) * kappa / (3. * kappa - 1.) * tgamma(5. / 3.) * hyp2F1*A_low_factor_Q;
		double A_Q_high_kappa= pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
             		((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
             		(2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) *A_high_Factor_Q;
		
		return pow((pow(A_Q_low_kappa, -A_x_Q) + pow(A_Q_high_kappa, -A_x_Q)), -1. / A_x_Q)*A_S_factor_Q;
		}
	double A_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		?
		}
	if(DF==kappa)
		return A_S_Q_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return A_Q_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}
double A_V(double theta_e, double n_e, double nu, double B, double theta_B) {
    	
	double A_S_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
		double A_low_factor_V=((77/(100*w))*pow((pow(sth,-114./50.)-1),223./500.)*pow(X_kappa,-7./20.)*pow(kappa,-7./10.));
		double A_high_factor_V=((143/10)*pow(w,-116./125.)*pow((pow(sth,-41./20.)-1),1./2.)*(pow(13,2)*
			pow(kappa,-8)+(13/2500)*kappa-(263/5000)+(47/(200*kappa)))*pow(X_kappa,-1./2.));
		double A_S_factor_V=(sgn(cos(theta)));
		double A_x_V=((61./50.)*pow(kappa,-142./125.)+(7./1000.));
		double A_V_low_kappa= pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) * (10. / 41.) *
            				2*M_PI / pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
            				(kappa - 1.) * kappa / (3. * kappa - 1.) * tgamma(5. / 3.) * hyp2F1*A_low_factor_V;
		double A_V_high_kappa= pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
             				((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
             				(2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) *A_high_Factor;

		
		return pow((pow(A_V_low_kappa, -A_x_V) + pow(A_V_high_kappa, -A_x_V)), -1. / A_x_V)*A_S_factor_V;
		}
	double A_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B) {
    		?
		}

	if(DF==kappa)
		return A_S_V_kappa(double theta_e, double n_e, double nu, double B, double theta_B);
	if(DF==thermal)
		return A_V_thermal(double theta_e, double n_e, double nu, double B, double theta_B);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////THERMAL//////////////////////////
// POLARIZED COEFFICIENTS
/////////////////////////
//kappa implementation started



// ABSORPTION
/////////////

// Use Kirchoff: a_nu = j_nu/B_nu

// ROTATION
///////////

// B.13 NOTE: make sure that, in B.13, Dexter does not mean X when he writes x
// (in the final term) Note that X = ... (B.8) UPDATE: it was apparently indeed
// a typo. while x = nu / nu_c UPDATE: it
// should be ln(X/120), not "ln x / 120". Thus I have greyed out the original,
// new one is below.
