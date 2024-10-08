#include<exception>

#include "defs.hpp"


real slopelimiter(int limiter,real r,real eta) {
	real phi=0.;
	if (limiter==0) {
		phi=0.;
	} else if (limiter==1) {
		phi=1.;
	} else if (limiter==2) {
		real intm0=MIN(1.,r);
		phi=MAX(0.,intm0);
	} else if (limiter==3) {
		phi=(abs(r)+r)/(1.+abs(r));
	} else if (limiter==4) {
		real intm0=2.*r;
		real intm1=MIN(intm0,1.);
		intm0=MIN(r,2.);
		real intm2=MAX(intm0,intm1);
		phi=MAX(0.,intm2);
	} else if (limiter==5) {
		real intm0=2.*r,intm1=0.5*(1.+r);
		real intm2=MIN(intm0,intm1);
		intm0=MIN(intm2,2.); 
		phi=MAX(0.,intm0);
	} else if (limiter==6) {
		real intm0=(1.+2.*r)/3.;
		real intm1=MIN(intm0,2.);
		intm0=2.*r;
		real intm2=MIN(intm0,intm1);
		phi=MAX(0.,intm2);
	} else if (limiter==7) {
		// eta is smoothness indicator for compact 3rd order reconstruction
		real intm0=MIN(2.*r,(2.+r)/3.);
		real intm1=MIN(intm0,1.6);
		real phi_hat=MAX(-0.5*r,intm1);
		intm0=MIN((2.+r)/3.,phi_hat);
		phi_hat=MAX(0.,intm0);
		if (eta<=1.-eps) {
		  phi=(2.+r)/3.;
		} else if (eta>=1.+eps) {
		  phi=phi_hat;
		} else {
		  phi=0.5*((1.-(eta-1.)/eps)*(2.+r)/3.+(1+(eta-1.)/eps)*phi_hat);
		}
	} else {
		cout<<"Error as invalid slope limiter chosen!"<<endl;
		throw exception();
	}

	return phi;
}


