#include "defs.hpp"

// AUSM+-up scheme (JCP 214, 2006)
void AUSMplus_up(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real aL=sqrt(Gamma*primL[3]/primL[0]);
	real ML=primL[1]/aL;
	real aR=sqrt(Gamma*primR[3]/primR[0]);
	real MR=primR[1]/aR;

	real MLplus,MRminus;
	real pLplus,pRminus;
	string p_polytype="2nd";
	if (abs(ML)<=1) {
		// Van Leer splitting for Mach number
		MLplus=0.25*pow(ML+1,2);
		if (p_polytype=="1st") {
			// 1st order polynomial expansions for pressure splitting
			pLplus=0.5*primL[3]*(1+ML);
		} else {
			// 2nd order polynomial expressions for pressure splitting
			pLplus=0.25*primL[3]*pow(ML+1,2)*(2-ML);
		}
	} else {
		MLplus=0.5*(ML+abs(ML));
		pLplus=0.5*primL[3]*(ML+abs(ML))/ML;
	}
	if (abs(MR)<=1) {
                MRminus=-0.25*pow(MR-1,2);
		if (p_polytype=="1st") {
			pRminus=0.5*primR[3]*(1-MR);
		} else {
                	pRminus=0.25*primR[3]*pow(MR-1,2)*(2+MR);
		}
        } else {
                MRminus=0.5*(MR-abs(MR));
               	pRminus=0.5*primR[3]*(MR-abs(MR))/MR;
        }	 
	real Mp=MLplus+MRminus;
	real pp=pLplus+pRminus;

	// Mach number weighted average
	ff[0]=0.5*Mp*(primL[0]*aL+primR[0]*aR);
        ff[1]=0.5*Mp*(primL[0]*primL[1]*aL+primR[0]*primR[1]*aR);
        ff[2]=0.5*Mp*(primL[0]*primL[2]*aL+primR[0]*primR[2]*aR);
        ff[3]=0.5*Mp*((0.5*primL[0]*pow(primL[1],2)+Gamma*primL[3]/(Gamma-1.))*aL +
	      (0.5*primR[0]*pow(primR[1],2)+Gamma*primR[3]/(Gamma-1.))*aR);

	// Numerical dissipation that facilitates upwinding
	ff[0]-=0.5*abs(Mp)*(-primL[0]*aL+primR[0]*aR);
        ff[1]-=0.5*abs(Mp)*(-primL[0]*primL[1]*aL+primR[0]*primR[1]*aR);
        ff[2]-=0.5*abs(Mp)*(-primL[0]*primL[2]*aL+primR[0]*primR[2]*aR);
        ff[3]-=0.5*abs(Mp)*(-(0.5*primL[0]*pow(primL[1],2)+Gamma*primL[3]/(Gamma-1.))*aL +
	      (0.5*primR[0]*pow(primR[1],2)+Gamma*primR[3]/(Gamma-1.))*aR);

	ff[1]+=pp;
}
