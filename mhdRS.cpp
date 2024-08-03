#include "defs.hpp"

real prim2u(const real Gamma,const int i,const real* pp);

/*-------------------MHD solver can be used for hydro too----------------*/
// Rusanov MHD solver
void prim2RUSA_mhd(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	// Left state
	real vnL=primL[1];
	real BnL=primL[4];
	real caL=sqrt((pow(primL[4],2)+pow(primL[5],2))/primL[0]);
	real aL=sqrt(Gamma*primL[3]/primL[0]);
	real canL=sqrt(BnL*BnL/primL[0]);
	real cfL=sqrt(0.5*(aL*aL+caL*caL)+0.5*sqrt(pow(aL*aL+caL*caL,2)-4*(aL*aL*canL*canL)));
	real eL=(primL[3]/(primL[0]*(Gamma-1)))+0.5*(pow(primL[1],2)+pow(primL[2],2))
		+0.5*(pow(primL[4],2)+pow(primL[5],2))/primL[0];
	eL=eL*primL[0];
	real HL=(eL+primL[3]+0.5*(pow(primL[4],2)+pow(primL[5],2)))/primL[0];
	real qL[6]={primL[0],primL[0]*primL[1],primL[0]*primL[2],eL,primL[4],primL[5]};

	// Right state
	real vnR=primR[1];
        real BnR=primR[4];
        real caR=sqrt((pow(primR[4],2)+pow(primR[5],2))/primR[0]);
        real aR=sqrt(Gamma*primR[3]/primR[0]);
        real canR=sqrt(BnR*BnR/primR[0]);
        real cfR=sqrt(0.5*(aR*aR+caR*caR)+0.5*sqrt(pow(aR*aR+caR*caR,2)-4*(aR*aR*canR*canR)));
        real eR=(primR[3]/(primR[0]*(Gamma-1)))+0.5*(pow(primR[1],2)+pow(primR[2],2))
                +0.5*(pow(primR[4],2)+pow(primR[5],2))/primR[0];
        eR=eR*primR[0];
        real HR=(eR+primR[3]+0.5*(pow(primR[4],2)+pow(primR[5],2)))/primR[0];
        real qR[6]={primR[0],primR[0]*primR[1],primR[0]*primR[2],eR,primR[4],primR[5]};

	// Left fluxes 
	real ptL=primL[3]+0.5*(pow(primL[4],2)+pow(primL[5],2));
	real FL[6]={primL[0]*vnL, 
		primL[0]*vnL*primL[1]+ptL-BnL*primL[4],
		primL[0]*vnL*primL[2]-BnL*primL[5],
		primL[0]*vnL*HL-BnL*(primL[1]*primL[4]+primL[2]*primL[5]),
		vnL*primL[4]-BnL*primL[1],
		vnL*primL[5]-BnL*primL[2]};
	// Right fluxes
	real ptR=primR[3]+0.5*(pow(primR[4],2)+pow(primR[5],2));
        real FR[6]={primR[0]*vnR,
                primR[0]*vnR*primR[1]+ptR-BnR*primR[4],
		primR[0]*vnR*primR[2]-BnR*primR[5],
                primR[0]*vnR*HR-BnR*(primR[1]*primR[4]+primR[2]*primR[5]),
                vnR*primR[4]-BnR*primR[1],
                vnR*primR[5]-BnR*primR[2]};	
	
	// Compute the Rusanov flux
	for (int i=0;i<nvar;i++) {
		ff[i]=0.5*(FL[i]+FR[i])-0.5*max(cfL+abs(primL[1]),cfR+abs(primR[1]))*(qR[i]-qL[i]);
	} 
}

