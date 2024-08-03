#include<exception>
#include "defs.hpp"

void prim2hllc_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2hll_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2roe_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2exact_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2RUSA_mhd(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
real slopelimiter(int limiter,real r,real eta);

// Computes the fluxes 
void fluxSolver(meshblock* dom,real (&f)[nvar][nx+2*nghosts][nbmax]) {
	real primL[nvar], primR[nvar], ff[nvar];

	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			for (int i=dom->nxminb;i<=dom->nxmax;i++) {
				// Cell edge values (left vs right)
				for (int ii=0;ii<nvar;ii++) {
					// High-res TVD slope limiter
					real dWp2, dWp1, dWm1, r, phi_l=0, phi_r=0;
					real eta=0.;
					if (dom->SLtype>0) {
						dWp2=dom->prim[ii][i+2][nb]-dom->prim[ii][i+1][nb];
						dWp1=dom->prim[ii][i+1][nb]-dom->prim[ii][i][nb];
						dWm1=dom->prim[ii][i][nb]-dom->prim[ii][i-1][nb];
						r=dWm1/(dWp1+eps);
						phi_l=slopelimiter(dom->SLtype,r,eta);
						r=dWp1/(dWp2+eps);
						phi_r=slopelimiter(dom->SLtype,1./r,eta);
					}
					// Compute cell edges (from Eq. 3.5 of Cada)					
					primL[ii]=dom->prim[ii][i][nb]+0.5*phi_l*dWp1;
					primR[ii]=dom->prim[ii][i+1][nb]-0.5*phi_r*dWp1;
				}
				// Calculate flux with left and right cell edge values
				if (nvar==4) {
					if (dom->RStype==0) {
						prim2hll_hydro(primL,primR,ff);
					} else if (dom->RStype==1) {
						prim2hllc_hydro(primL,primR,ff);
					} else if (dom->RStype==2) {	
						prim2roe_hydro(primL,primR,ff);
					} else if (dom->RStype==3) {
						prim2exact_hydro(primL,primR,ff);
					}
				} else {					
					prim2RUSA_mhd(primL,primR,ff);
				}
				for (int ii=0;ii<nvar;ii++) {
					f[ii][i][nb]=ff[ii];
				}
			}
		}
	}
}

// Computes the euler fluxes, one cell
void eulerfluxes(real (&pp)[nvar],real (&ff)[nvar]) {
	ff[0]=pp[0]*pp[1];	
	ff[1]=pp[0]*pow(pp[1],2)+pp[3];
	ff[2]=pp[0]*pp[1]*pp[2];
	ff[3]=pp[1]*(0.5*pp[0]*pow(pp[1],2)+Gamma*pp[3]/(Gamma-1.));
	if (nvar>4) {
		ff[1]+=0.5*(pow(pp[4],2)+pow(pp[5],2))-pow(pp[4],2);
		ff[3]+=pp[1]*(pow(pp[4],2)+pow(pp[5],2))-pp[4]*(pp[1]*pp[4]+0.*pp[5]);
		ff[4]=pp[1]*pp[4]-pp[4]*pp[1];
		ff[5]=pp[1]*pp[5]-0;
	}
}



