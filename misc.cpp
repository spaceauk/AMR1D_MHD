#include<iostream>
#include<exception>
#include<vector>
#include "defs.hpp"

real u2prim(const real Gamma,const int i,const real uu[nvar],string loc);

// Computes the primitives as a function of the Us in all ActiveBlocks
void update_prim(meshblock* dom) {
	for (int nb=0; nb<dom->lastActive; nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			for (int i=0; i<=dom->ntot; i++) {
				real uu[nvar],pp;
				for (int ii=0; ii<nvar; ii++) {
					uu[ii]=dom->u[ii][i][nb];
				}
				for (int ii=0; ii<nvar; ii++) {
					pp=u2prim(Gamma,ii,uu,"update prim");
					dom->prim[ii][i][nb]=pp;					
				}
			}
		}
	}
}

// Computes primitives from the conserved vars in a single cell
real u2prim(const real Gamma,const int i,const real uu[nvar],string loc) {
	real pp=0;
	if (i==0) {
		pp=uu[0];
		if (pp<0.) {
			cout<<loc<<") Negative density at i="<<i<<" where rho="<<pp<<endl;
			throw exception();
		}
	} else if (i==1 or i==2) {
		pp=uu[i]/uu[0];
	} else if (i==3) {
		pp=(uu[3]-0.5*(pow(uu[1],2)+pow(uu[2],2))/uu[0])*(Gamma-1.);
		if (nvar>4) pp-=0.5*(pow(uu[4],2)+pow(uu[5],2))*(Gamma-1.);
		if (pp<0.) {
			cout<<loc<<") Negative pressure at i="<<i<<" where P="<<pp<<endl;
			throw exception();
		}
	} else {
		pp=uu[i];
	}
	return pp;
}

// Computes the primitives as a function of the Us, only in one cell
real prim2u(const real Gamma,const int i,const real pp[nvar]) {
	real uu;
	if (i==0) {
		uu=pp[0];
	} else if (i==1 or i==2) {
		uu=pp[0]*pp[i];
	} else if (i==3) {
		uu=0.5*pp[0]*(pow(pp[1],2)+pow(pp[2],2))+pp[3]/(Gamma-1.);
		if (nvar>4) uu+=0.5*(pow(pp[4],2)+pow(pp[5],2));
	} else {
		uu=pp[i];
	}
	return uu;
}

int minloc(vector<real> V,int lenV) {
	int loc=0;
	real minval=abs(V[0]);
	for (int i=1; i<lenV; i++) {
		if (abs(V[i])<minval) {
			minval=abs(V[i]); loc=i;
		}
	}
	return loc;
}

void Roe_avg(real primL[nvar],real Hl,real primR[nvar],real Hr,real (&primRoe)[nvar],real &Hroe,real &aroe) { // Only for Hydro
	for (int k=1;k<nvar;k++) {
		primRoe[k]=(sqrt(primL[0])*primL[k]+sqrt(primR[0])*primR[k])/(sqrt(primL[0])+sqrt(primR[0]));
	}
	Hroe=(sqrt(primL[0])*Hl+sqrt(primR[0])*Hr)/(sqrt(primL[0])+sqrt(primR[0]));
	aroe=(Gamma-1)*(Hroe-0.5*MAG2D(primRoe[1],primRoe[2]));
}
