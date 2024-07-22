#include<exception>
#include "defs.hpp"

void prim2hllc_hydro(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2hll_hydro(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
void prim2RUSA_mhd(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]);
real prim2u(const real Gamma,const int i,const real* pp);
void eulerfluxes(real gamma,real (&pp)[nvar],real (&ff)[nvar]);

// Computes the fluxes 
void fluxSolver(meshblock* dom,real (&f)[nvar][nx+2][nbmax]) {
	real primL[nvar], primR[nvar], ff[nvar];

	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			for (int i=0;i<=nx;i++) {
				// Cell edge values (left vs right) - Godunov method (1st order)
				for (int ii=0;ii<nvar;ii++) {
					primL[ii]=dom->prim[ii][i][nb];
					primR[ii]=dom->prim[ii][i+1][nb];
				}
				// Calculate flux with left and right cell edge values
				if (nvar==4) {
					prim2hllc_hydro(Gamma,primL,primR,ff); 
					//prim2hll_hydro(Gamma,primL,primR,ff);
				} else {					
					prim2RUSA_mhd(Gamma,primL,primR,ff);
				}
				for (int ii=0;ii<nvar;ii++) {
					f[ii][i][nb]=ff[ii];
				}
			}
		}
	}
}

// Hydro---------------------------------------------------------------------------------
// Obtain the HLL fluxes
void prim2hll_hydro(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real sl,sr,csl,csr,sst,ek;
	real fl[nvar],fr[nvar],ust[nvar];

	csl=sqrt(gamma*primL[3]/primL[0]);
	csr=sqrt(gamma*primR[3]/primR[0]);

	sl=min(primL[1]-csl,primR[1]-csr);
	sr=max(primL[1]+csl,primR[1]+csr);

	if (sl>0.) {
		eulerfluxes(gamma,primL,ff);
	} else if (sr<0.) {
		eulerfluxes(gamma,primR,ff);
	} else if (sst>=0.) {
		eulerfluxes(gamma,primL,fl);
		eulerfluxes(gamma,primR,fr);
		for (int ii=0;ii<4;ii++) {
			real uL=prim2u(gamma,ii,primL);	
			real uR=prim2u(gamma,ii,primR);
			ff[ii]=(sr*fl[ii]-sl*fr[ii]+sl*sr*(uR-uL))/(sr-sl);
		}
	} 
}


// Obtain the HLLC fluxes
void prim2hllc_hydro(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real sl,sr,csl,csr,sst,ek;
	real fl[nvar],fr[nvar],ust[nvar];

	csl=sqrt(gamma*primL[3]/primL[0]);
	csr=sqrt(gamma*primR[3]/primR[0]);

	sl=min(primL[1]-csl,primR[1]-csr);
	sr=max(primL[1]+csl,primR[1]+csr);

	sst=(primR[3]-primL[3]+primL[0]*primL[1]*(sl-primL[1])
			      -primR[0]*primR[1]*(sr-primR[1]))
		/(primL[0]*(sl-primL[1])-primR[0]*(sr-primR[1]));

	if (sl>0.) {
		eulerfluxes(gamma,primL,ff);
	} else if (sr<0.) {
		eulerfluxes(gamma,primR,ff);
	} else if (sst>=0.) {
		ust[0]=primL[0]*(sl-primL[1])/(sl-sst);
		ust[1]=ust[0]*sst;

		ek=(0.5*primL[0]*pow(primL[1],2))+primL[3]/(gamma-1.);
		ust[3]=ust[0]*(ek/primL[0]+(sst-primL[1])*
				(sst+primL[3]/(primL[0]*(sl-primL[1]))));

		eulerfluxes(gamma,primL,fl);
		for (int ii=0;ii<4;ii++) {
			real uu=prim2u(gamma,ii,primL);	
			ff[ii]=fl[ii]+sl*(ust[ii]-uu);
		}
	} else if (sst<=0.) {
		ust[0]=primR[0]*(sr-primR[1])/(sr-sst);
		ust[1]=ust[0]*sst;

		ek=(0.5*primR[0]*pow(primR[1],2))+primR[3]/(gamma-1.);
		ust[3]=ust[0]*(ek/primR[0]+(sst-primR[1])*
				(sst+primR[3]/(primR[0]*(sr-primR[1]))));

		eulerfluxes(gamma,primR,fr);
		for (int ii=0;ii<4;ii++) {
			real uu=prim2u(gamma,ii,primR);
			ff[ii]=fr[ii]+sr*(ust[ii]-uu);
		}
	} else {
		cout<<"Error in HLLC: sl="<<sl<<", sr="<<sr<<", sst="<<sst<<endl;
		cout<<"-->   primL=["<<primL[0]<<","<<primL[1]<<","<<primL[3]<<"]; "
			"primR=["<<primR[0]<<","<<primR[1]<<","<<primR[3]<<"];"<<endl;
		throw exception();
	}
}


// Computes the euler fluxes, one cell
void eulerfluxes(real gamma,real (&pp)[nvar],real (&ff)[nvar]) {
	ff[0]=pp[0]*pp[1];	
	ff[1]=pp[0]*pow(pp[1],2)+pp[3];
	ff[2]=pp[0]*pp[1]*pp[2];
	ff[3]=pp[1]*(0.5*pp[0]*pow(pp[1],2)+gamma*pp[3]/(gamma-1.));
	if (nvar>4) {
		ff[1]+=0.5*(pow(pp[4],2)+pow(pp[5],2))-pow(pp[4],2);
		ff[3]+=pp[1]*(pow(pp[4],2)+pow(pp[5],2))-pp[4]*(pp[1]*pp[4]+0.*pp[5]);
		ff[4]=pp[1]*pp[4]-pp[4]*pp[1];
		ff[5]=pp[1]*pp[5]-0;
	}
}


/*-------------------MHD solver can be used for hydro too----------------*/
// Rusanov MHD solver
void prim2RUSA_mhd(real gamma,real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	// Left state
	real vnL=primL[1];
	real BnL=primL[4];
	real caL=sqrt((pow(primL[4],2)+pow(primL[5],2))/primL[0]);
	real aL=sqrt(gamma*primL[3]/primL[0]);
	real canL=sqrt(BnL*BnL/primL[0]);
	real cfL=sqrt(0.5*(aL*aL+caL*caL)+0.5*sqrt(pow(aL*aL+caL*caL,2)-4*(aL*aL*canL*canL)));
	real eL=(primL[3]/(primL[0]*(gamma-1)))+0.5*(pow(primL[1],2)+pow(primL[2],2))
		+0.5*(pow(primL[4],2)+pow(primL[5],2))/primL[0];
	eL=eL*primL[0];
	real HL=(eL+primL[3]+0.5*(pow(primL[4],2)+pow(primL[5],2)))/primL[0];
	real qL[6]={primL[0],primL[0]*primL[1],primL[0]*primL[2],eL,primL[4],primL[5]};

	// Right state
	real vnR=primR[1];
        real BnR=primR[4];
        real caR=sqrt((pow(primR[4],2)+pow(primR[5],2))/primR[0]);
        real aR=sqrt(gamma*primR[3]/primR[0]);
        real canR=sqrt(BnR*BnR/primR[0]);
        real cfR=sqrt(0.5*(aR*aR+caR*caR)+0.5*sqrt(pow(aR*aR+caR*caR,2)-4*(aR*aR*canR*canR)));
        real eR=(primR[3]/(primR[0]*(gamma-1)))+0.5*(pow(primR[1],2)+pow(primR[2],2))
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
