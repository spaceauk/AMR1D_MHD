#include<vector>
#include "defs.hpp"

real prim2u(const real Gamma,const int i,const real* pp);
void eulerfluxes(real (&pp)[nvar],real (&ff)[nvar]);
void Newton(real (&primL)[nvar],real (&primR)[nvar],real aL,real aR,real &pstar,real &ustar);
void PREFUN(real &F,real &FD,real p,real dk,real pk,real ak);
int minloc(vector<real> V,int lenV);
void Roe_avg(real primL[nvar],real Hl,real primR[nvar],real Hr,real (&primRoe)[nvar],real &Hroe, real &aroe);

// Hydro---------------------------------------------------------------------------------
// Obtain the HLL fluxes
void prim2hll_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real sl,sr,csl,csr,sst,ek;
	real fl[nvar],fr[nvar],ust[nvar];

	csl=sqrt(Gamma*primL[3]/primL[0]);
	csr=sqrt(Gamma*primR[3]/primR[0]);

	sl=min(primL[1]-csl,primR[1]-csr);
	sr=max(primL[1]+csl,primR[1]+csr);

	if (sl>0.) {
		eulerfluxes(primL,ff);
	} else if (sr<0.) {
		eulerfluxes(primR,ff);
	} else if (sst>=0.) {
		eulerfluxes(primL,fl);
		eulerfluxes(primR,fr);
		for (int ii=0;ii<4;ii++) {
			real uL=prim2u(Gamma,ii,primL);	
			real uR=prim2u(Gamma,ii,primR);
			ff[ii]=(sr*fl[ii]-sl*fr[ii]+sl*sr*(uR-uL))/(sr-sl);
		}
	} 
}


// Obtain the HLLC fluxes
void prim2hllc_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real sl,sr,csl,csr,sst,ek;
	real fl[nvar],fr[nvar],ust[nvar];

	csl=sqrt(Gamma*primL[3]/primL[0]);
	csr=sqrt(Gamma*primR[3]/primR[0]);

	sl=min(primL[1]-csl,primR[1]-csr);
	sr=max(primL[1]+csl,primR[1]+csr);

	sst=(primR[3]-primL[3]+primL[0]*primL[1]*(sl-primL[1])
			      -primR[0]*primR[1]*(sr-primR[1]))
		/(primL[0]*(sl-primL[1])-primR[0]*(sr-primR[1]));

	if (sl>0.) {
		eulerfluxes(primL,ff);
	} else if (sr<0.) {
		eulerfluxes(primR,ff);
	} else if (sst>=0.) {
		ust[0]=primL[0]*(sl-primL[1])/(sl-sst);
		ust[1]=ust[0]*sst;

		ek=(0.5*primL[0]*pow(primL[1],2))+primL[3]/(Gamma-1.);
		ust[3]=ust[0]*(ek/primL[0]+(sst-primL[1])*
				(sst+primL[3]/(primL[0]*(sl-primL[1]))));

		eulerfluxes(primL,fl);
		for (int ii=0;ii<4;ii++) {
			real uu=prim2u(Gamma,ii,primL);	
			ff[ii]=fl[ii]+sl*(ust[ii]-uu);
		}
	} else if (sst<=0.) {
		ust[0]=primR[0]*(sr-primR[1])/(sr-sst);
		ust[1]=ust[0]*sst;

		ek=(0.5*primR[0]*pow(primR[1],2))+primR[3]/(Gamma-1.);
		ust[3]=ust[0]*(ek/primR[0]+(sst-primR[1])*
				(sst+primR[3]/(primR[0]*(sr-primR[1]))));

		eulerfluxes(primR,fr);
		for (int ii=0;ii<4;ii++) {
			real uu=prim2u(Gamma,ii,primR);
			ff[ii]=fr[ii]+sr*(ust[ii]-uu);
		}
	} else {
		cout<<"Error in HLLC: sl="<<sl<<", sr="<<sr<<", sst="<<sst<<endl;
		cout<<"-->   primL=["<<primL[0]<<","<<primL[1]<<","<<primL[3]<<"]; "
			"primR=["<<primR[0]<<","<<primR[1]<<","<<primR[3]<<"];"<<endl;
		throw exception();
	}
}


// Roe solver which involves linearizing the nonlinear Riemann problem before solving it exactly
void prim2roe_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {	
	real fl[nvar],fr[nvar];

	// Left state
	real dl=primL[0], ul=primL[1], vl=primL[2], pl=primL[3];
	real al=sqrt(Gamma*pl/dl);
	real el=pl/(dl*(Gamma-1));
        real El=0.5*dl*(ul*ul+vl*vl)+dl*el;
	real Hl=(El+pl)/dl;

	// Right state
	real dr=primR[0], ur=primR[1], vr=primR[2], pr=primR[3];
	real ar=sqrt(Gamma*pr/dr);
	real er=pr/(dr*(Gamma-1));
        real Er=0.5*dr*(ur*ur+vr*vr)+dr*er;
	real Hr=(Er+pr)/dr;

	// Compute Roe-averaged quantitites
	real primRoe[nvar]; real Hroe, aroe;
	Roe_avg(primL,Hl,primR,Hr,primRoe,Hroe,aroe);

	// Compute Roe-averaged eigenvalues 
	real lambda[nvar];
	lambda[0]=primRoe[1]-aroe; 
	lambda[1]=primRoe[1]; lambda[2]=primRoe[2];
	lambda[3]=primRoe[1]+aroe;

	// Compute right eigenvectors based on Roe avg quantities
	real Kroe[nvar][nvar];
	Kroe[0][0]=1.;
	Kroe[1][0]=primRoe[1]-aroe;
	Kroe[2][0]=primRoe[2];
	Kroe[3][0]=Hroe-primRoe[1]*aroe;
	Kroe[0][1]=1.;
	Kroe[1][1]=primRoe[1];
	Kroe[2][1]=primRoe[2];
	Kroe[3][1]=0.5*MAG2D(primRoe[1],primRoe[2]);
	Kroe[0][2]=0.;
	Kroe[1][2]=0.;
	Kroe[2][2]=1.;
	Kroe[3][2]=primRoe[2];
	Kroe[0][3]=1.;
	Kroe[1][3]=primRoe[1]+aroe;
	Kroe[2][3]=primRoe[2];
	Kroe[3][3]=Hroe+primRoe[1]*aroe;

	// Compute wave strengths \alpha based on (11.68)-(11.70)
	real qL[nvar], qR[nvar];
	for (int k=0;k<nvar;k++) {
		qL[k]=prim2u(Gamma,k,primL);
		qR[k]=prim2u(Gamma,k,primR);
	}
	real alpha[nvar];
	alpha[2]=(qR[2]-qL[2])-primRoe[2]*(qR[0]-qL[0]);
	real bar_Du5=(qR[3]-qL[3])-((qR[2]-qL[2])-primRoe[2]*(qR[0]-qL[0]));
	alpha[1]=((Gamma-1)/pow(aroe,2))*((qR[0]-qL[0])*(Hroe-pow(primRoe[1],2))+primRoe[1]*(qR[1]-qL[1])-bar_Du5);
	alpha[0]=(1./(2.*aroe))*((qR[0]-qL[0])*(primRoe[1]+aroe)-(qR[1]-qL[1])-aroe*alpha[1]);
	alpha[3]=(qR[0]-qL[0])-(alpha[0]+alpha[1]);

	// Compute intercell flux using wavespeed and Roe eigenvectors computed above - use (11.29)
	eulerfluxes(primL,fl);
	eulerfluxes(primR,fr);
	for (int k=0;k<nvar;k++) {
		real summ=0.;
		for (int kk=0;kk<nvar;kk++) {
			summ+=alpha[kk]*abs(lambda[kk])*Kroe[k][kk];
		}
		ff[k]=0.5*(fl[k]+fr[k])-0.5*summ;
	}
}


// Exact Riemann Solver
void prim2exact_hydro(real (&primL)[nvar],real (&primR)[nvar],real (&ff)[nvar]) {
	real pstar, ustar;

	// Correct negative rho & p at cell edges
	primL[0]= primL[0]<=0. ? eps : primL[0];
	primL[3]= primL[3]<=0. ? eps : primL[3];
	primR[0]= primR[0]<=0. ? eps : primR[0];
	primR[3]= primR[3]<=0. ? eps : primR[3];
	// Compute sound speed and Mach number
        real al=sqrt(Gamma*primL[3]/primL[0]);
        real ar=sqrt(Gamma*primR[3]/primR[0]);
	real Ml=primL[1]/al;
	real Mr=primR[1]/ar;	
	
	real primf[nvar]; primf[2]=0.;
	real af;
	if (Ml>=1.0 and Mr>=1.0) { // Supersonic flow to right
		for (int k=0;k<nvar;k++) {primf[k]=primL[k];}
		af=al;
	} else if (Ml<=-1.0 and Mr<=-1.0) { // Supersonic flow to left
		for (int k=0;k<nvar;k++) {primf[k]=primR[k];}
		af=ar;
	} else {		
		// Check for cavity--------------------------------
		if ((primL[1]+Gamma6*al)<(primR[1]-Gamma6*ar)) {
			cout<<"Cavity occurs! Program stopped!"<<endl;
			throw exception();
		}
		// Call root-finding method for pstar--------------
		Newton(primL,primR,al,ar,pstar,ustar);
		// Calculate final states--------------------------	
		// (a) Density
		//     (i) Left running wave
		real dLstar,dRstar;
		if (pstar>=primL[3]) {
			dLstar=primL[0]*(1.+Gamma9*pstar/primL[3])/(Gamma9+pstar/primL[3]);
		} else if (pstar<primL[3]) {
			dLstar=primL[0]*pow(pstar/primL[3],1./Gamma);
		}
		//     (ii) Right running wave
		if (pstar>=primR[3]) {
			dRstar=primR[0]*(1.+Gamma9*pstar/primR[3])/(Gamma9+pstar/primR[3]);
		} else if (pstar<primR[3]) {
			dRstar=primR[0]*pow(pstar/primR[3],1./Gamma);
		}
		// (b) Speed of sound
		real aLstar=sqrt(Gamma*pstar/dLstar);
		real aRstar=sqrt(Gamma*pstar/dRstar);
		// Determine characteristics-----------------------
		real G1=primL[1]-al;
		real G2=ustar-aLstar;
		real G3=ustar+aRstar;
		real G4=primR[1]+ar;
		real G0=ustar;                   
		// Determine wave speeds---------------------------
		real w0,w1,w2,w3,w4;
		// (a) Entropy wave
		w0=G0;
		// (b) Left running wave
		if (pstar>=primL[3]) {
			w1=primL[1]-al*sqrt(1.+Gamma3*(pstar/primL[3]-1.));
			w2=w1;
		} else if (pstar<primL[3]) {
			w1=G1; w2=G2;
		}
		// (c) Right running wave
		if (pstar>=primR[3]) {
			w4=primR[1]+ar*sqrt(1.+Gamma3*(pstar/primR[1]-1.));
			w3=w4;
		} else if (pstar<primR[3]) {
			w4=G4; w3=G3;
		}
		// Wave vector:
		vector<real> w={w1,w2,w0,w3,w4};
		// Determine conditions on time axis---------------
		// (a) Wave closest to time axis
		int N, K = minloc(w,5);
		// Numbering of areas
		// inf|w1 -> N=0 LEFT INITIAL STATE
	        // w1|w2 -> N=1 LEFT EXPANSION FAN
	        // w2|w0 -> N=2 LEFT FINAL STATE
	        // w0|w3 -> N=3 RIGHT FINAL STATE
	        // w3|w4 -> N=4 RIGHT EXPANSIION FAN
	        // w4|inf -> N=5 RIGHT INITIAL STATE
	        // Region of time axis
		if (w[K]>=0) {
			N=K;
		} else if (w[K]<0 and w[K]!=w[K+1]) {
			N=K+1;
		} else if (w[K]<0 and w[K]==w[K+1]) {
			N=K+2;
		}
		// Determing final conditions on time axis
	       	switch (N) {
			case 0:
				primf[0]=primL[0];
				primf[1]=primL[1];
				primf[3]=primL[3];
				af=al;
				break;
			case 1:
				primf[1]=Gamma8*al+Gamma10*primL[1];
				af=Gamma8*al+Gamma10*primL[1];
				primf[3]=primL[3]*pow(af/al,Gamma2);
				primf[0]=primL[0]*pow(af/al,Gamma6);
				break;
			case 2:
				primf[0]=dLstar;
				primf[1]=ustar;
				primf[3]=pstar;
				af=aLstar;
				break;
			case 3:
				primf[0]=dRstar;
				primf[1]=ustar;
				primf[3]=pstar;
				af=aRstar;
				break;
			case 4:
				primf[1]=-Gamma8*ar+Gamma10*primR[1];
				af=Gamma8*ar-Gamma10*primR[1];
				primf[3]=primR[3]*pow(af/ar,Gamma2);
				primf[0]=primR[0]*pow(af/ar,Gamma6);
				break;
			case 5:
				primf[0]=primR[0];
				primf[1]=primR[1];
				primf[3]=primR[3];
				af=ar;
				break;
		}
	}
	// Calculate fluxes on cell faces------------------
	eulerfluxes(primf,ff);
}

void Newton(real (&primL)[nvar],real (&primR)[nvar],real aL,real aR,real &pstar,real &ustar) {
	// Iteration parameters
	int iter=20, QUSER=2.0;
	// Flow paramters
	pstar=0.0; ustar=0.0;
	real CUP=0.25*(primL[0]+primR[0])*(aL+aR);
	real PPV=0.5*(primL[3]+primR[3])+0.5*(primL[1]-primR[1])*CUP;
	PPV=MAX(eps,PPV);
	real pMIN=MIN(primL[3],primR[3]);
	real pMAX=MAX(primL[3],primR[3]);
	real qMAX=pMAX/pMIN;
	real pStart,pQ,uM,ptL,ptR;
	if (qMAX<=QUSER and (pMIN<=PPV and PPV<=pMAX)) {
		pStart=PPV;
	} else {
		if (PPV<pMIN) {
			pQ=pow(primL[3]/primR[3],Gamma1);
			uM=(pQ*primL[1]/aL+primR[1]/aR+(1./Gamma5)*(pQ-1.))/(pQ/aL+1./aR);
			ptL=1.+Gamma5*(primL[1]-uM)/aL;
			ptR=1.+Gamma5*(uM-primR[1])/aR;
			pStart=0.5*(primL[3]*pow(ptL,1./Gamma1)+primR[3]*pow(ptR,1./Gamma1));
		} else {
			real GEL=sqrt((Gamma8/primL[0])/(Gamma10*primL[3]+PPV));
			real GER=sqrt((Gamma8/primR[0])/(Gamma10*primR[3]+PPV));
			pStart=(GEL*primL[3]+GER*primR[3]-(primR[1]-primL[1]))/(GEL+GER);
		}
	}
	real pOld=pStart;
	real Udiff=primR[1]-primL[1];
	real Change=0.0, emax=1.E-6;
	real FL,FLD,FR,FRD;
	for (int i=0; i<iter; i++) {
		PREFUN(FL,FLD,pOld,primL[0],primL[3],aL);
		PREFUN(FR,FRD,pOld,primR[0],primR[3],aR);
		pstar=pOld-(FL+FR+Udiff)/(FLD+FRD);
		Change=2.0*abs((pstar-pOld)/(pstar+pOld));
		if (i>=iter-1 and Change>emax) {
			cout<<"Unable to converge during Newton iteration! Pressure change at "<<Change<<endl;
			throw exception();
		}
		if (pstar<0.0) {
			cout<<"Negative pressure of "<<pstar<<" encountered!"<<endl;
			pstar=eps;
		}
		pOld=pstar;		
	}
	ustar=0.5*(primL[1]+primR[1]+FR-FL);
}


void PREFUN(real &F,real &FD,real p,real dk,real pk,real ak) {
	if (p<=pk) {
		real PRAT=p/pk;
		F=Gamma6*ak*(pow(PRAT,Gamma1)-1.);
		FD=(1.0/(dk*ak))*pow(PRAT,-Gamma3);		
	}else {
		real AK=Gamma8/dk;
		real BK=Gamma10*pk;
		real QRT=sqrt(AK/(BK+p));
		F=(p-pk)*QRT;
		FD=(1.-0.5*(p-pk)/(BK+p))*QRT;
	}
}
