#include "defs.hpp"
int getLeft(meshblock* dom, const int nb);
int getRight(meshblock* dom, const int nb);
int getlevel(meshblock* dom, const int nb);

void setbcs(meshblock* dom,real (&u)[nvar][nx+2*nghosts][nbmax]) {
	int nbL,nbR;

	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			// Open boundary conditions
			// Left
			nbL=getLeft(dom,nb);
			for (int ii=0;ii<nvar;ii++) {
				if (nbL==-1) {
					for (int ng=0; ng<nghosts; ng++) {
						u[ii][ng][nb]=u[ii][dom->nxmin+nghosts-ng-1][nb];
					}
				} else {
					for (int ng=0; ng<nghosts; ng++) {
						int nblev=getlevel(dom,nb);
						int nbLlev=getlevel(dom,nbL);
						if (abs(nblev-nbLlev)>1) {cout<<"refine diff>1!"<<endl; throw exception();}
						if (nblev==nbLlev) {
							u[ii][ng][nb]=u[ii][dom->nxmax-nghosts+ng+1][nbL];
						} else if (nblev>nbLlev) { // More refined than left
							int ng2=dom->nxmax-(nghosts-ng-1)/2;
							u[ii][ng][nb]=u[ii][ng2][nbL];
						} else if (nblev<nbLlev) { // Coarser than left
							int ng2=dom->nxmax-2*(nghosts-ng-1);
							u[ii][ng][nb]=0.5*(u[ii][ng2-1][nbL]+u[ii][ng2][nbL]);
						}	
					}
				}
			}

			// Right
			nbR=getRight(dom,nb);
			for (int ii=0;ii<nvar;ii++) {
				if (nbR==-1) {
					for (int ng=0; ng<nghosts; ng++) {
						u[ii][dom->nxmaxb+ng][nb]=u[ii][dom->nxmax-ng][nb];
					}
				} else {
					for (int ng=0; ng<nghosts; ng++) {
						int nblev=getlevel(dom,nb);
						int nbRlev=getlevel(dom,nbR);						
						if (abs(nblev-nbRlev)>1) {cout<<"refine diff>1!"<<endl; throw exception();}
						if (nblev==nbRlev) {
							u[ii][dom->nxmaxb+ng][nb]=u[ii][dom->nxmin+ng][nbR];
						} else if (nblev>nbRlev) { // More refined than right
                                                        int ng2=dom->nxmin+ng/2;
                                                        u[ii][dom->nxmaxb+ng][nb]=u[ii][ng2][nbR];
                                                } else if (nblev<nbRlev) { // Coarser than right
                                                        int ng2=dom->nxmin+2*ng;
                                                        u[ii][dom->nxmaxb+ng][nb]=0.5*(u[ii][ng2][nbR]+u[ii][ng2+1][nbR]);
                                                }
					}
				}
			}
		} 
	}
}
