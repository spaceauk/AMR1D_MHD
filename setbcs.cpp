#include<iostream>
#include "defs.hpp"
int getLeft(meshblock &dom, const int nb);
int getRight(meshblock &dom, const int nb);

void setbcs(meshblock &dom,real (&u)[nvar][nx+2][nbmax]) {
	int nbL,nbR;

	for (int nb=0;nb<dom.lastActive;nb++) {
		if (dom.ActiveBlocks[nb]!=-1) {
			// Open boundary conditions
			// Left
			nbL=getLeft(dom,nb);
			for (int ii=0;ii<nvar;ii++) {
				if (nbL==-1) {
					u[ii][0][nb]=u[ii][1][nb];
				} else {
					u[ii][0][nb]=u[ii][nx][nbL];
				}
			}

			// Right
			nbR=getRight(dom,nb);
			for (int ii=0;ii<nvar;ii++) {
				if (nbR==-1) {
					u[ii][nx+1][nb]=u[ii][nx][nb];
				} else {
					u[ii][nx+1][nb]=u[ii][1][nbR];
				}
			}
		} 
	}
}
