#include<math.h>
#include<iostream>
#include "defs.hpp"
int getlevel(meshblock &dom,int nb);
int get_nb(meshblock &dom,int nb);
real getBlockPosition(meshblock &dom, const int nb);
void output(meshblock &dom,string ictype,const int itprint);
// Mark (set flag) blocks for refinement/coarsening based on then gradients of density & pressure
void FlagGrads(meshblock &dom) {
	int level;
	real gradP, gradrho, maxgrad=0.;

	for (int nb=0;nb<dom.lastActive;nb++) {
		if (dom.ActiveBlocks[nb]!=-1) {
			maxgrad=0.;
			level=getlevel(dom,nb);
			for (int i=1;i<=nx;i++) {
				gradrho=(dom.prim[0][i+1][nb]-dom.prim[0][i-1][nb])/(dom.prim[0][i][nb]*2.*dom.dx[level]);
				gradP=(dom.prim[3][i+1][nb]-dom.prim[3][i-1][nb])/(dom.prim[3][i][nb]*2.*dom.dx[level]);

				maxgrad=max(maxgrad,abs(gradrho));
				maxgrad=max(maxgrad,abs(gradP));
				if (MAG_FIELD_ENABLED) {
					real gradBx=(dom.prim[4][i+1][nb]-dom.prim[4][i-1][nb])/(dom.prim[4][i][nb]*2.*dom.dx[level]);
					real gradBy=(dom.prim[5][i+1][nb]-dom.prim[5][i-1][nb])/(dom.prim[5][i][nb]*2.*dom.dx[level]);
					real gradB=max(abs(gradBx),abs(gradBy));
					maxgrad=max(maxgrad,gradB);
				}
				if (maxgrad>=rThresh) {
					// Mark for refinement (only if it is not at max resolution alr) and stop checking
					if (level<nlevs-1) dom.FlagRefine[nb]=true;
					break;
				}
			}

			if (maxgrad<cThresh) dom.FlagCoarse[nb]=true;

		}
	}
}

// Mark blocks for refinement based on proximity criteria
void FlagProx(meshblock &dom) {
	int dadID, nbLeft, nbRight, myID;

	for (int ilevs=nlevs-1;ilevs>=0;ilevs--) {
		for (int nb=0;nb<dom.lastActive;nb++) {
			if (dom.ActiveBlocks[nb]!=-1) {
				// Check all blocks marked for refinement and mark the neighbour
				// for refinement (and inhibit coarsening) if it is at a lower
				// level, just inhibit refinement if it is at same level
				if (dom.FlagRefine[nb]) {
					dadID=dom.ActiveBlocks[nb];
					// (a) Same level neighbours
					//     left
					nbLeft=get_nb(dom,dadID-1);
					if (nbLeft!=-1) dom.FlagCoarse[nbLeft]=false;
					//     right
					nbRight=get_nb(dom,dadID+1);
					if (nbRight!=-1) dom.FlagCoarse[nbRight]=false;
					// (b) Lower level neighbours
					//     left
					nbLeft=get_nb(dom,dadID/2-1);
					if (nbLeft!=-1) {
						dom.FlagRefine[nbLeft]=true;
						dom.FlagCoarse[nbLeft]=false;
					}
					//      right
					nbRight=get_nb(dom,dadID/2+1);
					if (nbRight!=-1) {
						dom.FlagRefine[nbRight]=true;
						dom.FlagCoarse[nbRight]=false;
					}
				}

				if (dom.FlagCoarse[nb]) {
					myID=dom.ActiveBlocks[nb];
					// (a) Inhibit coarsening if finer neighbours are not set for coarsening finer neighbours
					nbLeft=get_nb(dom,myID*2-1);
					if (nbLeft!=-1) {
						if (!dom.FlagCoarse[nbLeft]) dom.FlagCoarse[nb]=false;
					}
					//      right
					nbRight=get_nb(dom,myID*2+2);
					if (nbRight!=-1) {
						if (!dom.FlagCoarse[nbRight]) dom.FlagCoarse[nb]=false;
					}

					// (b) Inhibit coarsening if sibling is not set for coarsening
					if (myID%2==0) { // if block is a 1st son
						nbRight=get_nb(dom,myID+1);
						if (!dom.FlagCoarse[nbRight]) dom.FlagCoarse[nb]=false;
					} else {
						nbLeft=get_nb(dom,myID-1);
						if (!dom.FlagCoarse[nbLeft]) dom.FlagCoarse[nb]=false;
					}
				}
			}
		}
	}
}
