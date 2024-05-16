#include<math.h>
#include<iostream>
#include<exception>
#include<iomanip>
#include "defs.hpp"

real u2prim(const real Gamma,const int i,const real uu[nvar]);
void FlagGrads(meshblock &dom);
void FlagProx(meshblock &dom);


// Initializes all things mesh------------------------------------------------------------------------
// For the ID, (if 7 levels are used)
// Level		minID                                 maxID
// 0      		1                                     1
// 1      		2                                     3
// 2      		4                                     7
// 3      		8                                     15
// 4      		16                                    31
// 5      		32                                    63
// 6      		64                                    127
// Note that minID is where x=0 and maxID is where x=1. Thus, if bID=4, level is 2+1 while xb=0 & x=0.
void init_mesh(meshblock &dom) {
	// Compute dx
	dom.dx[0]=xmax/nx;
	for (int i=1;i<nlevs;i++) {
		dom.dx[i]=dom.dx[i-1]/2.;
	}
	// Min and Max block ID per level
	for (int i=0;i<nlevs;i++) {
		dom.minID[i]=pow(2,i);    
		dom.maxID[i]=pow(2,i+1)-1;
	}
	for (int nb=0;nb<nbmax;nb++) {
		dom.FlagRefine[nb]=false;
		dom.FlagCoarse[nb]=false;
		dom.ActiveBlocks[nb]=-1;
	}
	dom.ActiveBlocks[0]=4 ;
	dom.ActiveBlocks[1]=10;
	dom.ActiveBlocks[2]=22;
        dom.ActiveBlocks[3]=23;
	dom.ActiveBlocks[4]=24;
        dom.ActiveBlocks[5]=25;
	dom.ActiveBlocks[6]=7 ;
        dom.ActiveBlocks[7]=13;
	dom.lastActive=8;	
}

// Obtains the level of a block with bID
int getlevel(meshblock &dom, const int nb) {
	int level=-1;
	int bID=dom.ActiveBlocks[nb];

	if (bID==-1) {
		level=-1;
	} else {
		for (int i=0;i<nlevs;i++) {
			if ((bID>=dom.minID[i]) && (bID<=dom.maxID[i])) {
				level=i; 
				return level;
			}
		}
	}
	return level;
}

// Obtains position of first cell in block
real getBlockPosition(meshblock &dom, const int nb) {
	real xb;
	int level;
	int bID=dom.ActiveBlocks[nb];
	if (bID==-1) {
		xb=-1.;
	} else {
		level=getlevel(dom,nb);
		xb=(bID-dom.minID[level])*nx*dom.dx[level];
	}
	return xb;
}

// Obtains nb from bID 
int get_nb(meshblock &dom, const int bID) {
	int nb=-1;
	for (int i=0; i<dom.lastActive; i++) {
		if (dom.ActiveBlocks[i]==bID) nb=i;
	}
	return nb;
}

// Obtains nb from the neighbour at the left
int getLeft(meshblock &dom, const int nb) {
	int nbLeft=-1;
	int selfID,fatherID,level,sonID;

	// Get bID of current block
	selfID=dom.ActiveBlocks[nb];

	// (a) Domain boundary
	level=getlevel(dom,nb);
	if (selfID==dom.minID[level]) return -1;
	// (b) Same level of refinement
	nbLeft=get_nb(dom,selfID-1);
	if (nbLeft!=-1) return nbLeft;
	// (c) Lower level of refinement
	fatherID=selfID/2;
	nbLeft=get_nb(dom,fatherID-1);
	if (nbLeft!=-1) return nbLeft;
	// (d) Higher level of refinement
	sonID=selfID*2;
	nbLeft=get_nb(dom,sonID-1);
	return nbLeft;
}

// Obtains nb from the neighbour at the right
int getRight(meshblock &dom, const int nb) {
	int nbRight=-1;
	int selfID,fatherID,level,sonID;

	// Get bID of current block
	selfID=dom.ActiveBlocks[nb];

	// (a) Domain boundary
        level=getlevel(dom,nb);
        if (selfID==dom.maxID[level]) return -1;
        // (b) Same level of refinement
        nbRight=get_nb(dom,selfID+1);
        if (nbRight!=-1) return nbRight;
        // (c) Lower level of refinement
        fatherID=selfID/2;
        nbRight=get_nb(dom,fatherID+1);
        if (nbRight!=-1) return nbRight;
        // (d) Higher level of refinement
        sonID=selfID*2;
        nbRight=get_nb(dom,sonID+2); 
        return nbRight;
}

// Refines block nb to next level, refinement is applied both to the primitives & conserved variables
void refineBlock(meshblock &dom, const int dadNb) {
	int dadID, son1ID, son2ID, son1nb, son2nb;

	dadID=dom.ActiveBlocks[dadNb];
	son1ID=dadID*2;
	son2ID=son1ID+1;
	// change the bID of the parent to that of the first son
	dom.ActiveBlocks[dadNb]=son1ID;
	son1nb=dadNb;
	// Activate the bID of the second son in the first available spot
	// (increase lastActive if needed)
	son2nb=-1;
	for (int i=0;i<nbmax;i++){
		if (dom.ActiveBlocks[i]==-1) {
			dom.ActiveBlocks[i]=son2ID;
			son2nb=i;		
			if (i>dom.lastActive-1) dom.lastActive=i+1;  
			break;
		}
	}
	if (son2nb==-1) {
		cout<<"No. of blocks used exceeded max no. of blocks allowed! Aborting..."<<endl;
		throw std::exception();
	}

	// Copy data from dad to sons and update primitives
	// |          dad          |
	// |    son1   |    son2   |
	// |      (b)<-|->(a)      |
	// (a) Second son (son2) - ascending order (->)
	real intmArray[nvar];
	for (int i=0;i<nvar;i++) {intmArray[i]=0.;}
	for (int i=0;i<=nx+1;i++) {
		real pp;
		for (int ii=0;ii<nvar;ii++) {
			dom.u[ii][i][son2nb]=dom.u[ii][(i+1+nx)/2][dadNb];
			intmArray[ii]=dom.u[ii][i][son2nb];
		}
		for (int ii=0;ii<nvar;ii++) {
			pp=u2prim(Gamma,ii,intmArray);
			dom.prim[ii][i][son2nb]=pp;
		}
	}
	// (b) First son (son1) - descending order (<-)
	for (int i=nx+1;i>=0;i--){
		real pp;
		for (int ii=0;ii<nvar;ii++) {
			dom.u[ii][i][son1nb]=dom.u[ii][(i+1)/2][dadNb];
			intmArray[ii]=dom.u[ii][i][son1nb];
		}
		for (int ii=0;ii<nvar;ii++) {
			pp=u2prim(Gamma,ii,intmArray);
			dom.prim[ii][i][son1nb]=pp;
		}
	}
	// Reset refining flag
	dom.FlagRefine[dadNb]=false;
}

// Coarsens block nb to previous level, coarsening is applied both to the primitives & conservatives variables
void coarseBlock(meshblock &dom, int son1nb) {
	int son2ID, son1ID, dadID, son2nb, dadNb;

	son1ID=dom.ActiveBlocks[son1nb];
	if (son1ID%2==0) { 
		son2ID=son1ID+1;
		son2nb=get_nb(dom,son2ID);
	} else {
		son1ID-=1;
		son2ID=son1ID+1;
		son2nb=son1nb;
		son1nb=get_nb(dom,son1ID);
	}
	dadID=son1ID/2;

	// Use the memory space of first son as target
	dadNb=son1nb;
	// (a) 1st half
	//     The ghost cell does not need to be modified.
	real intmArray[nvar];
	for (int i=0;i<nvar;i++) {intmArray[i]=0.;}
	for (int i=1;i<=nx/2;i++) {
		real pp;
		for (int ii=0;ii<nvar;ii++) {
			dom.u[ii][i][dadNb]=0.5*(dom.u[ii][2*i-1][son1nb]+dom.u[ii][2*i][son1nb]);
			intmArray[ii]=dom.u[ii][i][dadNb];
		}
		for (int ii=0;ii<nvar;ii++) {
			pp=u2prim(Gamma,ii,intmArray);
			dom.prim[ii][i][dadNb]=pp;
		}
	}
	// (b) 2nd half
	for (int i=1;i<=nx/2;i++) {
		real pp;
		for (int ii=0;ii<nvar;ii++) {
			dom.u[ii][i+nx/2][dadNb]=0.5*(dom.u[ii][2*i-1][son2nb]+dom.u[ii][2*i][son2nb]);
			intmArray[ii]=dom.u[ii][i+nx/2][dadNb];
		} 
		for (int ii=0;ii<nvar;ii++) {
			pp=u2prim(Gamma,ii,intmArray);
			dom.prim[ii][i+nx/2][dadNb]=pp;
		}
	}
	// (c) Ghost cell at the right
	real pp;
	for (int ii=0;ii<nvar;ii++) {
		dom.u[ii][nx+1][dadNb]=dom.u[ii][nx+1][son2nb];
		intmArray[ii]=dom.u[ii][nx+1][dadNb];
	}
	for (int ii=0;ii<nvar;ii++) {
		pp=u2prim(Gamma,ii,intmArray);
		dom.prim[ii][nx+1][dadNb]=pp;
	}

	// Update ActiveBlocks
	dom.ActiveBlocks[dadNb]=dadID;
	dom.ActiveBlocks[son2nb]=-1;

	// Update Coarsening Flags
	dom.FlagCoarse[son1nb]=false;
	dom.FlagCoarse[son2nb]=false;
}

// Updates (refines & coarsens) mesh
void update_mesh(meshblock &dom) {
	int width=10;
	// Mark by physical criteria
	FlagGrads(dom);
	// Mark by proximity
	FlagProx(dom);
	
	// Proceed with refinement of marked blocks
	for (int i=0;i<dom.lastActive;i++) {
		if (dom.ActiveBlocks[i]!=-1) {
			if (dom.FlagCoarse[i]) {
				cout<<setw(17)<<"Coarsening: bID="<<setw(width)<<dom.ActiveBlocks[i]<<", nb="<<setw(width)<<i<<", lastActive="<<setw(width)<<dom.lastActive<<endl;
				coarseBlock(dom,i);
			}
			if (dom.FlagRefine[i]) {
				cout<<setw(17)<<"Refining: bID="<<setw(width)<<dom.ActiveBlocks[i]<<", nb="<<setw(width)<<i<<", lastActive="<<setw(width)<<dom.lastActive<<endl;
				refineBlock(dom,i);
			}
		}
	}
}
