#include<exception>
#include<string>
#include<iostream>
#include<math.h>
#include "defs.hpp"
int getlevel(meshblock &dom, const int nb);
real getBlockPosition(meshblock &dom, const int nb);

void initconds(meshblock &dom,real time,real tprint,real itprint,string ictype) {
	int level;
	real xb;
	real xp=0.5;

	// Sweep all ActiveBlocks
	for (int nb=0;nb<dom.lastActive;nb++) {
		// Get the level and position of first cell in each block
		if (dom.ActiveBlocks[nb]!=-1) {
		level=getlevel(dom,nb);
		xb=getBlockPosition(dom,nb);		
		for (int i=0;i<=nx+1;i++) {
			real x=xb+(i-0.5)*dom.dx[level];
			dom.u[2][i][nb]=0.; // For 1D, this is usually zero.
			if (ictype=="TC1") {
				if (i==0 && nb==0) {
					MAG_FIELD_ENABLED ? cout<<"(MHD) " : cout<<"(Hydro) "; 
					cout<<"Test case 1: Shock tube"<<endl;
				}
				if (x<xp) {
					dom.u[0][i][nb]=1.0;
					dom.u[1][i][nb]=0.0;
					dom.u[3][i][nb]=1.0/(Gamma-1.);
					if (MAG_FIELD_ENABLED) dom.u[5][i][nb]=1.0;					
				} else {
					dom.u[0][i][nb]=0.125;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=0.1/(Gamma-1.);
					if (MAG_FIELD_ENABLED) dom.u[5][i][nb]=-1.0;
				}
				if (MAG_FIELD_ENABLED) {
					dom.u[4][i][nb]=0.75; 
					dom.u[3][i][nb]+=0.5*(pow(dom.u[4][i][nb],2)+pow(dom.u[5][i][nb],2));
				}				
			} else if (ictype=="TC2") {
                                if (i==0 && nb==0) cout<<"(Hydro) Test case 2: Double rarefaction waves"<<endl;
                                if (x<xp) {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=2.0;
                                        dom.u[3][i][nb]=1.4/(Gamma-1.);
                                } else {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=-2.0;
                                        dom.u[3][i][nb]=0.4/(Gamma-1.);
                                }
				if (nvar>4) {dom.u[4][i][nb]=0.; dom.u[5][i][nb]=0.;}
			} else if (ictype=="TC3") {
                                if (i==0 && nb==0) cout<<"(Hydro) Test case 3: Ostrowski"<<endl;
                                if (x<xp) {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=1000.0/(Gamma-1.);
                                } else {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=0.01/(Gamma-1.);
                                }
				if (nvar>4) {dom.u[4][i][nb]=0.; dom.u[5][i][nb]=0.;}
			} else if (ictype=="TC4") {
                                if (i==0 && nb==0) cout<<"(Hydro) Test case 4: Ostrowski"<<endl;
                                if (x<xp) {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=0.01/(Gamma-1.);
                                } else {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=100.0/(Gamma-1.);
                                }
				if (nvar>4) {dom.u[4][i][nb]=0.; dom.u[5][i][nb]=0.;}
			} else if (ictype=="TC5") {
                                if (i==0 && nb==0) cout<<"(Hydro) Test case 5: Ostrowski"<<endl;
                                if (x<xp) {
                                        dom.u[0][i][nb]=5.99924;
                                        dom.u[1][i][nb]=19.5975;
                                        dom.u[3][i][nb]=460.0950/(Gamma-1.);
                                } else {
                                        dom.u[0][i][nb]=5.99242;
                                        dom.u[1][i][nb]=-6.19633;
                                        dom.u[3][i][nb]=46.0950/(Gamma-1.);
                                }
				if (nvar>4) {dom.u[4][i][nb]=0.; dom.u[5][i][nb]=0.;}
			} else {
				cout<<"Wrong testcase selected! Ending program..."<<endl;
				throw exception();
			}

			if ((x-0.5*dom.dx[level]<xp) && (x+0.5*dom.dx[level]>xp)){
				if (ictype=="TC1") {
					dom.u[0][i][nb]=1.125/2.;
					dom.u[1][i][nb]=0.0;
					dom.u[3][i][nb]=1.1/2./(Gamma-1.);
					if (MAG_FIELD_ENABLED) {
						dom.u[5][i][nb]=0.;
						dom.u[3][i][nb]+=0.5*(pow(dom.u[4][i][nb],2)+pow(dom.u[5][i][nb],2));
					}
				} else if (ictype=="TC2") {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=0.4/(Gamma-1.);
                                } else if (ictype=="TC3") {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=1000.01/2./(Gamma-1.);
                                } else if (ictype=="TC4") {
                                        dom.u[0][i][nb]=1.0;
                                        dom.u[1][i][nb]=0.0;
                                        dom.u[3][i][nb]=100.01/2./(Gamma-1.);
                                } else if (ictype=="TC5") {
                                        dom.u[0][i][nb]=5.99583;
                                        dom.u[1][i][nb]=6.700585;
                                        dom.u[3][i][nb]=506.989/2./(Gamma-1.);
                                }
			}
		}
		}
	}

	// Reset the counters and time to 0
	time=0.;
	tprint=0.;
	itprint=0.;
}
