#include<math.h>
#include "defs.hpp"

int getlevel(meshblock* dom,int nb);
void setbcs(meshblock* dom,real (&u)[nvar][nx+2*nghosts][nbmax]);
void fluxSolver(meshblock* dom,real (&f)[nvar][nx+2*nghosts][nbmax]);
void update_prim(meshblock* dom);

// Compute the timestep allowed by the CFL
real timestep(meshblock* dom) {
	real dt;
	// Courant number = 0.9
	real cfl=0.9;
	real del=1.e30, cs;
	int level;
	
	// Sweep all active blocks
	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			level=getlevel(dom,nb);
			for (int i=0;i<=dom->ntot;i++) {
				// Sound speed
				cs=sqrt(Gamma*dom->prim[3][i][nb]/dom->prim[0][i][nb]);
				if (nvar>4) {
					// Speed of Alven waves
					real ca=sqrt((pow(dom->prim[4][i][nb],2)+pow(dom->prim[5][i][nb],2))/dom->prim[0][i][nb]);
					real cax=sqrt(pow(dom->prim[4][i][nb],2)/dom->prim[0][i][nb]);
					real cay=sqrt(pow(dom->prim[5][i][nb],2)/dom->prim[0][i][nb]);
					// Speed of magnetosonic waves
					real cfx=sqrt(0.5*(cs*cs+ca*ca)+0.5*sqrt(pow(cs*cs+ca*ca,2)-4*(cs*cs*cax*cax)));
					real cfy=sqrt(0.5*(cs*cs+ca*ca)+0.5*sqrt(pow(cs*cs+ca*ca,2)-4*(cs*cs*cay*cay)));
					cs=max(abs(dom->prim[1][i][nb])+cfx,abs(dom->prim[2][i][nb])+cfy);
				} else {cs+=abs(dom->prim[1][i][nb]);}
				del=min(del,dom->dx[level]/cs);
			}	
		}
	}
	dt=cfl*del;
	return dt;
}

// Integrate from t to t+dt with the method of Lax
void tstep(meshblock* dom,const real dt,const real time) {
	real dtx;
	int level;
	real f[nvar][nx+2*nghosts][nbmax];
	
	// Obtain fluxes
	fluxSolver(dom,f);

	// RK2 1st timestep	
	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			level=getlevel(dom,nb);
			dtx=dt/dom->dx[level];
			for (int i=0;i<=dom->ntot;i++){
				for (int ii=0;ii<nvar;ii++){
					dom->up[ii][i][nb]=dom->u[ii][i][nb];
					if (i>=dom->nxmin and i<=dom->nxmax) {
						dom->u[ii][i][nb]=dom->u[ii][i][nb]-dtx*(f[ii][i][nb]-f[ii][i-1][nb]);
					}
				}
			}
		}	
	}
	
	// Boundary conditions to the U^n+1
	setbcs(dom,dom->u);
	
	// Find updated flux
	update_prim(dom);
	fluxSolver(dom,f);

	// RK2 2nd timestep
	for (int nb=0;nb<dom->lastActive;nb++) {
                if (dom->ActiveBlocks[nb]!=-1) {
                        level=getlevel(dom,nb);
                        dtx=dt/dom->dx[level];
                        for (int i=dom->nxmin;i<=dom->nxmax;i++){
                                for (int ii=0;ii<nvar;ii++){
                                        dom->u[ii][i][nb]=0.5*(dom->u[ii][i][nb]+dom->up[ii][i][nb]-dtx*(f[ii][i][nb]-f[ii][i-1][nb]));
                                }
                        }
                }
        }

	// Boundary conditions to the U^n+1
        setbcs(dom,dom->u);

}
