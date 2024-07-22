#include<string>
#include<exception>   // call error exception
#include<fstream>     // input and output
#include<iomanip>     // std::setw		 
#include "defs.hpp"
#include "mesh.hpp"

void initconds(meshblock* dom,real time,real tprint,real itprint,string ictype);
void update_prim(meshblock* dom);
void output(meshblock* dom,string ictype,const int itprint);
real timestep(meshblock* dom);
void tstep(meshblock* dom,const real dt,const real time);
void plotdata(RiemannInv* dom, const int itprint);
real prim2u(const real Gamma,const int i,const real pp[nvar]);

int main() {
	real             time=0., dt=0.;
	real             tmax;
	const real       dtprint=0.01;
	real             tprint=0.;
	int              itprint=0;
	string           ictype="TC1";
	RiemannInv* dom1 = new RiemannInv();
	if (MAG_FIELD_ENABLED) {
		delete dom1;
		meshblock* dom1 = new meshblock();
	}
	dom1->solvertype();

	if (ictype=="TC1") {
		tmax=0.20;
	} else if (ictype=="TC2") {
                tmax=0.15;
        } else if (ictype=="TC3") {
                tmax=0.012;
        } else if (ictype=="TC4") {
                tmax=0.035;
        } else if (ictype=="TC5") {
                tmax=0.035;
        }

	// Initializes AMR
	init_mesh(dom1);

	// Generates initial conditions
	initconds(dom1,time,tprint,itprint,ictype);

	int chk_mesh=0;
	cout<<"Block position: ";
	for (int nb=0;nb<dom1->lastActive;nb++) {
		real xb=getBlockPosition(dom1, nb);
		if (dom1->ActiveBlocks[nb]!=-1) {
			cout<<xb<<"; ";
			for (int i=0;i<nx;i++) {
				chk_mesh+=1;
			}
		}
	}
	cout<<endl<<"Initial meshsize="<<chk_mesh<<endl;
	// Updates the conservative quantities
	real tmpprim[nvar];
	for (int nb=0; nb<dom1->lastActive; nb++) {
                if (dom1->ActiveBlocks[nb]!=-1) {
                        for (int i=0; i<=nx+1; i++) {
				for (int k=0; k<nvar; k++) {tmpprim[k]=dom1->prim[k][i][nb];}
				for (int k=0; k<nvar; k++) {dom1->u[k][i][nb]=prim2u(Gamma,k,tmpprim);}
			}
		}
	}
	
	// Main loop - iterate until maximum time is reached
	int count=0;
	while (time<tmax and count<maxiter) {		
		// Output at tprint intervals
		if (time>=tprint) {
			cout<<setw(3)<<count<<") t="<<time<<", tmax="<<tmax<<", dt="<<dt<<", itprint="<<itprint<<endl;
			output(dom1,ictype,itprint);
			tprint+=dtprint;
			itprint+=1;
		}

		// Calculate Riemann invariants ($u \pm 2*a/(\gamma-1)$ here only valid for hydro + isentropic)
                if (!MAG_FIELD_ENABLED) {
        	        dom1->numt=count+1;
                	dom1->tsave[count]=time;
        	        for (int nb=0; nb<dom1->lastActive; nb++) {
	                        if (dom1->ActiveBlocks[nb]!=-1) {
                	                for (int i=0; i<nx+2; i++) {
        	                                real a=sqrt(Gamma*dom1->prim[2][i][nb]/dom1->prim[0][i][nb]);
						if (dom1->prim[2][i][nb]<0. or dom1->prim[0][i][nb]<0.) {
							cout<<"Negative pressure or density! P="<<dom1->prim[2][i][nb]
							    <<", \rho="<<dom1->prim[0][i][nb]<<" @ i="<<i<<", nb="<<nb<<endl;}	
	                                        dom1->Jplus[i][nb][count]=dom1->prim[1][i][nb]+2*a/(Gamma-1.);
                                        	dom1->Jminus[i][nb][count]=dom1->prim[1][i][nb]-2*a/(Gamma-1.);
						dom1->s[i][nb][count]=log(dom1->prim[2][i][nb]/pow(dom1->prim[0][i][nb],Gamma));
                                	}
                        	}
                	}
		}

		// Obtain the $\delta t$ allowed by the CFL criterium
		dt=timestep(dom1);

		// Integrate u from t to t+dt
		tstep(dom1,dt,time);

		// Updates the primitives
		update_prim(dom1);
		
		// Updates mesh (b4 advancing the solution at every time step, mesh is updated first)
		update_mesh(dom1);

		// Time counter increases
		time+=dt;
		count+=1;	
	}

	// Plot data
	plotdata(dom1,itprint);
}
