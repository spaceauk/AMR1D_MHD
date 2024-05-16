#include<math.h>
#include<string>
#include<iostream>    // std::cout
#include<exception>   // call error exception
#include<fstream>     // input and output
#include<iomanip>     // std::setw		 
#include "defs.hpp"
#include "mesh.hpp"

void initconds(meshblock &dom,real time,real tprint,real itprint,string ictype);
void update_prim(meshblock &dom);
void output(meshblock &dom,string ictype,const int itprint);
real timestep(meshblock &dom);
void tstep(meshblock &dom,const real dt,const real time);
void plotdata(meshblock &dom, const int itprint);


int main() {
	real             time=0., dt=0.;
	real             tmax;
	const real       dtprint=0.01;
	real             tprint=0.;
	int              itprint=0;
	string           ictype="TC1";
	meshblock dom1;

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
	for (int nb=0;nb<dom1.lastActive;nb++) {
		real xb=getBlockPosition(dom1, nb);
		if (dom1.ActiveBlocks[nb]!=-1) {
			cout<<xb<<"; ";
			for (int i=0;i<nx;i++) {
				chk_mesh+=1;
			}
		}
	}
	cout<<endl<<"Initial meshsize="<<chk_mesh<<endl;
	// Updates the primitives
	update_prim(dom1);
	
	// Main loop - iterate until maximum time is reached
	int count=1;
	while (time<tmax and count<1500) {		
		// Output at tprint intervals
		if (time>=tprint) {
			cout<<setw(3)<<count<<") t="<<time<<", tmax="<<tmax<<", dt="<<dt<<", itprint="<<itprint<<endl;
			output(dom1,ictype,itprint);
			tprint+=dtprint;
			itprint+=1;
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
