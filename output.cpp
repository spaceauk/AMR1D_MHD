#include<fstream>
#include<sstream>
#include<iomanip>   // setw()
#include "defs.hpp"

int getlevel(meshblock* dom,const int nb);
real getBlockPosition(meshblock* dom,const int nb);

// Output to file
void output(meshblock* dom,string ictype,const int itprint) {
	int level;
	real xb;

	// Open output file
	stringstream ss;
	ss<<setw(3)<<setfill('0')<<itprint;
	string id=ss.str();
	string fname="./data/"+ictype+"_t"+id+".dat";
	
	ofstream Wdata;
	int width=20;
	Wdata.open(fname);
	// Write x and rho, u & P
	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			level=getlevel(dom,nb);
			xb=getBlockPosition(dom,nb);	
			for (int i=dom->nxmin;i<=dom->nxmax;i++) {
				real x=xb+(i-(0.5+nghosts/2))*dom->dx[level];
				Wdata<<setw(width)<<x<<setw(width)
				<<dom->prim[0][i][nb]<<setw(width)
				<<dom->prim[1][i][nb]<<setw(width)
				<<dom->prim[3][i][nb];
				if (nvar>4) {
					Wdata<<setw(width)<<dom->prim[2][i][nb]<<
					setw(width)<<dom->prim[4][i][nb]<<
					setw(width)<<dom->prim[5][i][nb];	
				}
				Wdata<<endl;
			}
		}
	}
	Wdata.close();
}

void plotdata(RiemannInv* dom,const int itprint) {
	real xb,x;
	int level;
	stringstream ss;
	string fname="data/finalresult.dat";
	ofstream Wdata;
	int width=20;
	Wdata.open(fname);
	for (int nb=0;nb<dom->lastActive;nb++) {
		if (dom->ActiveBlocks[nb]!=-1) {
			level=getlevel(dom,nb);
			xb=getBlockPosition(dom,nb);
			for (int i=dom->nxmin;i<=dom->nxmax;i++) {
				x=xb+(i-(0.5+nghosts/2))*dom->dx[level];
				Wdata<<setw(width)<<x<<setw(width)
                	                <<dom->prim[0][i][nb]<<setw(width)
                        	        <<dom->prim[1][i][nb]<<setw(width)
                                	<<dom->prim[3][i][nb];
				if (nvar>4) {
                                        Wdata<<setw(width)<<dom->prim[2][i][nb]<<
					setw(width)<<dom->prim[4][i][nb]<<
                                        setw(width)<<dom->prim[5][i][nb];
                                }
                                Wdata<<endl;
			}
		}
	}
	Wdata.close();

	// Save data for Riemann invariant to be plotted
	if (!MAG_FIELD_ENABLED) {
		fname="data/char_curves.dat";
		Wdata.open(fname);
		for (int nt=0; nt<dom->numt; nt++) {
			for (int nb=0; nb<dom->lastActive; nb++) {
				if (dom->ActiveBlocks[nb]!=-1) {
				level=getlevel(dom,nb);
				xb=getBlockPosition(dom,nb);
				for (int i=dom->nxmin; i<=dom->nxmax; i++) {
					x=xb+(i-(0.5+nghosts/2))*dom->dx[level];
					Wdata<<setw(width)<<dom->tsave[nt]
					<<setw(width)<<x
   					<<setw(width)<<dom->Jplus[i][nb][nt]
					<<setw(width)<<dom->Jminus[i][nb][nt]
					<<setw(width)<<dom->s[i][nb][nt]<<endl;					
				}
			}
			}
		}
	}
	Wdata.close();

	// Construct the command to plot data
	string command="python3 plotdata.py ";
	cout<<command<<endl;
	system(command.c_str());

}
