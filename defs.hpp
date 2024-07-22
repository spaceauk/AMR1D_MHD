#include<math.h>
#include<iostream>

#define eps 1.0E-8
using real=float;
using namespace std;

// Global variables
const bool MAG_FIELD_ENABLED=false;
const int nvar= MAG_FIELD_ENABLED? 7 : 4; 
const int nx=5;
const int nlevs=7;                               // Base level is at 4 which mesh is initiated thus, cannot go below 4.
const int nbmax=pow(2,nlevs-1)+1;                // Initial # of blocks used here is 8.
const int maxiter = 1000;				  
const int ntmax= MAG_FIELD_ENABLED? 1 : maxiter; // Number of timesteps to be saved for plotting Riemann invariants.

const real xmax=1.;
const real Gamma=2.0;

const real rThresh=0.5, cThresh=0.05;

// Declarations for class file (domain setup)
class meshblock {
	public:
		real u[nvar][nx+2][nbmax];
		real prim[nvar][nx+2][nbmax];
		real dx[nlevs];
		int minID[nlevs], maxID[nlevs], ActiveBlocks[nbmax];
		int lastActive=1; 
		bool FlagRefine[nbmax], FlagCoarse[nbmax];
		virtual void solvertype() {
			cout<<"MHD solver used! Includes Bx & By."<<endl;
		}
};

class RiemannInv: public meshblock {
	public:
		// Include Riemann invariants for contour plot (only for Hydro)
                int numt;
                real tsave[ntmax];
                real Jplus[nx+2][nbmax][ntmax], Jminus[nx+2][nbmax][ntmax];
                real s[nx+2][nbmax][ntmax]; // Entropy
		void solvertype() override {
			cout<<"Hydro solver used! Includes Riemann invariant for contour plotting."<<endl;
		}
};
