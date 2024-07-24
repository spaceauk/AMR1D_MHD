#include<math.h>
#include<iostream>

#define eps 1.0E-8
using real=float;
using namespace std;

// Essential functions
#define MIN(a,b) ((a>=b)? b : a)
#define MAX(a,b) ((a>=b)? a : b)

// Global variables
const bool MAG_FIELD_ENABLED=true;
const int nvar= MAG_FIELD_ENABLED? 7 : 4; 
const int nx=6;                                  // Use even #
const int nlevs=7;                               // Base level is at 4 which mesh is initiated thus, cannot go below 4.
const int nbmax=pow(2,nlevs-1)+1;                // Initial # of blocks used here is 8.
const int maxiter = 500;				  
const int ntmax= MAG_FIELD_ENABLED? 1:maxiter/2; // Number of timesteps to be saved for plotting Riemann invariants.
const int nghosts=2;                             // If 4 ghost cells used, nx must be at least 8 to avoid using boundary of refined cells.

const real xmax=1.;
const real Gamma=2.0;

const real rThresh=0.5, cThresh=0.05;

// Declarations for class file (domain setup)
class meshblock {
	private:
		string RStype_MHD[1]={"Rusanov"};
		string SLname[8]={"1st","No SL","minmod","VanLeer","Superbee","MC","Koren","Cada3rd"};
	public:
		real u[nvar][nx+2*nghosts][nbmax];
		real up[nvar][nx+2*nghosts][nbmax];
		real prim[nvar][nx+2*nghosts][nbmax];
		real dx[nlevs];
		int minID[nlevs], maxID[nlevs], ActiveBlocks[nbmax];
		int lastActive=1; 
		bool FlagRefine[nbmax], FlagCoarse[nbmax];
		// Parameters
		int ntot, nxmin, nxmax, nxminb, nxmaxb, nx2;
		// Slope limiter type
                int SLtype=5; // 0) Godunov 1st-order; 1) No SL; 2)...
		// Riemann Solver type 
		int RStype=0; // MHD: 0) Rusanov;  
		
		void setParams();
		virtual void solvertype() {			
			cout<<"MHD Riemann solver ("<<RStype_MHD[RStype]<<") used! Includes Bx & By."<<endl;
			cout<<"-> Slope limiter used is "<<SLname[SLtype]<<endl;
		}
};

class RiemannInv: public meshblock {
	private:
		string RStype_hydro[2]={"HLL","HLLC"};
	public:
		// Include Riemann invariants for contour plot (only for Hydro)
                int numt;
                real tsave[ntmax];
                real Jplus[nx+2*nghosts][nbmax][ntmax], Jminus[nx+2*nghosts][nbmax][ntmax];
                real s[nx+2*nghosts][nbmax][ntmax]; // Entropy
		
		void solvertype() override {
			RStype=1;
			cout<<"Hydro Riemann solver ("<<RStype_hydro[RStype]
				<<") used! Includes Riemann invariant for contour plotting."<<endl;			
		}
};
