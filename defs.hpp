#define eps 1.0E-8
using real=float;
using namespace std;

// Global variables
const bool MAG_FIELD_ENABLED=false;
const int nvar=4; 
const int nx=5;
const int nlevs=7;
const int nbmax=64;
const int ntmax=1500; // Number of timesteps to be saved for plotting Riemann invariants

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
		// Include Riemann invariants for contour plot
		int numt;
		real tsave[ntmax];
		real Jplus[nx+2][nbmax][ntmax], Jminus[nx+2][nbmax][ntmax];
		real s[nx+2][nbmax][ntmax]; // Entropy
};

