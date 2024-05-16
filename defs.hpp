#define eps 1.0E-8
using real=float;
using namespace std;

// Global variables
const bool MAG_FIELD_ENABLED=true;
const int nvar=6; 
const int nx=5;
const int nlevs=7;
const int nbmax=64;

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
};

// Try binary tree
class AMRmesh {
	public:
		int nb;
		real u[nvar][nx+2][nbmax];
		real prim[nvar][nx+2][nbmax];
		real dx[nlevs];
		int minID[nlevs], maxID[nlevs], ActiveBlocks[nbmax];
		int lastActive=1;
		bool FlagRefine[nbmax], FlagCoarse[nbmax];
		AMRmesh* left;
		AMRmesh* right;
		AMRmesh(int val) {
			nb=val;
			left=nullptr;
			right=nullptr;
		}
};
