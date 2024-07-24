#include "defs.hpp"

void meshblock::setParams() {
	ntot=nx+2*nghosts-1;
	nxmin=nghosts; nxmax=ntot-nghosts;
	nxminb=nxmin-1; nxmaxb=nxmax+1;
	nx2=ntot/2;
}
