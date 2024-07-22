int init_mesh(meshblock* dom);
int getlevel(meshblock* dom, const int nb);
real getBlockPosition(meshblock* dom, const int nb);
int get_nb(meshblock* dom, const int bID);
int getLeft(meshblock* dom, const int nb);
int getRight(meshblock* dom, const int nb);
void refineBlock(meshblock* dom, const int dadNb);
void coarseBlock(meshblock* dom, int son1nb);
void update_mesh(meshblock* dom);
