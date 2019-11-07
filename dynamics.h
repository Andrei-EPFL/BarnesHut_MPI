#include<vector>
#include "mystruct.h"

int compute_force(MyNode *node, MyParticle particle, double *fx, double *fy, double *fz);
int compute_force_partially(MyNode *node, MyParticle particle, double *fx, double *fy, double *fz, std::vector<std::vector<MyParticle>> &mat_particles_send, int puf);

