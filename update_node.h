#include "mystruct.h"

void add_particle(MyNode *node, MyParticle newparticle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index);
MyNode* initialize_node(MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index);
void free_node(MyNode *node);
void update_child_node(MyNode *node, MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index);
void update_current_node(MyNode *node, MyParticle particle);