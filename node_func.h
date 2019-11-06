#include "mystruct.h"
#include <vector>

MyNode_val transfer_info(MyNode *node);
void serialize(MyNode *node, std::vector<MyNode_val> &vect, int depth);
MyNode* newNode(MyNode_val local_node_val);
void deSerialize(MyNode *&node, std::vector<MyNode_val> &vect);
int numNodesHeightK(MyNode *root, int k);
void flagParticlesToNode(MyNode *root, int k, std::vector<MyParticle> &particles);