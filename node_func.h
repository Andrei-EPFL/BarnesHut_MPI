#include "mystruct.h"
#include <vector>

MyNode_val transfer_info(MyNode *node);
void serialize(MyNode *node, std::vector<MyNode_val> &vect, int depth);
MyNode* newNode(MyNode_val local_node_val);
void deSerialize(MyNode *&node, std::vector<MyNode_val> &vect);

void numNodeDepthKandLeaves(MyNode *root, int k, int *n_nodes);
void linkParticlesNodepRank(MyNode *root, int k, std::vector<MyParticle> &particles, int p, int *ln, unsigned int *n_nodes);
bool compareByprank(const MyParticle &a, const MyParticle &b);

int numNodesHeightDepthFlag(MyNode *root, int *n);
int numNodesHeightKandFlags(MyNode *root, int k, int *n);
int numNodesElemHeightK(MyNode *root, int k, int *n);
int numNodesElemFlagged(MyNode *root, int *n);