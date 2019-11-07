#include <vector>
#include <iostream>
#include "node_func.h"

MyNode_val transfer_info(MyNode *node)
{
    MyNode_val local_node_val;
    
    local_node_val.elements = node->elements;
    local_node_val.index = node->index;
    local_node_val.depthflag = node->depthflag;
    local_node_val.proc_rank = node->proc_rank;
    local_node_val.totalmass = node->totalmass;
    
    local_node_val.COM_x = node->COM_x;
    local_node_val.COM_y = node->COM_y;
    local_node_val.COM_z = node->COM_z;
    
    local_node_val.COM_vx = node->COM_vx;
    local_node_val.COM_vy = node->COM_vy;
    local_node_val.COM_vz = node->COM_vz;

    local_node_val.bound_min_x = node->bound_min_x;
    local_node_val.bound_min_y = node->bound_min_y;
    local_node_val.bound_min_z = node->bound_min_z;

    local_node_val.bound_max_x = node->bound_max_x;
    local_node_val.bound_max_y = node->bound_max_y;
    local_node_val.bound_max_z = node->bound_max_z;

    local_node_val.nwf = (node->nwf ? true : false);
    local_node_val.nef = (node->nef ? true : false);
    local_node_val.swf = (node->swf ? true : false);
    local_node_val.sef = (node->sef ? true : false);

    local_node_val.nwb = (node->nwb ? true : false);
    local_node_val.neb = (node->neb ? true : false);
    local_node_val.swb = (node->swb ? true : false);
    local_node_val.seb = (node->seb ? true : false);

    return local_node_val;
}


void serialize(MyNode *node, std::vector<MyNode_val> &vect, int depth) 
{ 
    if (node == NULL) 
    { 
        std::cout<<"There is no root node\n"; 
    }     
    else
    {
        if(depth != 0)
        {
            MyNode_val local_node_val;
            local_node_val = transfer_info(node);
            vect.push_back(local_node_val);

            if(node->nwf) {serialize(node->nwf, vect, depth-1);} 
            if(node->nef) {serialize(node->nef, vect, depth-1);}
            if(node->swf) {serialize(node->swf, vect, depth-1);} 
            if(node->sef) {serialize(node->sef, vect, depth-1);} 
            if(node->nwb) {serialize(node->nwb, vect, depth-1);} 
            if(node->neb) {serialize(node->neb, vect, depth-1);}
            if(node->swb) {serialize(node->swb, vect, depth-1);} 
            if(node->seb) {serialize(node->seb, vect, depth-1);} 
        }
        else if(depth == 0)
        {
            MyNode_val local_node_val;
            local_node_val = transfer_info(node);
            vect.push_back(local_node_val);
        }
    }    
}
  


MyNode* newNode(MyNode_val local_node_val)
{
    MyNode *node = new MyNode;
    
    node->elements = local_node_val.elements;
    node->index = local_node_val.index;
    node->depthflag = local_node_val.depthflag;
    node->proc_rank = local_node_val.proc_rank;
    
    node->totalmass = local_node_val.totalmass;

    node->COM_x = local_node_val.COM_x;
    node->COM_y = local_node_val.COM_y;
    node->COM_z = local_node_val.COM_z;
    
    node->COM_vx = local_node_val.COM_vx;
    node->COM_vy = local_node_val.COM_vy;
    node->COM_vz = local_node_val.COM_vz;

    node->bound_min_x = local_node_val.bound_min_x;
    node->bound_min_y = local_node_val.bound_min_y;
    node->bound_min_z = local_node_val.bound_min_z;

    node->bound_max_x = local_node_val.bound_max_x;
    node->bound_max_y = local_node_val.bound_max_y;
    node->bound_max_z = local_node_val.bound_max_z;

    node->nwf = NULL;
    node->nef = NULL;
    node->swf = NULL;
    node->sef = NULL;
    node->nwb = NULL;
    node->neb = NULL;
    node->swb = NULL;
    node->seb = NULL;
    return node;
}


void deSerialize(MyNode *&node, std::vector<MyNode_val> &vect) 
{ 
    if(!vect.empty())
    {
        MyNode_val local_node_val = *vect.begin();
        //std::cout<<"\n\nTrial "<<vect[0].elements<<"\n\n";
        node = newNode(*vect.begin());
        vect.erase(vect.begin()); 
        if(local_node_val.depthflag == 0)
        {
            if(local_node_val.nwf) {deSerialize(node->nwf, vect);} 
            if(local_node_val.nef) {deSerialize(node->nef, vect);}
            if(local_node_val.swf) {deSerialize(node->swf, vect);} 
            if(local_node_val.sef) {deSerialize(node->sef, vect);} 
            if(local_node_val.nwb) {deSerialize(node->nwb, vect);} 
            if(local_node_val.neb) {deSerialize(node->neb, vect);}
            if(local_node_val.swb) {deSerialize(node->swb, vect);} 
            if(local_node_val.seb) {deSerialize(node->seb, vect);} 
        }
    }
} 

void numNodeDepthKandLeaves(MyNode *root, int k, int *n_nodes)
{
    //This function flags the nodes at depth k and the leaves at depth < k
    //Returns by reference the number of flagged nodes.
    if(root == NULL){std::cout<<"There is no root node\n";}
    else
    {   
        if(k == 0)
        {
            *n_nodes = *n_nodes + 1;
            root->depthflag = 1;          
        }
        else
        {
            if(root->nwf){numNodeDepthKandLeaves(root->nwf, k-1, n_nodes);}
            if(root->nef){numNodeDepthKandLeaves(root->nef, k-1, n_nodes);}
            if(root->swf){numNodeDepthKandLeaves(root->swf, k-1, n_nodes);}
            if(root->sef){numNodeDepthKandLeaves(root->sef, k-1, n_nodes);}
            if(root->nwb){numNodeDepthKandLeaves(root->nwb, k-1, n_nodes);}
            if(root->neb){numNodeDepthKandLeaves(root->neb, k-1, n_nodes);}
            if(root->swb){numNodeDepthKandLeaves(root->swb, k-1, n_nodes);}
            if(root->seb){numNodeDepthKandLeaves(root->seb, k-1, n_nodes);}
            if(root->nwf == NULL && root->nef == NULL && root->swf ==NULL && root->sef == NULL && root->nwb == NULL && root->neb == NULL && root->swb ==NULL && root->seb == NULL)
            {
                *n_nodes = *n_nodes + 1;
                root->depthflag = 1; 
            }
        }
    }
}

void linkParticlesNodepRank(MyNode *root, int k, std::vector<MyParticle> &particles, int p, int *ln, unsigned int *n_nodes)
{   
    if(root == NULL){std::cout<<"There is no root node\n";}
    else
    {
        if(*ln > 0 && root->proc_rank<0)
        {
            if(k == 0)
            {
                root->proc_rank = p;
                *ln = *ln - 1;
                *n_nodes = *n_nodes + root->elements;
                for(unsigned i = 0; i < particles.size(); i++)
                {
                    if(particles[i].x <= root->bound_max_x && particles[i].x > root->bound_min_x && 
                    particles[i].y <= root->bound_max_y && particles[i].y > root->bound_min_y &&
                    particles[i].z <  root->bound_max_z && particles[i].z >= root->bound_min_z)
                        {
                            particles[i].proc_rank = p;
                            particles[i].node_index = root->index;            
                        }               
                }   
            }
            else
            {
                if(root->nwf){linkParticlesNodepRank(root->nwf, k-1, particles, p, ln, n_nodes);}
                if(root->nef){linkParticlesNodepRank(root->nef, k-1, particles, p, ln, n_nodes);}
                if(root->swf){linkParticlesNodepRank(root->swf, k-1, particles, p, ln, n_nodes);}
                if(root->sef){linkParticlesNodepRank(root->sef, k-1, particles, p, ln, n_nodes);}
                if(root->nwb){linkParticlesNodepRank(root->nwb, k-1, particles, p, ln, n_nodes);}
                if(root->neb){linkParticlesNodepRank(root->neb, k-1, particles, p, ln, n_nodes);}
                if(root->swb){linkParticlesNodepRank(root->swb, k-1, particles, p, ln, n_nodes);}
                if(root->seb){linkParticlesNodepRank(root->seb, k-1, particles, p, ln, n_nodes);}
                if(root->nwf == NULL && root->nef == NULL && root->swf ==NULL && root->sef == NULL && root->nwb == NULL && root->neb == NULL && root->swb ==NULL && root->seb == NULL)
                {
                    root->proc_rank = p;
                    *ln = *ln - 1;
                    *n_nodes = *n_nodes + root->elements;
                    for(unsigned i = 0; i < particles.size(); i++)
                    {
                        if(particles[i].x <= root->bound_max_x && particles[i].x > root->bound_min_x && 
                        particles[i].y <= root->bound_max_y && particles[i].y > root->bound_min_y &&
                        particles[i].z <  root->bound_max_z && particles[i].z >= root->bound_min_z)
                        {
                            particles[i].proc_rank = p;
                            particles[i].node_index = root->index;            
                        }               
                    }
                }
            }
        }
    }
}

bool compareByprank(const MyParticle &a, const MyParticle &b)
{
    return a.proc_rank < b.proc_rank;
}

///////////////// The following functions are useful for debugging

int numNodesHeightDepthFlag(MyNode *root, int *n)
{
    //This functions returns by reference the total number of elements from all flagged nodes and
    //returns the depth of the flagged nodes.
    if(root == NULL) return 0;
    if(root->depthflag == 1)
    {
        *n = *n+ root->elements;
        std::cout<<"Index of the node "<<root->index<<"; The no of elem of that node "<<root->elements<<"; The assigned process "<<root->proc_rank<<std::endl;
        return 1; 
    }
    return numNodesHeightDepthFlag(root->nwf, n) + numNodesHeightDepthFlag(root->nef, n) + numNodesHeightDepthFlag(root->swf, n) + numNodesHeightDepthFlag(root->sef, n)+ numNodesHeightDepthFlag(root->nwb, n) + numNodesHeightDepthFlag(root->neb, n) + numNodesHeightDepthFlag(root->swb, n) + numNodesHeightDepthFlag(root->seb, n);
}

int numNodesHeightKandFlags(MyNode *root, int k, int *n)
{
    //This function flags the nodes that are exactly at the depth k
    //Returns the number of nodes at precisely depth k
    //Returns by reference the number of elements at depth k
    if(root == NULL) return 0; 
    if(k == 0)
    {
        root->depthflag=1;
        *n = *n+ root->elements;
        return 1;
    }
    return numNodesHeightKandFlags(root->nwf, k-1, n) + numNodesHeightKandFlags(root->nef, k-1, n) + numNodesHeightKandFlags(root->swf, k-1, n) + numNodesHeightKandFlags(root->sef, k-1, n)+ numNodesHeightKandFlags(root->nwb, k-1, n) + numNodesHeightKandFlags(root->neb, k-1, n) + numNodesHeightKandFlags(root->swb, k-1, n) + numNodesHeightKandFlags(root->seb, k-1, n);
}

int numNodesElemHeightK(MyNode *root, int k, int *n)
{
    //This function returns the number of nodes precisely at depth k
    //This function returns by reference the number of elements at depth k
    if(root == NULL) return 0; 
    if(k == 0)
    {
        *n = *n+ root->elements;
        return 1;  
    }
    return numNodesElemHeightK(root->nwf, k-1, n) + numNodesElemHeightK(root->nef, k-1, n) + numNodesElemHeightK(root->swf, k-1, n) + numNodesElemHeightK(root->sef, k-1, n)+ numNodesElemHeightK(root->nwb, k-1, n) + numNodesElemHeightK(root->neb, k-1, n) + numNodesElemHeightK(root->swb, k-1, n) + numNodesElemHeightK(root->seb, k-1, n);
}

int numNodesElemFlagged(MyNode *root, int *n)
{
    //This function returns the number of elements from all flagged nodes.
    //Returns the number of flagged nodes.
    if(root == NULL) return 0;
    if(root->depthflag==1)
    {
        *n = *n+ root->elements;
        return 1;
    }
    return numNodesElemFlagged(root->nwf, n) + numNodesElemFlagged(root->nef, n) + numNodesElemFlagged(root->swf, n) + numNodesElemFlagged(root->sef, n)+ numNodesElemFlagged(root->nwb, n) + numNodesElemFlagged(root->neb, n) + numNodesElemFlagged(root->swb, n) + numNodesElemFlagged(root->seb, n);
}

