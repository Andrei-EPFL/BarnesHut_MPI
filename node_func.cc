#include <vector>
#include <iostream>
#include "node_func.h"

MyNode_val transfer_info(MyNode *node)
{
    MyNode_val local_node_val;
    
    local_node_val.elements = node->elements;
    local_node_val.index = node->index;
    local_node_val.depthflag = node->depthflag;
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
    
    // Read next item from file. If theere are no more items or next 
    // item is marker, then return 
    //if(!vect)
    //{
    //    std::cout<<"The vector is empty;\n"
    //} 
    // Else create node with this item and recur for children 
    if(vect.size())
    {
        MyNode_val local_node_val = vect[0];
        //std::cout<<"\n\nTrial "<<vect[0].elements<<"\n\n";
        node = newNode(vect[0]);
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


int numNodesHeightK(MyNode *root, int k)
{
    if(root == NULL) return 0; //if the tree is empty return 0
    if(k == 0)
    {
        root->depthflag=1;
        return 1; //if k = 0, then the root is the only node to return 
    }
    return numNodesHeightK(root->nwf, k-1) + numNodesHeightK(root->nef, k-1) + numNodesHeightK(root->swf, k-1) + numNodesHeightK(root->sef, k-1)+ numNodesHeightK(root->nwb, k-1) + numNodesHeightK(root->neb, k-1) + numNodesHeightK(root->swb, k-1) + numNodesHeightK(root->seb, k-1);
}

void flagParticlesToNode(MyNode *root, int k, std::vector<MyParticle> &particles)
{
    if(root != NULL)
    {   
        if(k == 0)
        {
            for(unsigned i = 0; i < particles.size(); i++)
            {
                if(particles[i].x <= root->bound_max_x && particles[i].x > root->bound_min_x && 
                   particles[i].y <= root->bound_max_y && particles[i].y > root->bound_min_y &&
                   particles[i].z <  root->bound_max_y && particles[i].z >= root->bound_min_y)
                   {particles[i].node_index = root->index;}               
            }
        }
        else
        {
            flagParticlesToNode(root->nwf, k-1, particles);
            flagParticlesToNode(root->nef, k-1, particles);
            flagParticlesToNode(root->swf, k-1, particles);
            flagParticlesToNode(root->sef, k-1, particles);
            flagParticlesToNode(root->nwb, k-1, particles);
            flagParticlesToNode(root->neb, k-1, particles);
            flagParticlesToNode(root->swb, k-1, particles);
            flagParticlesToNode(root->seb, k-1, particles);
        }
    }
}