#include "update_node.h"
#include <iostream>

MyNode* initialize_node(MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index)
{
    MyNode *node = new MyNode;
    node->elements = 1;
    node->index = *index;
    *index = *index + 1;
    node->depthflag = 0;
    node->proc_rank = -1;
    node->totalmass = particle.mass;
    node->COM_x = particle.x;
    node->COM_y = particle.y;
    node->COM_z = particle.z;          

    node->COM_vx = particle.vx;
    node->COM_vy = particle.vy;           
    node->COM_vz = particle.vz;           
    

    node->bound_min_x = bound_min_x;
    node->bound_max_x = bound_max_x;
    node->bound_min_y = bound_min_y;
    node->bound_max_y = bound_max_y;
    node->bound_min_z = bound_min_z;
    node->bound_max_z = bound_max_z;


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

void free_node(MyNode *node)
{   
    if(node->nwf){free_node(node->nwf);}
    if(node->nef){free_node(node->nef);}
    if(node->swf){free_node(node->swf);}
    if(node->sef){free_node(node->sef);}
    if(node->nwb){free_node(node->nwb);}
    if(node->neb){free_node(node->neb);}
    if(node->swb){free_node(node->swb);}
    if(node->seb){free_node(node->seb);}
    if(node){delete(node);}
}

void update_child_node(MyNode *node, MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index)
{
    double x = particle.x;
    double y = particle.y;
    double z = particle.z;

    double x_center = (bound_min_x + bound_max_x)/2.;
    double y_center = (bound_min_y + bound_max_y)/2.;
    double z_center = (bound_min_z + bound_max_z)/2.;
    //Behind
    if(x<=x_center && y>y_center && z<z_center &&
    x>bound_min_x && y<=bound_max_y && z>=bound_min_z)
    {
        if(node->nwb) {add_particle(node->nwb, particle, bound_min_x, x_center, y_center, bound_max_y, bound_min_z, z_center, index);}
        else {node->nwb = initialize_node(particle, bound_min_x, x_center, y_center, bound_max_y,bound_min_z, z_center, index);}
    }
    if(x<=x_center && y<=y_center&& z<z_center
    &&x>bound_min_x && y>bound_min_y && z>=bound_min_z) 
    {
        if(node->swb) {add_particle(node->swb, particle, bound_min_x, x_center, bound_min_y, y_center, bound_min_z, z_center, index );}
        else {node->swb = initialize_node(particle, bound_min_x, x_center, bound_min_y, y_center, bound_min_z, z_center, index);}
    }
    if(x>x_center && y>y_center&& z<z_center
    && x<=bound_max_x && y<=bound_max_y && z>=bound_min_z)
    {
        if(node->neb) {add_particle(node->neb, particle, x_center, bound_max_x, y_center, bound_max_y, bound_min_z, z_center, index);}
        else {node->neb = initialize_node(particle, x_center, bound_max_x, y_center, bound_max_y, bound_min_z, z_center, index);}
    }
    if(x>x_center && y<=y_center&& z<z_center&&
    x<=bound_max_x && y>bound_min_y && z>=bound_min_z) 
    {
        if(node->seb) {add_particle(node->seb, particle, x_center, bound_max_x, bound_min_y, y_center, bound_min_z, z_center, index );}
        else {node->seb = initialize_node(particle, x_center, bound_max_x, bound_min_y, y_center, bound_min_z, z_center, index );}
    }
    //Front
    if(x<=x_center && y>y_center && z>=z_center&&
    x>bound_min_x && y<=bound_max_y && z<bound_max_z)
    {
        if(node->nwf) {add_particle(node->nwf, particle, bound_min_x, x_center, y_center, bound_max_y, z_center, bound_max_z, index);}
        else {node->nwf = initialize_node(particle, bound_min_x, x_center, y_center, bound_max_y, z_center, bound_max_z, index);}
    }
    if(x<=x_center && y<=y_center&& z>=z_center
    &&x>bound_min_x && y>bound_min_y && z<bound_max_z) 
    {
        if(node->swf) {add_particle(node->swf, particle, bound_min_x, x_center, bound_min_y, y_center, z_center, bound_max_z, index );}
        else {node->swf = initialize_node(particle, bound_min_x, x_center, bound_min_y, y_center, z_center, bound_max_z, index);}
    }
    if(x>x_center && y>y_center&& z>=z_center
    &&x<=bound_max_x && y<=bound_max_y && z<bound_max_z)
    {
        if(node->nef) {add_particle(node->nef, particle, x_center, bound_max_x, y_center, bound_max_y, z_center, bound_max_z, index);}
        else {node->nef = initialize_node(particle, x_center, bound_max_x, y_center, bound_max_y, z_center, bound_max_z, index);}
    }
    if(x>x_center && y<=y_center&& z>=z_center
    &&x<=bound_max_x && y>bound_min_y && z<bound_max_z) 
    {
        if(node->sef) {add_particle(node->sef, particle, x_center, bound_max_x, bound_min_y, y_center, z_center, bound_max_z, index );}
        else {node->sef = initialize_node(particle, x_center, bound_max_x, bound_min_y, y_center, z_center, bound_max_z, index );}
    }
}

void update_current_node(MyNode *node, MyParticle particle)
{
    node->elements++;
    node->COM_x = (node->totalmass * node->COM_x + particle.mass * particle.x)/(node->totalmass+particle.mass);
    node->COM_y = (node->totalmass * node->COM_y + particle.mass * particle.y)/(node->totalmass+particle.mass);
    node->COM_z = (node->totalmass * node->COM_z + particle.mass * particle.z)/(node->totalmass+particle.mass);
    
    node->COM_vx = (node->totalmass * node->COM_vx + particle.mass * particle.vx)/(node->totalmass+particle.mass);
    node->COM_vy = (node->totalmass * node->COM_vy + particle.mass * particle.vy)/(node->totalmass+particle.mass);
    node->COM_vz = (node->totalmass * node->COM_vz + particle.mass * particle.vz)/(node->totalmass+particle.mass);
    
    node->totalmass+=particle.mass;
}

void add_particle(MyNode *node, MyParticle newparticle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y, double bound_min_z, double bound_max_z, int *index)
{
    
    if(!node)
    {       
        if(newparticle.outside == false)
        {
            std::cout << "This is the !node branch of the if clauses\n";      
        }
        else
        {
            std::cout << "The particle is placed outside the boundaries(Node\n";
        }
    }
    else
    {
        //std::cout<<"This is the other branch of the if clause\n";
        if(newparticle.outside == false)
        {
            if(node->elements == 1)
            {
                MyParticle oldparticle;
                //std::cout<<"This node has one element\n";
                oldparticle.x = node->COM_x;
                oldparticle.y = node->COM_y;
                oldparticle.z = node->COM_z;
                oldparticle.vx = node->COM_vx;
                oldparticle.vy = node->COM_vy;
                oldparticle.vz = node->COM_vz;
                oldparticle.mass = node->totalmass;
                
                update_child_node(node, oldparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, index);
            }
            update_current_node(node, newparticle);
            update_child_node(node, newparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, index);
        }
        else
        {
            //std::cout << "The particle is placed outside the boundaries\n";
            //std::cout <<newparticle.outside<<std::endl;
        }
    }
}

void add_particle_locally(MyNode *node, MyParticle newparticle, int prank, int *index)
{
    
    if(!node)
    {       
        std::cout << "ERROR: This is odd!! The root should exists because it was broadcasted to every process\n";    
    }
    else
    {
        if(newparticle.outside == false)
        {

            if(newparticle.node_index != node->index)
            {
                if(node->nwf) {add_particle_locally(node->nwf, newparticle, prank, index);}
                if(node->nef) {add_particle_locally(node->nef, newparticle, prank, index);}
                if(node->swf) {add_particle_locally(node->swf, newparticle, prank, index);}
                if(node->sef) {add_particle_locally(node->sef, newparticle, prank, index);}
                if(node->nwb) {add_particle_locally(node->nwb, newparticle, prank, index);}
                if(node->neb) {add_particle_locally(node->neb, newparticle, prank, index);}
                if(node->swb) {add_particle_locally(node->swb, newparticle, prank, index);}
                if(node->seb) {add_particle_locally(node->seb, newparticle, prank, index);}
            }
            else
            {
                if(node->depthflag==1 && node->elements>1)
                {
                    update_child_node(node, newparticle, node->bound_min_x, node->bound_max_x, node->bound_min_y, node->bound_max_y, node->bound_min_z, node->bound_max_z, index);
                }
                else 
                
                {
                    if(node->depthflag==0){std::cout<<"ERROR: This should not have occured\n";}
                    if(node->elements<=0){std::cout<<"ERROR: This should not have occured\n";}
                    if(node->elements==1 && (node->COM_x != newparticle.x || node->COM_y != newparticle.y || node->COM_z != newparticle.z )){std::cout<<"ERROR: This should not have occured\n";}
                }
            }

        }
        else
        {
            //std::cout << "The particle is placed outside the boundaries\n";
            //std::cout <<newparticle.outside<<std::endl;
        }
    }
}
