#include "update_node.h"
#include <iostream>

MyNode* initialize_node(MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y)
{
    MyNode *node = new MyNode;
    node->elements = 1;
    node->totalmass = particle.mass;
    node->COM_x = particle.x;
    node->COM_y = particle.y;           

    node->COM_vx = particle.vx;
    node->COM_vy = particle.vy;           
    
    node->bound_min_x = bound_min_x;
    node->bound_max_x = bound_max_x;
    node->bound_min_y = bound_min_y;
    node->bound_max_y = bound_max_y;

    node->nw = NULL;
    node->ne = NULL;
    node->sw = NULL;
    node->se = NULL;
    return node;
}

void free_node(MyNode *node)
{   
    free_node(node->nw);
    free_node(node->ne);
    free_node(node->sw);
    free_node(node->se);
    delete(node);
}

void update_child_node(MyNode *node, MyParticle particle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y)
{
    double x = particle.x;
    double y = particle.y;
    double x_center = (bound_min_x + bound_max_x)/2.;
    double y_center = (bound_min_y + bound_max_y)/2.;

    if(x<=x_center && y>y_center)
    {
        if(node->nw) {add_particle(node->nw, particle, bound_min_x, x_center, y_center, bound_max_y);}
        else {node->nw = initialize_node(particle, bound_min_x, x_center, y_center, bound_max_y);}
    }
    if(x<=x_center && y<=y_center) 
    {
        if(node->sw) {add_particle(node->sw, particle, bound_min_x, x_center, bound_min_y, y_center );}
        else {node->sw = initialize_node(particle, bound_min_x, x_center, bound_min_y, y_center);}
    }
    if(x>x_center && y>y_center)
    {
        if(node->ne) {add_particle(node->ne, particle, x_center, bound_max_x, y_center, bound_max_y);}
        else {node->ne = initialize_node(particle, x_center, bound_max_x, y_center, bound_max_y);}
    }
    if(x>x_center && y<=y_center) 
    {
        if(node->se) {add_particle(node->se, particle, x_center, bound_max_x, bound_min_y, y_center );}
        else {node->se = initialize_node(particle, x_center, bound_max_x, bound_min_y, y_center );}
    }
}

void update_current_node(MyNode *node, MyParticle particle)
{
    node->elements++;
    node->COM_x = (node->totalmass * node->COM_x + particle.mass * particle.x)/(node->totalmass+particle.mass);
    node->COM_y = (node->totalmass * node->COM_y + particle.mass * particle.y)/(node->totalmass+particle.mass);
    node->COM_vx = (node->totalmass * node->COM_vx + particle.mass * particle.vx)/(node->totalmass+particle.mass);
    node->COM_vy = (node->totalmass * node->COM_vy + particle.mass * particle.vy)/(node->totalmass+particle.mass);
    node->totalmass+=particle.mass;
}

void add_particle(MyNode *node, MyParticle newparticle, double bound_min_x, double bound_max_x, double bound_min_y, double bound_max_y)
{
    double x = newparticle.x;
    double y = newparticle.y;
    
    if(!node)
    {       
        if(y>=bound_min_y && y<=bound_max_y && x>=bound_min_x && x<=bound_max_x)
        {
            std::cout << "This is the !node branch of the if clauses\n";
            
        }
        else
        {
            std::cout << "The particle is placed outside the boundaries\n";
        }
    }
    else
    {
        //std::cout<<"This is the other branch of the if clause\n";
        if(y>=bound_min_y && y<=bound_max_y && x>=bound_min_x && x<=bound_max_x)
        {
            if(node->elements == 1)
            {
                MyParticle oldparticle;
                //std::cout<<"This node has one element\n";
                oldparticle.x = node->COM_x;
                oldparticle.y = node->COM_y;
                oldparticle.vx = node->COM_vx;
                oldparticle.vy = node->COM_vy;
                oldparticle.mass = node->totalmass;
                
                update_child_node(node, oldparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y);
            }
            update_current_node(node, newparticle);
            update_child_node(node, newparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        }
        else
        {
            std::cout << "The particle is placed outside the boundaries\n";
        }
    }
}
