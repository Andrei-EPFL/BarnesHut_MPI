#include <iostream>
#include <cmath>
#include "dynamics.h"

int compute_force(MyNode *node, MyParticle particle, double *fx, double *fy)
{
    const float theta = 0.5;
    //const double G = 6.672e-11;
    const double G = 1;
    if(!node) {std::cout<<"There is no node, so Bye Bye\n"; return 0;}
    *fx = *fy = 0;
    // Compute the distance between the COM of the node and the particle
    double distance = std::sqrt(std::pow(node->COM_x-particle.x, 2) + 
                              std::pow(node->COM_y-particle.y, 2));
    // If this distance is zero, it means that the node contains only that particle. This would be autointeraction.
    if(distance==0) {std::cout<<"The particle is the node itself so the distance = 0\n"; return 1;}
    
    //Compute the average size of the quadrant (the node bounding area).
    double quadrant_avg_size = ((node->bound_max_x-node->bound_min_x) + 
		                (node->bound_max_y-node->bound_min_y)) / 2.;

    /*If the ratio below is smaller than theta, the node is far away such that it can be treated
    as a particle with totalmass at position COM_x and COM_y. Thus the computation of the force can be done*/
    std::cout<<quadrant_avg_size/distance<<std::endl;
    if(quadrant_avg_size/distance < theta)
    {
        double distance3 = std::pow(distance, 3);
        *fx = G * particle.mass * node->totalmass * (node->COM_x - particle.x)/distance3;
        *fy = G * particle.mass * node->totalmass * (node->COM_y - particle.y)/distance3;
    }
    else
    {
        /*Else, we need to examinate the children of this node to compute the force.*/
        double child_fx = 0; double child_fy = 0;
        if (node->nw)
        {
            compute_force(node->nw, particle, &child_fx, &child_fy);
            *fx += child_fx;
            *fy += child_fy;
        }
       child_fx = 0; child_fy = 0;
        if (node->ne)
        {
            compute_force(node->ne, particle, &child_fx, &child_fy);
            *fx += child_fx;
            *fy += child_fy;
        }
        child_fx = 0; child_fy = 0;
        if (node->sw)
        {
            compute_force(node->sw, particle, &child_fx, &child_fy);
            *fx += child_fx;
            *fy += child_fy;
        }
        child_fx = 0; child_fy = 0;
        if (node->se)
        {
            compute_force(node->se, particle, &child_fx, &child_fy);
            *fx += child_fx;
            *fy += child_fy;
        }
    }
    return 1;
}