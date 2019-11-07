#include <iostream>
#include <cmath>
#include <vector>
#include "dynamics.h"
#define THETA_DEF 0.01

int compute_force(MyNode *node, MyParticle particle, double *fx, double *fy, double *fz)
{
    const float theta = THETA_DEF;
    //const double G = 6.672e-11; //m, kg, s
    const double G = 4.49e-11;//kpc, Msun, MegaYear
    const float epsilon = 0.1; //kpc
    if(!node) {std::cout<<"There is no node, so Bye Bye\n"; return 0;}
    *fx = *fy = *fz = 0;
    // Compute the distance between the COM of the node and the particle
    double distance = std::sqrt(std::pow(node->COM_x-particle.x, 2) + 
                              std::pow(node->COM_y-particle.y, 2) + 
                              std::pow(node->COM_z-particle.z, 2));
    // If this distance is zero, it means that the node contains only that particle. This would be autointeraction.
    if(distance==0) {/*std::cout<<"The particle is the node itself so the distance = 0\n";*/ return 1;}
    
    //Compute the average size of the quadrant (the node bounding area).
    double quadrant_avg_size = ((node->bound_max_x-node->bound_min_x) + 
		                (node->bound_max_y-node->bound_min_y)+
                        (node->bound_max_z-node->bound_min_z)) / 3.;

    distance = std::sqrt(std::pow(distance,2) + std::pow(epsilon,2));
    /*If the ratio below is smaller than theta, the node is far away such that it can be treated
    as a particle with totalmass at position COM_x and COM_y. Thus the computation of the force can be done*/
    if(node->elements == 1)
    {
        double distance3 = std::pow(distance, 3);
        *fx = G * particle.mass * node->totalmass * (node->COM_x - particle.x)/distance3;
        *fy = G * particle.mass * node->totalmass * (node->COM_y - particle.y)/distance3;
        *fz = G * particle.mass * node->totalmass * (node->COM_z - particle.z)/distance3;
    }
    else if(node->elements > 1)
    {
        if(quadrant_avg_size/distance < theta)
        {
            double distance3 = std::pow(distance, 3);
            *fx = G * particle.mass * node->totalmass * (node->COM_x - particle.x)/distance3;
            *fy = G * particle.mass * node->totalmass * (node->COM_y - particle.y)/distance3;
            *fz = G * particle.mass * node->totalmass * (node->COM_z - particle.z)/distance3;
        }
        else
        {
            /*Else, we need to examinate the children of this node to compute the force.*/
            double child_fx = 0; double child_fy = 0; double child_fz = 0;
            if (node->nwf)
            {
                compute_force(node->nwf, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->nef)
            {
                compute_force(node->nef, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->swf)
            {
                compute_force(node->swf, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->sef)
            {
                compute_force(node->sef, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->nwb)
            {
                compute_force(node->nwb, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->neb)
            {
                compute_force(node->neb, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->swb)
            {
                compute_force(node->swb, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->seb)
            {
                compute_force(node->seb, particle, &child_fx, &child_fy, &child_fz );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
        }
    }
    else {std::cout<<"This node has no elements\n";}
    return 1;
}


int compute_force_partially(MyNode *node, MyParticle particle, double *fx, double *fy, double *fz, std::vector<std::vector<MyParticle>> &mat_particles_send)
{
    const float theta = THETA_DEF;
    //const double G = 6.672e-11; //m, kg, s
    const double G = 4.49e-11;//kpc, Msun, MegaYear
    const float epsilon = 0.1; //kpc
    *fx = *fy = *fz = 0;
    
    if(!node) {std::cout<<"There is no node, so Bye Bye\n"; return 0;}
    
    // Compute the distance between the COM of the node and the particle
    double distance = std::sqrt(std::pow(node->COM_x-particle.x, 2) + 
                              std::pow(node->COM_y-particle.y, 2) + 
                              std::pow(node->COM_z-particle.z, 2));
    
    // If this distance is zero, it means that the node contains only that particle. This would be autointeraction.
    if(distance==0) {/*std::cout<<"The particle is the node itself so the distance = 0\n";*/ return 1;}
    
    //Compute the average size of the quadrant (the node bounding area).
    double quadrant_avg_size = ((node->bound_max_x-node->bound_min_x) + 
		                (node->bound_max_y-node->bound_min_y)+
                        (node->bound_max_z-node->bound_min_z)) / 3.;

    if(node->elements == 1)
    {
        distance = std::sqrt(std::pow(distance,2) + std::pow(epsilon,2));
        double distance3 = std::pow(distance, 3);
        *fx = G * particle.mass * node->totalmass * (node->COM_x - particle.x)/distance3;
        *fy = G * particle.mass * node->totalmass * (node->COM_y - particle.y)/distance3;
        *fz = G * particle.mass * node->totalmass * (node->COM_z - particle.z)/distance3;
    }
    else if(node->elements > 1)
    {
        /*If the ratio below is smaller than theta, the node is far away such that it can be treated
        as a particle with totalmass at position COM_x and COM_y. Thus the computation of the force can be done*/
        if(quadrant_avg_size/distance < theta)
        {   
            distance = std::sqrt(std::pow(distance,2) + std::pow(epsilon,2));  
            double distance3 = std::pow(distance, 3);
            *fx = G * particle.mass * node->totalmass * (node->COM_x - particle.x)/distance3;
            *fy = G * particle.mass * node->totalmass * (node->COM_y - particle.y)/distance3;
            *fz = G * particle.mass * node->totalmass * (node->COM_z - particle.z)/distance3;
        }
        else
        {
            /*Else, we need to examinate the children of this node to compute the force.*/
            double child_fx = 0; double child_fy = 0; double child_fz = 0;
            if (node->nwf)
            {
                compute_force_partially(node->nwf, particle, &child_fx, &child_fy, &child_fz, mat_particles_send);
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->nef)
            {
                compute_force_partially(node->nef, particle, &child_fx, &child_fy, &child_fz, mat_particles_send);
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->swf)
            {
                compute_force_partially(node->swf, particle, &child_fx, &child_fy, &child_fz, mat_particles_send);
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->sef)
            {
                compute_force_partially(node->sef, particle, &child_fx, &child_fy, &child_fz, mat_particles_send);
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->nwb)
            {
                compute_force_partially(node->nwb, particle, &child_fx, &child_fy, &child_fz, mat_particles_send );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->neb)
            {
                compute_force_partially(node->neb, particle, &child_fx, &child_fy, &child_fz, mat_particles_send);
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->swb)
            {
                compute_force_partially(node->swb, particle, &child_fx, &child_fy, &child_fz, mat_particles_send );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            child_fx = child_fy = child_fz = 0;
            if (node->seb)
            {
                compute_force_partially(node->seb, particle, &child_fx, &child_fy, &child_fz, mat_particles_send );
                *fx += child_fx;
                *fy += child_fy;
                *fz += child_fz;
            }
            if(node->depthflag == 1)
            {
                if(node->nwf != NULL || node->nef != NULL || node->swf != NULL || node->sef != NULL || node->nwb != NULL || node->neb != NULL || node->swb != NULL || node->seb != NULL)
                {
                    std::cout<<"BRANZASTRICATA\n";
                }
                else
                {
                    std::cout<<"CASCAVALAFUMAT\n";
                }
            }
            
            if(node->nwf == NULL && node->nef == NULL && node->swf == NULL && node->sef == NULL && node->nwb == NULL && node->neb == NULL && node->swb == NULL && node->seb == NULL)
            {
                std::cout<<"The flag of the node is "<< node->depthflag <<"\n";
            } 
        }
    }
    else {std::cout<<"This node has no elements\n";}
    return 2;

    //std::cout<<"pufPUFPUFPUF: "<<puf<<" "<<node->depthflag<<" "<<node->index<<std::endl;
    //if(node->depthflag == 1) { std::cout<<"ma caca pe teot ce se poate\n";} //return 5;}
    
    //if(node->depthflag == 1) {return 3;}
    //if(node->index == 95) {return 3;}
    
     /*if(node->nwf != NULL || node->nef != NULL || node->swf != NULL || node->sef != NULL || node->nwb != NULL || node->neb != NULL || node->swb != NULL || node->seb != NULL)
            {
                std::cout<<"branzaSTRICATA\n";
            }
            std::cout<<"CACATPEBATCACATPEBAT \n";
        */


            /*if(node->depthflag == 1)
            {  
                std::cout<<"rank of process of the node " <<node->proc_rank<<std::endl; 
                    //std::cout<<"this has flag1\n";
             //       mat_particles_send[node->proc_rank].push_back(particle);
            }
            if(node->nwf == NULL && node->nef == NULL && node->swf == NULL && node->sef == NULL && node->nwb == NULL && node->neb == NULL && node->swb == NULL && node->seb == NULL)
            {
                std::cout<<"branza\n";
                return 0;
            }*/
}
