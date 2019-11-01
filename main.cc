#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "update_node.h"
#include "dynamics.h"

int main()
{
    //Declartion of variables
    std::ifstream infile;
    std::ofstream ofile;
    MyNode *root = NULL;
    
    int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 2500;
    double bound_min_y = 0.;
    double bound_max_y = 2500;

    std::vector<MyParticle> particles;
    MyParticle tmpparticle;

    infile.open("sd.dat", std::ios::in);
    float z, vz;
    //Initialisation of the root node
    infile>>tmpparticle.x>>tmpparticle.y>>z>>tmpparticle.vx>>tmpparticle.vy>>vz>>tmpparticle.mass;
    root = initialize_node(tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y);
    particles.push_back(tmpparticle);
    while(infile>>tmpparticle.x)
    {
        infile>>tmpparticle.y>>z>>tmpparticle.vx>>tmpparticle.vy>>vz>>tmpparticle.mass;
        particles.push_back(tmpparticle);
        add_particle(root, tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y);
    }
    infile.close();

    n = root->elements;
    std::cout<<"The root node has "<<root->elements << " elements"<<std::endl;
    std::cout<<"The particles vector has " << particles.size() << " particles" << std::endl;
    
    //Declaration of variables for the actual computation
    double fx = 0., fy = 0.;
    double ax = 0., ay = 0.;
    float dt = 1.;
    
    //Computation of new positions
    for(int step = 0; step<24; step++)
    {
        fx = fy = 0;
        if(step%10==0)
        {
            std::string filename="output"+std::to_string(step)+".dat";
            ofile.open(filename, std::ios::out);
        }
        for(int i = 0; i < n; i++)
        {
            compute_force(root, particles[i], &fx, &fy);
            ax = fx/particles[i].mass;
            ay = fy/particles[i].mass;
            particles[i].vx += ax * dt;
            particles[i].vy += ay * dt;
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
            
            //This is allows a particle to go out of the boundary, by returning it on the other side of the box.
            if(particles[i].x < bound_min_x){particles[i].x = bound_max_x - (bound_min_x - particles[i].x);}
            if(particles[i].y < bound_min_y){particles[i].y = bound_max_y - (bound_min_y - particles[i].y);}
            if(particles[i].x > bound_max_x){particles[i].x = bound_min_x + (particles[i].x - bound_max_x);}
            if(particles[i].y > bound_max_y){particles[i].y = bound_min_y + (particles[i].y - bound_max_y);}

            fx = fy = 0;
            if(step%10==0){ofile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].vx<<" "<<particles[i].vy<<" "<<std::endl;}
        }
        if(step%10==0){ofile.close();}

        free_node(root);
        root = NULL;
        
        root = initialize_node(particles[0], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        for(int p = 1; p < n; p++)
        {        
            add_particle(root, particles[p], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        }
    }

    std::cout<<"End of program"<< std::endl;
    return 0;
}



