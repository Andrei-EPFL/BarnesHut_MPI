#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "update_node.h"
#include "dynamics.h"

//#include <mpi.h>

int main()
{
    //MPI_Init();
    //int prank, psize;
    //MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    //MPI_Comm_size(MPI_COMM_WORLD, &psize);

    //Declartion of variables
    std::ifstream infile;
    std::ofstream ofile;
    MyNode *root = NULL;
    
    int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 500;
    double bound_min_y = 0.;
    double bound_max_y = 500;

    std::vector<MyParticle> particles;
    MyParticle tmpparticle;

    infile.open("disk.txt", std::ios::in);
    //Initialisation of the root node
    infile>>tmpparticle.x>>tmpparticle.y>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.mass;
    root = initialize_node(tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y);
    particles.push_back(tmpparticle);
    while(infile>>tmpparticle.x)
    {
        infile>>tmpparticle.y>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.mass;
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
    float dt = 0.01;
    
    //Computation of new positions
    ofile.open("diskout.txt", std::ios::out);
    for(int step = 0; step<500; step++)
    {
        fx = fy = 0.;
      
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
            //if(particles[i].x < bound_min_x){particles[i].x = bound_max_x - (bound_min_x - particles[i].x);}
            //if(particles[i].y < bound_min_y){particles[i].y = bound_max_y - (bound_min_y - particles[i].y);}
            //if(particles[i].x > bound_max_x){particles[i].x = bound_min_x + (particles[i].x - bound_max_x);}
            //if(particles[i].y > bound_max_y){particles[i].y = bound_min_y + (particles[i].y - bound_max_y);}
            if(particles[i].x < bound_min_x){particles.erase(particles.begin()+i); n--;i--;}
            if(particles[i].y < bound_min_y){particles.erase(particles.begin()+i); n--;i--;}
            if(particles[i].x > bound_max_x){particles.erase(particles.begin()+i); n--;i--;}
            if(particles[i].y > bound_max_y){particles.erase(particles.begin()+i); n--;i--;}

 
            fx = fy = 0.;
            ofile<<particles[i].x<<" "<<particles[i].y<<std::endl;
        }
        ofile<<step<<std::endl;
       
        free_node(root);
        root = NULL;
        root = initialize_node(particles[0], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        for(int p = 1; p < n; p++)
        {        
            add_particle(root, particles[p], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        }
    }
    ofile.close();

    std::cout<<"End of program"<< std::endl;
    //MPI_Finalize();
    return 0;
}



