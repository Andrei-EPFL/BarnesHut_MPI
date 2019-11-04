#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <numeric>
#include "update_node.h"
#include "dynamics.h"
#include <omp.h>

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

int main()
{
    auto nthreads = omp_get_max_threads();
    //Declartion of variables
    std::ifstream infile;
    std::ofstream ofile;
    MyNode *root = NULL;
    
    int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 500;
    double bound_min_y = 0.;
    double bound_max_y = 500;
    double bound_min_z = 0.;
    double bound_max_z = 500;

    std::vector<MyParticle> particles;
    MyParticle tmpparticle;
    
    std::vector<double> times(1000);
    std::vector<double> times2(1000);

    auto t0 = clk::now();
    infile.open("./input/disk.txt", std::ios::in);
    //Initialisation of the root node
    infile>>tmpparticle.x>>tmpparticle.y>>tmpparticle.z>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.vz>>tmpparticle.mass;
    tmpparticle.outside = false;
    root = initialize_node(tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
    particles.push_back(tmpparticle);
    while(infile>>tmpparticle.x)
    {
        infile>>tmpparticle.y>>tmpparticle.z>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.vz>>tmpparticle.mass;
        tmpparticle.outside = false;
        particles.push_back(tmpparticle);
        add_particle(root, tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
    }

    if(root->nwf){std::cout<<root->nwf->elements<<std::endl;}
    if(root->nef){std::cout<<root->nef->elements<<std::endl;}
    if(root->swf){std::cout<<root->swf->elements<<std::endl;}
    if(root->sef){std::cout<<root->sef->elements<<std::endl;}
    if(root->nwb){std::cout<<root->nwb->elements<<std::endl;}
    if(root->neb){std::cout<<root->neb->elements<<std::endl;}
    if(root->swb){std::cout<<root->swb->elements<<std::endl;}
    if(root->seb){std::cout<<root->seb->elements<<std::endl;}

    infile.close();
    auto t1 = clk::now();

    n = root->elements;
    std::cout<<"The first creation of the tree took "<<second(t1 - t0).count() << " seconds"<<std::endl;
    std::cout<<"The root node has "<<root->elements << " elements"<<std::endl;
    std::cout<<"The particles vector has " << particles.size() << " particles" << std::endl;
    
    //Declaration of variables for the actual computation
    double fx = 0., fy = 0., fz = 0;
    double ax = 0., ay = 0., az = 0;
    float dt = 0.1;
    int Nsteps = 1000;

    //Computation of new positions
    ofile.open("./output/diskout.txt", std::ios::out);
    for(int step = 0; step<Nsteps; step++)
    {
        fx = fy = fz = 0.;
        auto t2 = clk::now();
#pragma omp parallel for private(fx, fy, fz, ax, ay, az) 
        for(int i = 0; i < n; i++)
        {
            if(particles[i].outside == false)
            {
                compute_force(root, particles[i], &fx, &fy, &fz);
                ax = fx/particles[i].mass;
                ay = fy/particles[i].mass;
                az = fz/particles[i].mass;
                
                particles[i].vx += ax * dt;
                particles[i].vy += ay * dt;
                particles[i].vz += az * dt;
                
                particles[i].x += particles[i].vx * dt;
                particles[i].y += particles[i].vy * dt;
                particles[i].z += particles[i].vz * dt;

                //This is allows a particle to go out of the boundary, by returning it on the other side of the box.
                //if(particles[i].x < bound_min_x){particles[i].x = bound_max_x - (bound_min_x - particles[i].x);}
                //if(particles[i].y < bound_min_y){particles[i].y = bound_max_y - (bound_min_y - particles[i].y);}
                //if(particles[i].x > bound_max_x){particles[i].x = bound_min_x + (particles[i].x - bound_max_x);}
                //if(particles[i].y > bound_max_y){particles[i].y = bound_min_y + (particles[i].y - bound_max_y);}
                if(particles[i].x < bound_min_x || particles[i].y < bound_min_y || particles[i].x > bound_max_x || particles[i].y > bound_max_y || particles[i].z > bound_max_z || particles[i].z < bound_min_z)
                {
                    particles[i].x = bound_min_x - 100; 
                    particles[i].y = bound_min_y - 100;
                    particles[i].z = bound_min_z - 100;
                    particles[i].vx = particles[i].vy = particles[i].vz = particles[i].mass = 0;
                    particles[i].outside = true;
                }
            }
            fx = fy = fz = 0.;
        }
        auto t3 = clk::now();
        times[step] = second (t3-t2).count();
        free_node(root);
        root = NULL;
        
        int p1 = 0;
        while(particles[p1].outside == true){p1++;}

        root = initialize_node(particles[p1], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        ofile<<particles[p1].x<<" "<<particles[p1].y<<" "<<particles[p1].z<<std::endl;
            
        for(int p = p1+1; p < n; p++)
        { 
            ofile<<particles[p].x<<" "<<particles[p].y<<" "<<particles[p].z<<std::endl;       
            add_particle(root, particles[p], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        }
        ofile<<step<<std::endl;
        auto t4 = clk::now();
        times2[step] = second (t4-t3).count();
    }
    second elapsed = clk::now() - t0;
    ofile.close();

    std::cout<<"There were "<< Nsteps <<" steps, each of "<< dt <<" Myr"<<std::endl;
    std::cout<<"The remaining number of particles in the particles vector is= "<<n <<std::endl;
    std::cout<<"The number of particles in the tree is= "<<root->elements <<std::endl;
    std::cout<<"The large loop with steps takes "<<elapsed.count() << " seconds using "<<nthreads<<" threads"<<std::endl;
    std::cout<<"The average time to compute forces for all particles in one step is "<< std::accumulate(times.begin(), times.end(), 0.0)/times.size() << " seconds"<<std::endl;
    std::cout<<"The average time to empty and refill the root node is "<< std::accumulate(times2.begin(), times2.end(), 0.0)/times2.size() << " seconds"<<std::endl;
    std::cout<<"End of program"<< std::endl;
    //MPI_Finalize();
    return 0;
}



