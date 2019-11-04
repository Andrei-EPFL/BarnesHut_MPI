#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>
#include "update_node.h"
#include "dynamics.h"

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

#include <mpi.h>

int main()
{
    MPI_Init(NULL, NULL);
    int prank, psize;
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

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
    
    int blk_length[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint address[9];
    MPI_Get_address(&tmpparticle, &address[0]);
    MPI_Get_address(&tmpparticle.x, &address[1]);
    MPI_Get_address(&tmpparticle.y, &address[2]);
    MPI_Get_address(&tmpparticle.z, &address[3]);
    MPI_Get_address(&tmpparticle.vx, &address[4]);
    MPI_Get_address(&tmpparticle.vy, &address[5]);
    MPI_Get_address(&tmpparticle.vz, &address[6]);
    MPI_Get_address(&tmpparticle.mass, &address[7]);
    MPI_Get_address(&tmpparticle.outside, &address[8]);

    MPI_Aint displs[8];
    displs[0] = MPI_Aint_diff(address[1], address[0]);
    displs[1] = MPI_Aint_diff(address[2], address[0]);
    displs[2] = MPI_Aint_diff(address[3], address[0]);
    displs[3] = MPI_Aint_diff(address[4], address[0]);
    displs[4] = MPI_Aint_diff(address[5], address[0]);
    displs[5] = MPI_Aint_diff(address[6], address[0]);
    displs[6] = MPI_Aint_diff(address[7], address[0]);
    displs[7] = MPI_Aint_diff(address[8], address[0]);
    
    MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LOGICAL};
    MPI_Datatype MyParticle_mpi_t;
    MPI_Type_create_struct(8, blk_length, displs, types, &MyParticle_mpi_t);
    MPI_Type_commit(&MyParticle_mpi_t);

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
    std::vector<MyParticle> particles_out(n);

    std::cout<<"The first creation of the tree took "<<second(t1 - t0).count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The root node has "<<root->elements << " elements"<<std::endl;
    std::cout<<"The particles vector has " << particles.size() << " particles" << std::endl;
    
    //Declaration of variables for the actual computation
    double fx = 0., fy = 0., fz = 0;
    double ax = 0., ay = 0., az = 0;
    float dt = 0.1;

    auto ln = n/psize + (prank < n % psize ? 1 : 0);
    auto i_start = prank * ln + (prank < n % psize ? 0 : n % psize);
    auto i_end = i_start + ln;
    std::array<int, 10> recvcounts;
    std::array<int, 10> displs_data;
    for(int c = 0; c<psize; c++)
    {
        recvcounts[c] = n/psize + (c < n % psize ? 1 : 0);
        displs[c] = c * ln + (c < n % psize ? 0 : n % psize);
    }
    
    //Computation of new positions
    if(prank==0){ofile.open("./output/diskout.txt", std::ios::out);}
    for(int step = 0; step<10; step++)
    {
        fx = fy = fz = 0.;
      
        for(int i = i_start; i < i_end; i++)
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
            MPI_Gatherv(&particles[i_start], ln, MyParticle_mpi_t, &particles_out, &recvcounts[0], &displs_data[0], MyParticle_mpi_t, 0, MPI_COMM_WORLD );
            if(prank==0){ofile<<particles[i].x<<" "<<particles[i].y<<" "<<particles[i].z<<std::endl;}
        }
        if(prank==0){ofile<<step<<std::endl;}
        free_node(root);
        root = NULL;
        root = initialize_node(particles[0], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        for(int p = 1; p < n; p++)
        {        
            add_particle(root, particles[p], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        }
    }
    if(prank==0){ofile.close();}
    second elapsed = clk::now() - t0;
    std::cout<<"The remaining number of particles in the particles vector is= "<<n <<std::endl;
    std::cout<<"The number of particles in the tree is= "<<root->elements <<std::endl;
    std::cout<<"The large loop with steps takes "<<elapsed.count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"End of program"<< std::endl;
    MPI_Finalize();
    return 0;
}



