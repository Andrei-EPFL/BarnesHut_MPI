#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>
#include "update_node.h"
#include "dynamics.h"
#include "node_func.h"

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

    //Declartion (and initialization) of some variables
    int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 500;
    double bound_min_y = 0.;
    double bound_max_y = 500;
    double bound_min_z = 0.;
    double bound_max_z = 500;
    
    std::ofstream ofile;
    MyNode *root = NULL;
    
    MyNode_val tmpnodeval = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 

    std::vector<MyParticle> particles_v;
    MyParticle tmpparticle = {0, 0, 0, 0, 0, 0, 0, 0};
    
    //Creation of a new MPI data type related to the MyParticle struct
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
    for(int add = 0; add<8; add++)
    {
        displs[add] = MPI_Aint_diff(address[add+1], address[0]);    
    }
    
    MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LOGICAL};
    MPI_Datatype MyParticle_mpi_t;
    MPI_Type_create_struct(8, blk_length, displs, types, &MyParticle_mpi_t);
    MPI_Type_commit(&MyParticle_mpi_t);

    //Creation of a new MPI data type related to MyNode_val struct;

    int blk_length_node[24] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint address_node[25];
    MPI_Get_address(&tmpnodeval, &address_node[0]);
    MPI_Get_address(&tmpnodeval.elements, &address_node[1]);
    MPI_Get_address(&tmpnodeval.index, &address_node[2]);
    MPI_Get_address(&tmpnodeval.depthflag, &address_node[3]);
    MPI_Get_address(&tmpnodeval.totalmass, &address_node[4]);
    MPI_Get_address(&tmpnodeval.COM_x, &address_node[5]);
    MPI_Get_address(&tmpnodeval.COM_y, &address_node[6]);
    MPI_Get_address(&tmpnodeval.COM_z, &address_node[7]);
    MPI_Get_address(&tmpnodeval.COM_vx, &address_node[8]);
    MPI_Get_address(&tmpnodeval.COM_vy, &address_node[9]);
    MPI_Get_address(&tmpnodeval.COM_vz, &address_node[10]);
    MPI_Get_address(&tmpnodeval.bound_min_x, &address_node[11]);
    MPI_Get_address(&tmpnodeval.bound_max_x, &address_node[12]);
    MPI_Get_address(&tmpnodeval.bound_min_y, &address_node[13]);
    MPI_Get_address(&tmpnodeval.bound_max_y, &address_node[14]);
    MPI_Get_address(&tmpnodeval.bound_min_z, &address_node[15]);
    MPI_Get_address(&tmpnodeval.bound_max_z, &address_node[16]);
    MPI_Get_address(&tmpnodeval.nwf, &address_node[17]);
    MPI_Get_address(&tmpnodeval.nef, &address_node[18]);
    MPI_Get_address(&tmpnodeval.swf, &address_node[19]);
    MPI_Get_address(&tmpnodeval.sef, &address_node[20]);
    MPI_Get_address(&tmpnodeval.nwb, &address_node[21]);
    MPI_Get_address(&tmpnodeval.neb, &address_node[22]);
    MPI_Get_address(&tmpnodeval.swb, &address_node[23]);
    MPI_Get_address(&tmpnodeval.seb, &address_node[24]);
    

    MPI_Aint displs_node[24];
    for(int add = 0; add<24; add++)
    {
        displs_node[add] = MPI_Aint_diff(address_node[add+1], address_node[0]);    
    }
    
    MPI_Datatype types_node[24] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, 
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,   
                MPI_LOGICAL, MPI_LOGICAL, MPI_LOGICAL, MPI_LOGICAL,
                MPI_LOGICAL, MPI_LOGICAL, MPI_LOGICAL, MPI_LOGICAL};

    MPI_Datatype MyNode_val_mpi_t;
    MPI_Type_create_struct(24, blk_length_node, displs_node, types_node, &MyNode_val_mpi_t);
    MPI_Type_commit(&MyNode_val_mpi_t);


    auto t0 = clk::now();
    if(prank == 0)
    {
        std::ifstream infile;
        int index = 0;
        infile.open("./input/disk.txt", std::ios::in);
        //Initialisation of the root node
        infile>>tmpparticle.x>>tmpparticle.y>>tmpparticle.z>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.vz>>tmpparticle.mass;
        tmpparticle.outside = false;
        root = initialize_node(tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, &index);
        particles_v.push_back(tmpparticle);
        while(infile>>tmpparticle.x)
        {
            infile>>tmpparticle.y>>tmpparticle.z>>tmpparticle.vx>>tmpparticle.vy>>tmpparticle.vz>>tmpparticle.mass;
            tmpparticle.outside = false;
            particles_v.push_back(tmpparticle);
            add_particle(root, tmpparticle, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, &index);
        }
        infile.close();
        std::cout<<"The final number of nodes in the root tree (the whole tree) is "<<index<<std::endl;
    }
    auto t1 = clk::now();
    n = root->elements;
    
    std::cout<<"The first creation of the tree took "<<second(t1 - t0).count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The root node has "<<root->elements << " elements for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The particles vector has " << particles_v.size() << " particles for process "<<prank<<" from the total of "<<psize<<std::endl;
    


    /////Tests/////
    int depth = 5;

    std::cout<<"There are "<<numNodesHeightK(root, depth)<<" nodes at depth "<<depth<<std::endl;
    
    std::vector<MyNode_val> serializedNode;
    serialize(root, serializedNode, depth+1);
    std::cout<<"The serialized node vector has " << serializedNode.size() << " particles for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[0].index <<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[1].index <<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[2].index <<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[3].index <<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[5].index <<std::endl;
    std::cout<<"The index of the serialized node vector " << serializedNode[30].index <<std::endl;
    
    MyNode *test_root = NULL;
    deSerialize(test_root, serializedNode);

    std::cout<<"There are "<<root->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<root->nwb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<root->swb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<root->swb->neb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<numNodesHeightK(root, 4)<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<numNodesHeightK(root, 5)<<" nodes at depth "<<depth<<std::endl;
    

    std::cout<<"There are "<<test_root->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<test_root->nwb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<test_root->swb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<test_root->swb->neb->elements<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<numNodesHeightK(test_root, 4)<<" nodes at depth "<<depth<<std::endl;
    std::cout<<"There are "<<numNodesHeightK(test_root, 5)<<" nodes at depth "<<depth<<std::endl;
    

    //////////////////
    

    //Declaration of variables for the actual computation
    /*double fx = 0., fy = 0., fz = 0;
    double ax = 0., ay = 0., az = 0;
    float dt = 0.1;

    auto ln = n/psize + (prank < n % psize ? 1 : 0);
    auto i_start = prank * ln + (prank < n % psize ? 0 : n % psize);
    auto i_end = i_start + ln;
    int recvcounts[psize];
    int displs_data[psize];
    MyParticle particles[n], particles_out[n];
    
    for(int aux = 0; aux < n; aux ++)
    {
        particles[aux] = particles_v[aux];
    }
    
    for(int c = 0; c<psize; c++)
    {
        recvcounts[c] = n/psize + (c < n % psize ? 1 : 0);
        displs_data[c] = c * recvcounts[c] + (c < n % psize ? 0 : n % psize);
        if(prank==0)
        {
            std::cout<<recvcounts[c]<< " " << displs_data[c]<<std::endl;
        }  
    }
    
    //Computation of new positions
    if(prank==0){ofile.open("./output/diskout.txt", std::ios::out);}
    for(int step = 0; step<1000; step++)
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
        }
        MPI_Gatherv(particles+i_start, ln, MyParticle_mpi_t, particles_out, recvcounts, displs_data, MyParticle_mpi_t, 0, MPI_COMM_WORLD );
        
        for (int aux = 0; aux < n; aux++){particles[aux]=particles_out[aux];}
        
        MPI_Bcast(particles, n, MyParticle_mpi_t, 0, MPI_COMM_WORLD);
        
        free_node(root);
        root = NULL;
        if(prank==0){ofile<<particles[0].x<<" "<<particles[0].y<<" "<<particles[0].z<<std::endl;}    
        root = initialize_node(particles[0], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        for(int p = 1; p < n; p++)
        {    
            if(prank==0){ofile<<particles[p].x<<" "<<particles[p].y<<" "<<particles[p].z<<std::endl;}    
            add_particle(root, particles[p], bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z);
        }
        if(prank==0){ofile<<step<<std::endl;}
        
    }
    if(prank==0){ofile.close();}*/
    second elapsed = clk::now() - t0;
    std::cout<<"The remaining number of particles in the particles vector is= "<<n <<"for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The number of particles in the tree is= "<<root->elements << "for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"The large loop with steps takes "<<elapsed.count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<"End of program"<< std::endl;

    MPI_Type_free(&MyParticle_mpi_t);
    MPI_Finalize();
    return 0;
}



