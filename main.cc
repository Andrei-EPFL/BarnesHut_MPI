#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <chrono>
#include "update_node.h"
#include "dynamics.h"
#include "node_func.h"

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

#include <mpi.h>

#define DEPTH_DEF 3
int main()
{
    MPI_Init(NULL, NULL);
    int prank, psize;
    MPI_Comm_rank(MPI_COMM_WORLD, &prank);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);

    //Declartion (and initialization) of some variables
    MyNode *root = NULL;
    std::vector<MyNode_val> serializedNode_v;
    std::vector<MyParticle> particles_v;
    
    int n_serializednode = 0;     
    double bound_min_x = 0.;
    double bound_max_x = 500;
    double bound_min_y = 0.;
    double bound_max_y = 500;
    double bound_min_z = 0.;
    double bound_max_z = 500; 
    MyNode_val tmpnodeval[2]; 
    MyParticle tmpparticle[2];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //Creation of a new MPI data type related to the MyParticle struct
    const int n_elements = 10;
    int blk_length[n_elements] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint address[n_elements+1];
    MPI_Get_address(&tmpparticle[0], &address[0]);
    MPI_Get_address(&tmpparticle[0].x, &address[1]);
    MPI_Get_address(&tmpparticle[0].y, &address[2]);
    MPI_Get_address(&tmpparticle[0].z, &address[3]);
    MPI_Get_address(&tmpparticle[0].vx, &address[4]);
    MPI_Get_address(&tmpparticle[0].vy, &address[5]);
    MPI_Get_address(&tmpparticle[0].vz, &address[6]);
    MPI_Get_address(&tmpparticle[0].mass, &address[7]);
    MPI_Get_address(&tmpparticle[0].node_index, &address[8]);
    MPI_Get_address(&tmpparticle[0].proc_rank, &address[9]);
    MPI_Get_address(&tmpparticle[0].outside, &address[10]);

    MPI_Aint displs[n_elements];
    for(int add = 0; add<n_elements; add++)
    {
        displs[add] = MPI_Aint_diff(address[add+1], address[0]);    
    }
    
    MPI_Datatype types[n_elements] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                             MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
                             MPI_DOUBLE, MPI_INT, MPI_INT, MPI_CXX_BOOL};

    MPI_Datatype MyParticle_mpi_temp_t;
    MPI_Type_create_struct(n_elements, blk_length, displs, types, &MyParticle_mpi_temp_t);
    MPI_Datatype MyParticle_mpi_t;

    MPI_Aint extent, address_second;
    MPI_Get_address(tmpparticle + 1, &address_second);
    extent = MPI_Aint_diff(address_second, address[0]);
    MPI_Type_create_resized(MyParticle_mpi_temp_t, 0, extent, &MyParticle_mpi_t);
    MPI_Type_commit(&MyParticle_mpi_t);

    //Creation of a new MPI data type related to MyNode_val struct;
    const int n_elements_node = 25;
    int blk_length_node[n_elements_node] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint address_node[n_elements_node+1];
    MPI_Get_address(&tmpnodeval[0], &address_node[0]);
    MPI_Get_address(&tmpnodeval[0].elements, &address_node[1]);
    MPI_Get_address(&tmpnodeval[0].index, &address_node[2]);
    MPI_Get_address(&tmpnodeval[0].depthflag, &address_node[3]);
    MPI_Get_address(&tmpnodeval[0].proc_rank, &address_node[4]);
    MPI_Get_address(&tmpnodeval[0].totalmass, &address_node[5]);
    MPI_Get_address(&tmpnodeval[0].COM_x, &address_node[6]);
    MPI_Get_address(&tmpnodeval[0].COM_y, &address_node[7]);
    MPI_Get_address(&tmpnodeval[0].COM_z, &address_node[8]);
    MPI_Get_address(&tmpnodeval[0].COM_vx, &address_node[9]);
    MPI_Get_address(&tmpnodeval[0].COM_vy, &address_node[10]);
    MPI_Get_address(&tmpnodeval[0].COM_vz, &address_node[11]);
    MPI_Get_address(&tmpnodeval[0].bound_min_x, &address_node[12]);
    MPI_Get_address(&tmpnodeval[0].bound_max_x, &address_node[13]);
    MPI_Get_address(&tmpnodeval[0].bound_min_y, &address_node[14]);
    MPI_Get_address(&tmpnodeval[0].bound_max_y, &address_node[15]);
    MPI_Get_address(&tmpnodeval[0].bound_min_z, &address_node[16]);
    MPI_Get_address(&tmpnodeval[0].bound_max_z, &address_node[17]);
    MPI_Get_address(&tmpnodeval[0].nwf, &address_node[18]);
    MPI_Get_address(&tmpnodeval[0].nef, &address_node[19]);
    MPI_Get_address(&tmpnodeval[0].swf, &address_node[20]);
    MPI_Get_address(&tmpnodeval[0].sef, &address_node[21]);
    MPI_Get_address(&tmpnodeval[0].nwb, &address_node[22]);
    MPI_Get_address(&tmpnodeval[0].neb, &address_node[23]);
    MPI_Get_address(&tmpnodeval[0].swb, &address_node[24]);
    MPI_Get_address(&tmpnodeval[0].seb, &address_node[25]);
    
    MPI_Aint displs_node[n_elements_node];
    for(int add = 0; add<n_elements_node; add++)
    {
        displs_node[add] = MPI_Aint_diff(address_node[add+1], address_node[0]);    
    }
    
    MPI_Datatype types_node[n_elements_node] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, 
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,   
                MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_CXX_BOOL,
                MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_CXX_BOOL, MPI_CXX_BOOL};

    MPI_Datatype MyNode_val_mpi_temp_t;
    MPI_Type_create_struct(n_elements_node, blk_length_node, displs_node, types_node, &MyNode_val_mpi_temp_t);
    MPI_Datatype MyNode_val_mpi_t;

    MPI_Aint extent_node, address_second_node;
    MPI_Get_address(tmpnodeval + 1, &address_second_node);
    extent_node = MPI_Aint_diff(address_second_node, address_node[0]);
    MPI_Type_create_resized(MyNode_val_mpi_temp_t, 0, extent_node, &MyNode_val_mpi_t);
    MPI_Type_commit(&MyNode_val_mpi_t);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    auto t0 = clk::now();
    
    std::vector<int> displs_part_local = {0};
    std::vector<int> count_part_local;
    std::vector<MyParticle> particles_v_local;
    std::ofstream ofile;
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ///// Sequential part begin
    if(prank == 0)
    {
        std::ifstream infile;
        int index_local = 0;
        int depth_local = DEPTH_DEF;
        int ln_local = 0;
        int n_nodes_depth_k_leaves = 0;

        MyParticle particle_local;
        MyNode *root_local = NULL;

        //Initialisation of the root_local node and particles_v_local
        infile.open("./input/disk.txt", std::ios::in);
        infile>>particle_local.x>>particle_local.y>>particle_local.z>>particle_local.vx>>particle_local.vy>>particle_local.vz>>particle_local.mass;
        particle_local.outside = false; particle_local.node_index = -1; particle_local.proc_rank = -1;
        root_local = initialize_node(particle_local, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, &index_local);
        particles_v_local.push_back(particle_local);
        
        while(infile>>particle_local.x)
        {
            infile>>particle_local.y>>particle_local.z>>particle_local.vx>>particle_local.vy>>particle_local.vz>>particle_local.mass;
            particles_v_local.push_back(particle_local);
            add_particle(root_local, particle_local, bound_min_x, bound_max_x, bound_min_y, bound_max_y, bound_min_z, bound_max_z, &index_local);
        }
        infile.close();
        /////////////////////////////////

        numNodeDepthKandLeaves(root_local, depth_local, &n_nodes_depth_k_leaves); //Flags the nodes that are at depth=depth_local and the nodes that are leaves having a depth<depth_local.
        
        //Debugging
        int n_elem_flagged =0 ;
        std::cout<<"There are " << numNodesElemFlagged(root_local, &n_elem_flagged)<<" flagged nodes ";
        std::cout<<"and " <<n_elem_flagged<<" elements\n";
        
        int n_elem_depthk =0 ;
        std::cout<<"There are " << numNodesElemHeightK(root_local, depth_local, &n_elem_depthk)<<"  nodes at depth "<< depth_local;
        std::cout<<" and " <<n_elem_depthk<<" elements\n";
        
        int n_elem_flagged_1 =0 ;
        std::cout<<"There are " << numNodesHeightDepthFlag(root_local, &n_elem_flagged_1) <<" flagged nodes ";
        std::cout<<"and " << n_elem_flagged_1<< " elements\n";
        //////////

        unsigned int n_tmp_part_count = 0;
        unsigned int n_tmp_part_displs = 0;
        std::cout<<"The number of nodes per process: ";
        for(int p = 0; p < psize; p++)
        {
            ln_local = n_nodes_depth_k_leaves/psize + (p < n_nodes_depth_k_leaves % psize ? 1 : 0);
            std::cout<<ln_local<<" ";
            linkParticlesNodepRank(root_local, depth_local, particles_v_local, p, &ln_local, &n_tmp_part_count);
            
            count_part_local.push_back(n_tmp_part_count);
            n_tmp_part_displs = n_tmp_part_count + n_tmp_part_displs;
            displs_part_local.push_back(n_tmp_part_displs);
            n_tmp_part_count = 0;
        }
        std::cout<<"\n";

        displs_part_local.erase(displs_part_local.end()-1);
        if(displs_part_local.size() != (unsigned int)psize && count_part_local.size() != (unsigned int)psize){std::cout<<"ERROR: The number of displacements is not correct\n";}
        if(n_tmp_part_displs != particles_v_local.size()) {std::cout<<"ERROR: There is an issue with the repartization of nodes (and particles) to processes\n";}
       
        std::sort(particles_v_local.begin(), particles_v_local.end(), compareByprank);
        
        serialize(root_local, serializedNode_v, depth_local);
        n_serializednode = serializedNode_v.size();
        
        std::cout<<"\n\nThe final number of nodes in the root tree (the whole tree) is "<<index_local<<std::endl;
        std::cout<<"The root node has "<<root_local->elements << " elements for process "<<prank<<" from the total of "<<psize<<std::endl;
        std::cout<<"The particles vector has " << particles_v_local.size() << " particles for process "<<prank<<" from the total of "<<psize<<std::endl; 
        std::cout<<"There are "<<n_nodes_depth_k_leaves<<" nodes and leaves at depth "<<depth_local<<std::endl;        
        std::cout<<"The serialized node vector has " << n_serializednode << " nodes for process "<<prank<<" from the total of "<<psize<<std::endl<<std::endl;     
        
        //Debugging
        n_elem_flagged =0 ;
        std::cout<<"There are " << numNodesElemFlagged(root_local, &n_elem_flagged)<<" flagged nodes ";
        std::cout<<"and " <<n_elem_flagged<<" elements\n";
        
        n_elem_depthk =0 ;
        std::cout<<"There are " << numNodesElemHeightK(root_local, depth_local, &n_elem_depthk)<<"  nodes at depth "<< depth_local;
        std::cout<<" and " <<n_elem_depthk<<" elements\n";
        
        n_elem_flagged_1 =0 ;
        std::cout<<"There are " << numNodesHeightDepthFlag(root_local, &n_elem_flagged_1) <<" flagged nodes ";
        std::cout<<"and " << n_elem_flagged_1<< " elements\n\n\n";
        //////////

        free_node(root_local);
        root_local = NULL;
    }
    /////// END OF MOST OF THE SEQUENTIAL PART
    ////////////////////////////////////////////////////////////////////////////////////////////////

    auto t1 = clk::now();
    std::cout<<prank<<": The first creation of the tree took "<<second(t1 - t0).count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;

    MPI_Bcast(&n_serializednode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (prank != 0)
    {
        serializedNode_v.resize(n_serializednode);
        count_part_local.resize(psize);
    }
    MPI_Bcast(serializedNode_v.data(), serializedNode_v.size(), MyNode_val_mpi_t, 0, MPI_COMM_WORLD);
    
    std::cout<<prank<<": The serialized node vector has " << n_serializednode << " elements for process "<<prank<<" from the total of "<<psize<<std::endl; 
    std::cout<<prank<<": The serialized node vector has " << serializedNode_v.size() << " elements for process "<<prank<<" from the total of "<<psize<<std::endl;    
    
    deSerialize(root, serializedNode_v);
    std::cout<<prank<<": The serialized node vector after DeSerialization has " << serializedNode_v.size() << " elements for process "<<prank<<" from the total of "<<psize<<std::endl;
    
    //Debugging
    int n_elem_flagged =0 ;
    std::cout<<prank<<": There are " << numNodesElemFlagged(root, &n_elem_flagged)<<" flagged nodes ";
    std::cout<<"and " <<n_elem_flagged<<" elements\n";
        
    int n_elem_depthk =0 ;
    std::cout<<prank<<": There are " << numNodesElemHeightK(root, DEPTH_DEF, &n_elem_depthk)<<"  nodes at depth "<< DEPTH_DEF;
    std::cout<<"and " <<n_elem_depthk<<" elements\n";
        
    int n_elem_flagged_1 =0 ;
    std::cout<<prank<<": There are " << numNodesHeightDepthFlag(root, &n_elem_flagged_1) <<" flagged nodes ";
    std::cout<<"and " << n_elem_flagged_1<< " elements\n";
    //////////

    //MPI_Bcast(count_part_local.data(), count_part_local.size(), MPI_INT, 0, MPI_COMM_WORLD);
    //particles_v.resize(count_part_local[prank]);
    //MPI_Scatterv(particles_v_local.data(), count_part_local.data(), displs_part_local.data(), MyParticle_mpi_t, particles_v.data(), count_part_local[prank], MyParticle_mpi_t, 0, MPI_COMM_WORLD);

/*
    std::cout<<prank<<": The number of particles for process "<<prank<<" from the total of "<<psize << " is = " << count_part_local[prank] <<std::endl;    
    std::cout<<prank<<": The number of particles for process "<<prank<<" from the total of "<<psize << " is = " << particles_v.size() <<std::endl;
    std::cout<<prank<<": The root node has "<<root->elements << " elements for process "<<prank<<" from the total of "<<psize<<std::endl;
  */  
    //root contains the common tree; particles_v is a vector containing particles depending on the process; 
    
    //Compute the local tree starting from the common tree.
  /*  int index = 10000;
    for(unsigned int i = 0; i < particles_v.size(); i++)
    {
        add_particle_locally(root, particles_v[i], prank, &index);
    }
    
    //int ntmp= 0;
    //std::cout<<prank<<": PUUUUTEPUUUTEPUUUTE "<<numNodesHeightDepthFlag(root, &ntmp)<<" and "<<ntmp<<std::endl;
    //ntmp = 0;
   // std::cout<<prank<<": "<<numNodesHeightK_tmp(root, 4, &ntmp)<< " and " << ntmp<<std::endl;
    ////////////////////////////////////////////////////// 
    
    //Declaration of variables for the actual computation
    std::vector<double> fx(particles_v.size());
    std::vector<double> fy(particles_v.size());
    std::vector<double> fz(particles_v.size());
    
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    double ax = 0., ay = 0., az = 0;
    float dt = 0.1;

    std::vector<std::vector<int>> mat_size_particles(psize);
    for(unsigned i = 0; i < mat_size_particles.size(); i++)
    {
        mat_size_particles[i].resize(psize);
    }

    std::vector<std::vector<MyParticle>> mat_particles_send(psize);
  */  //std::vector<std::vector<MyParticle>> mat_particles_recv(psize);

    
    //if(prank==0){ofile.open("./output/diskout.txt", std::ios::out);}
 /*   for(int step = 0; step<1; step++)
    {
        //Computation of forces 
        std::fill(fx.begin(), fx.end(), 0);
        std::fill(fy.begin(), fy.end(), 0);
        std::fill(fz.begin(), fz.end(), 0);     
        for(unsigned int i = 0; i < particles_v.size(); i++)
        {
            if(particles_v[i].outside == false)
            {
                //std::cout<<particles_v[i].node_index<<" "<<particles_v[i].proc_rank<<std::endl;
                //if(compute_force_partially(root, particles_v[i], &fx[i], &fy[i], &fz[i], mat_particles_send)==10) {std::cout<<"futaifutai"<<std::endl;}
                //compute_force_partially(root, particles_v[i], &fx[i], &fy[i], &fz[i], mat_particles_send);
            }   
        }
*/
        /*for(int p = 0; p < psize; p++)
        {
            mat_size_particles[prank][p] = mat_particles_send[p].size();
            std::cout<< mat_particles_send[p].size()<<" ";
        }
        std::cout<<std::endl;*/
        /*MPI_Bcast(mat_size_particles[prank].data(), psize, MPI_INT, prank, MPI_COMM_WORLD);
    
        for(int j = 0; j < psize; j++)
        {        
            for(int p = 0; p < psize; p++)
            {
                std::cout<<mat_size_particles[j][p]<<" ";
            }
            std::cout<<std::endl;
        }*/
        //do communications and finish the computation of forces coming from nodes from other processes.
    
        //Computation of new positions
        /*for(unsigned int i = 0; i < particles_v.size(); i++)
        {
            ax = fx[i]/particles_v[i].mass;
            ay = fy[i]/particles_v[i].mass;
            az = fz[i]/particles_v[i].mass;
                
            particles_v[i].vx += ax * dt;
            particles_v[i].vy += ay * dt;
            particles_v[i].vz += az * dt;
                
            particles_v[i].x += particles_v[i].vx * dt;
            particles_v[i].y += particles_v[i].vy * dt;
            particles_v[i].z += particles_v[i].vz * dt;

            if(particles_v[i].x < bound_min_x || particles_v[i].y < bound_min_y || particles_v[i].x > bound_max_x || particles_v[i].y > bound_max_y || particles_v[i].z > bound_max_z || particles_v[i].z < bound_min_z)
            {
                particles_v[i].x = bound_min_x - 100; 
                particles_v[i].y = bound_min_y - 100;
                particles_v[i].z = bound_min_z - 100;
                particles_v[i].vx = particles_v[i].vy = particles_v[i].vz = particles_v[i].mass = 0;
                particles_v[i].outside = true;
            }
        }
        
        if(prank==0)
        {    
            for(unsigned int i = 0; i < particles_v.size(); i++)
            {ofile<<particles_v[i].x<<" "<<particles_v[i].y<<" "<<particles_v[i].z<<std::endl;}
            ofile<<step<<std::endl;
        }*/
    //}        
    //if(prank==0){ofile.close();}

    second elapsed = clk::now() - t0;
    std::cout<<prank<<": The remaining number of particles in the particles vector is= "<<particles_v.size() <<", for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<prank<<": The large loop with steps takes "<<elapsed.count() << " seconds for process "<<prank<<" from the total of "<<psize<<std::endl;
    std::cout<<prank<<": End of program"<< std::endl<<std::endl;

    MPI_Type_free(&MyParticle_mpi_t);
    MPI_Type_free(&MyNode_val_mpi_t);
    MPI_Finalize();
    return 0;
}



