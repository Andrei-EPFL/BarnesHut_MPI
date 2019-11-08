#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric> 
#include <chrono>
#include "update_node.h"
#include "dynamics.h"
#include "node_func.h"

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using time_point = std::chrono::time_point<clk>;

#include <mpi.h>

#define DEPTH_DEF 5
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
        
        /*//Debugging
        int n_elem_flagged =0 ;
        std::cout<<"There are " << numNodesElemFlagged(root_local, &n_elem_flagged)<<" flagged nodes ";
        std::cout<<"and " <<n_elem_flagged<<" elements\n";
        
        int n_elem_depthk =0 ;
        std::cout<<"There are " << numNodesElemHeightK(root_local, depth_local, &n_elem_depthk)<<"  nodes at depth "<< depth_local;
        std::cout<<" and " <<n_elem_depthk<<" elements\n";
        
        int n_elem_flagged_1 =0 ;
        std::cout<<"There are " << numNodesHeightDepthFlag(root_local, &n_elem_flagged_1) <<" flagged nodes ";
        std::cout<<"and " << n_elem_flagged_1<< " elements\n\n\n";
        //////////
        */
        free_node(root_local);
        root_local = NULL;
    }
    /////// END OF MOST OF THE SEQUENTIAL PART
    ////////////////////////////////////////////////////////////////////////////////////////////////

    auto t1 = clk::now();
    std::cout<<"prank="<<prank<<"-"<<psize<<": The first creation of the tree took "<<second(t1 - t0).count() << " seconds"<<std::endl;

    //// Send a global tree to everyone.
    MPI_Bcast(&n_serializednode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (prank != 0)
    {
        serializedNode_v.resize(n_serializednode);
        count_part_local.resize(psize);
    }
    MPI_Bcast(serializedNode_v.data(), serializedNode_v.size(), MyNode_val_mpi_t, 0, MPI_COMM_WORLD);
    
    //std::cout<<"prank="<<prank<<"-"<<psize<<": The serialized node vector has " << n_serializednode << " elements"<<std::endl; 
    //std::cout<<"prank="<<prank<<"-"<<psize<<": The serialized node vector has " << serializedNode_v.size() << " elements"<<std::endl;    
    
    deSerialize(root, serializedNode_v);
    //std::cout<<"prank="<<prank<<"-"<<psize<<": The serialized node vector after DeSerialization has " << serializedNode_v.size() << " elements"<<std::endl;
    ////

    ////Send particles to everyone.
    MPI_Bcast(count_part_local.data(), count_part_local.size(), MPI_INT, 0, MPI_COMM_WORLD);
    particles_v.resize(count_part_local[prank]);
    MPI_Scatterv(particles_v_local.data(), count_part_local.data(), displs_part_local.data(), MyParticle_mpi_t, particles_v.data(), count_part_local[prank], MyParticle_mpi_t, 0, MPI_COMM_WORLD);

    //std::cout<<"prank="<<prank<<"-"<<psize<<": The number of particles is = " << count_part_local[prank] <<std::endl;    
    //std::cout<<"prank="<<prank<<"-"<<psize<<": The number of particles is = " << particles_v.size() <<std::endl;
    //std::cout<<"prank="<<prank<<"-"<<psize<<": The root node has "<<root->elements << " elements"<<std::endl;
    ////
    ////root contains the common (global) tree; particles_v is a vector containing particles depending on the process; 
    
    ////Compute the local tree starting from the common tree.
    int index = 10000;
    for(unsigned int i = 0; i < particles_v.size(); i++)
    {
        add_particle_locally(root, particles_v[i], prank, &index);
    }
     
    ////////////////////////////////////////////////////// 
    //Declaration of variables for the actual computation
    std::vector<double> fx(particles_v.size());
    std::vector<double> fy(particles_v.size());
    std::vector<double> fz(particles_v.size());
    
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(fz.begin(), fz.end(), 0);

    std::vector<std::vector<double>> fx_send(psize);
    std::vector<std::vector<double>> fy_send(psize);
    std::vector<std::vector<double>> fz_send(psize);
    
    std::vector<std::vector<double>> fx_recv(psize);
    std::vector<std::vector<double>> fy_recv(psize);
    std::vector<std::vector<double>> fz_recv(psize);
    
    double ax = 0., ay = 0., az = 0;
    float dt = 0.01;

    std::vector<std::vector<int>> mat_size_particles(psize);
    for(unsigned i = 0; i < mat_size_particles.size(); i++)
    {
        mat_size_particles[i].resize(psize);
    }
    std::vector<std::vector<MyParticle>> mat_particles_send(psize);
    std::vector<std::vector<MyParticle>> mat_particles_recv(psize);

    std::vector<int> counts_size_recv(psize);
    std::fill(counts_size_recv.begin(), counts_size_recv.end(), psize);

    std::vector<int> displacements_size(psize);
    for(unsigned i = 0; i < displacements_size.size(); i++)
    {displacements_size[i]= psize*i;}
    std::vector<int> sizes_recv(psize*psize);
    
    std::vector<MPI_Request> request_send_part;
    std::vector<MPI_Request> request_recv_part;
    MPI_Request request_send_part_1;
    MPI_Request request_recv_part_1;
        
    std::vector<MPI_Request> request_send_force;
    std::vector<MPI_Request> request_recv_force;
    MPI_Request request_send_force_1;
    MPI_Request request_recv_force_1;

    ////// The begining of the main computational part
    if(prank==0){ofile.open("./output/diskout.txt", std::ios::out);}
    for(int step = 0; step<1000; step++)
    {
        //Computation of forces 
        std::fill(fx.begin(), fx.end(), 0);
        std::fill(fy.begin(), fy.end(), 0);
        std::fill(fz.begin(), fz.end(), 0);
        for(int p = 0; p < psize; p++) {mat_particles_send[p].clear(); mat_particles_send[p].resize(0);}  
        for(unsigned int i = 0; i < particles_v.size(); i++)
        {
            if(particles_v[i].outside == false)
            {
                compute_force_partially(root, particles_v[i], &fx[i], &fy[i], &fz[i], mat_particles_send);
            }   
        }
        ////Sending to everyone the sizes that they expect to get and send.
        for(int p = 0; p < psize; p++) {mat_size_particles[p].clear();mat_size_particles[p].resize(psize);}
        for(int p = 0; p < psize; p++)
        {
            mat_size_particles[prank][p] = mat_particles_send[p].size();
        }
        
        MPI_Gatherv(mat_size_particles[prank].data(), mat_size_particles[prank].size(), MPI_INT, sizes_recv.data(), counts_size_recv.data(), displacements_size.data(), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(sizes_recv.data(), sizes_recv.size(), MPI_INT, 0, MPI_COMM_WORLD);
    
        for(int j = 0; j < psize; j++)
        { 
            for(int p = 0; p < psize; p++)
            {
                mat_size_particles[j][p] = sizes_recv[j*psize+p];
                //std::cout<<mat_size_particles[j][p]<<" ";
            }
            //std::cout<<std::endl;
        }
        ////
        ////Transfer of particles
        ////
        int n_recv = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[p][prank]!=0)
            {   
                n_recv++;
                mat_particles_recv[p].resize(mat_size_particles[p][prank]);
                MPI_Irecv(mat_particles_recv[p].data(), mat_size_particles[p][prank], MyParticle_mpi_t, p, 0, MPI_COMM_WORLD, &request_recv_part_1);
                request_recv_part.push_back(request_recv_part_1);
            }
        }
        
        int n_send = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[prank][p]!=0)
            {
                n_send++;
                MPI_Isend(mat_particles_send[p].data(), mat_size_particles[prank][p], MyParticle_mpi_t, p, 0, MPI_COMM_WORLD, &request_send_part_1);  
                request_send_part.push_back(request_send_part_1);
            }
        }
        
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_recv-"<<n_recv<<std::endl;
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_send-"<<n_send<<std::endl;

        MPI_Waitall(request_send_part.size(), request_send_part.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(request_recv_part.size(), request_recv_part.data(), MPI_STATUSES_IGNORE);
        ////
        ////Computation of forces for the received particles;
        ////
        for (int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[p][prank]!=0)
            {
                fx_send[p].resize(mat_particles_recv[p].size());
                fy_send[p].resize(mat_particles_recv[p].size());
                fz_send[p].resize(mat_particles_recv[p].size());

                for(unsigned int i = 0; i < mat_particles_recv[p].size(); i++)
                {
                    std::fill(fx_send[p].begin(), fx_send[p].end(), 0);
                    std::fill(fy_send[p].begin(), fy_send[p].end(), 0);
                    std::fill(fz_send[p].begin(), fz_send[p].end(), 0);
                    if(mat_particles_recv[p][i].outside == false)
                    {
                        compute_force(root, mat_particles_recv[p][i], &fx_send[p][i], &fy_send[p][i], &fz_send[p][i]);
                    }   
                }
            }            
        }
        ////
        ////Transfer of forces
        ////
        //// fx
        int n_send_fx = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[p][prank]!=0)
            {   
                n_send_fx++;
                MPI_Isend(fx_send[p].data(), mat_size_particles[p][prank], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_send_force_1);  
                request_send_force.push_back(request_send_force_1);                
            }
        }

        int n_recv_fx = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[prank][p]!=0)
            {
                n_recv_fx++;
                fx_recv[p].resize(mat_size_particles[prank][p]);
                MPI_Irecv(fx_recv[p].data(), mat_size_particles[prank][p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_recv_force_1);
                request_recv_force.push_back(request_recv_force_1);
            }
        }
        
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_recv-"<<n_recv_fx<<std::endl;
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_send-"<<n_send_fx<<std::endl;

        MPI_Waitall(request_send_force.size(), request_send_force.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(request_recv_force.size(), request_recv_force.data(), MPI_STATUSES_IGNORE);

        //// fy
        int n_send_fy = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[p][prank]!=0)
            {   
                n_send_fy++;
                MPI_Isend(fy_send[p].data(), mat_size_particles[p][prank], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_send_force_1);  
                request_send_force.push_back(request_send_force_1);                
            }
        }

        int n_recv_fy = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[prank][p]!=0)
            {
                n_recv_fy++;
                fy_recv[p].resize(mat_size_particles[prank][p]);
                MPI_Irecv(fy_recv[p].data(), mat_size_particles[prank][p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_recv_force_1);
                request_recv_force.push_back(request_recv_force_1);
            }
        }
        
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_recv-"<<n_recv_fy<<std::endl;
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_send-"<<n_send_fy<<std::endl;

        MPI_Waitall(request_send_force.size(), request_send_force.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(request_recv_force.size(), request_recv_force.data(), MPI_STATUSES_IGNORE);

        //// fz
        int n_send_fz = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[p][prank]!=0)
            {   
                n_send_fz++;
                MPI_Isend(fz_send[p].data(), mat_size_particles[p][prank], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_send_force_1);  
                request_send_force.push_back(request_send_force_1);                
            }
        }

        int n_recv_fz = 0;
        for(int p = 0; p < psize; p++)
        {
            if(p!=prank && mat_size_particles[prank][p]!=0)
            {
                n_recv_fz++;
                fz_recv[p].resize(mat_size_particles[prank][p]);
                MPI_Irecv(fz_recv[p].data(), mat_size_particles[prank][p], MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &request_recv_force_1);
                request_recv_force.push_back(request_recv_force_1);
            }
        }
        
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_recv-"<<n_recv_fz<<std::endl;
        //std::cout<<"prank="<<prank<<"-"<<psize<<": n_send-"<<n_send_fz<<std::endl;

        MPI_Waitall(request_send_force.size(), request_send_force.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(request_recv_force.size(), request_recv_force.data(), MPI_STATUSES_IGNORE);
        
        /// Computation of the total force;
        for(unsigned int i = 0; i < particles_v.size(); i++)
        {
            if(particles_v[i].outside == false)
            {
                double fx_tmp = 0.;
                double fy_tmp = 0.;
                double fz_tmp = 0.;
                
                for(int p = 0; p < psize; p++)
                {
                    if(p!=prank && mat_size_particles[prank][p]!=0)
                    {
                        fx_tmp+= fx_recv[p][i];
                        fy_tmp+= fy_recv[p][i];
                        fz_tmp+= fz_recv[p][i];
                    }
                }
                fx[i]+= fx_tmp;
                fy[i]+= fy_tmp;
                fz[i]+= fz_tmp;
            }   
        //std::cout<<"THE FORCES at step= "<<step<< " for particle " <<i<< " are "<<fx[i]<<" "<<fy[i]<<" "<<fz[i]<<std::endl;
        }

        
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        //Computation of new positions
        for(unsigned int i = 0; i < particles_v.size(); i++)
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
        }
    }        
    if(prank==0){ofile.close();}
    ////// The ending of the main computational part

    second elapsed = clk::now() - t1;
    std::cout<<"prank="<<prank<<"-"<<psize<<": The remaining number of particles in the particles vector is= "<<particles_v.size() << std::endl;
    std::cout<<"prank="<<prank<<"-"<<psize<<": The large loop with steps takes "<<elapsed.count() << " seconds"<<std::endl;
    std::cout<<"prank="<<prank<<"-"<<psize<<": End of program"<< std::endl<<std::endl;

    MPI_Type_free(&MyParticle_mpi_t);
    MPI_Type_free(&MyNode_val_mpi_t);
    MPI_Finalize();
    return 0;
}



