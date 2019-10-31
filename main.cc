#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "update_node.h"
#include "dynamics.h"

int main(int argc, char *argv[])
{
    if(argc!=2) return 0;
    std::ifstream infile;
    MyNode *root = NULL;
    int MM = atoi(argv[1]);
    int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 1000e6;
    double bound_min_y = 0.;
    double bound_max_y = 1000e6;

    infile.open("2particles.dat", std::ios::in);
    
    std::vector<MyParticle> particles;
    MyParticle tmpparticle;

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
    
    double fx = 0., fy = 0.;
    
    double ax = 0., ay = 0.;
    float dt = 3600.;
    
    std::ofstream outfile;
    outfile.open("output.dat", std::ios::out);
    
    

    std::cout<<"\nInitial position: "<<std::endl;
    std::cout<<"Earth (x, y)"<< std::endl;
    std::cout<<particles[0].x<<" "<<particles[0].y<<std::endl;

    std::cout<<"Moon (x, y)"<< std::endl;
    std::cout<<particles[1].x<<" "<<particles[1].y<<std::endl;  

    for(int step = 0; step<24*27; step++)
    {
        fx = fy = 0.;
        for(int i = 0; i < n; i++)
        {
            compute_force(root, particles[i], &fx, &fy);
            ax = fx/particles[i].mass;
            ay = fy/particles[i].mass;
            
            particles[i].vx += ax * dt;
            particles[i].vy += ay * dt;
            //std::cout<<ax<<" "<<ay<<std::endl;
            particles[i].x += particles[i].vx * dt;
            particles[i].y += particles[i].vy * dt;
            fx = fy = 0;     
        }
        //if(step%20==0)
        //{
            std::cout<<"\nStep: "<<step<<std::endl;
            std::cout<<"Earth (x, y)"<< std::endl;
            std::cout<<particles[0].x<<" "<<particles[0].y<<std::endl;

            std::cout<<"Moon (x, y)"<< std::endl;
            outfile<<particles[1].x<<" "<<particles[1].y<<std::endl;    
        //}
        free_node(root);
        root = NULL;
        root = initialize_node(particles[0], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
        add_particle(root, particles[1], bound_min_x, bound_max_x, bound_min_y, bound_max_y);
    }
    //std::cout<<root->elements<<std::endl;
    //std::cout<<root->sw->elements<<std::endl;
    //compute_force(root, particles[MM], &fx, &fy);
    //std::cout<<fx<<" "<<fy<<std::endl;
    std::cout<<MM<<std::endl;
    std::cout<<"end of program"<< std::endl;
    return 0;
}



