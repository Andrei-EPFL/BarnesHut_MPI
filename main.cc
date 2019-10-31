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
    //int n = 0;
    double bound_min_x = 0.;
    double bound_max_x = 128.;
    double bound_min_y = 0.;
    double bound_max_y = 128.;

    infile.open("particles.dat", std::ios::in);
    
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
    //n = root->elements;
    std::cout<<"The root node has "<<root->elements << " elements"<<std::endl;
    std::cout<<"The particles vector has " << particles.size() << " particles" << std::endl;
    
    double fx = 0., fy = 0.;
    
    /*double ax = 0., ay = 0.;
    float dt = 1.;
    for(int i = 0; i < n; i++)
    {
        compute_force(root, particles[i], &fx, &fy);
        ax = fx/particles[i].mass;
        ay = fy/particles[i].mass;

        particles[i].vx = ax * dt;
        particles[i].vy = ay * dt;

        particles[i].x = particles[i].vx * dt;
        particles[i].y = particles[i].vy * dt;     
    }*/

    compute_force(root, particles[MM], &fx, &fy);
    std::cout<<fx<<" "<<fy<<std::endl;
    std::cout<<root->nw->elements;
    std::cout<<root->se->elements;
    std::cout<<root->ne->elements;
    std::cout<<root->ne->ne->elements;    
    std::cout<<"end of program"<< std::endl;
    return 0;
}



