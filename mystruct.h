#ifndef MYSTRUCT_H
#define MYSTRUCT_H

struct MyParticle
{
    double x, y, z, vx, vy, vz;
    double mass;
};

struct MyNode
{
    int elements;
    double totalmass, COM_x, COM_y, COM_z, COM_vx, COM_vy, COM_vz;
    
    double bound_min_x;
    double bound_max_x;
    
    double bound_min_y;
    double bound_max_y;

    double bound_min_z;
    double bound_max_z;


    MyNode *nwf, *nef, *swf, *sef;
    MyNode *nwb, *neb, *swb, *seb;
};

#endif