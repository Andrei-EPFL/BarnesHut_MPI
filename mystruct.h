#ifndef MYSTRUCT_H
#define MYSTRUCT_H

struct MyParticle
{
    double x, y, vx, vy;
    double mass;
};

struct MyNode
{
    int elements;
    double totalmass, COM_x, COM_y, COM_vx, COM_vy;
    
    double bound_min_x;
    double bound_max_x;
    
    double bound_min_y;
    double bound_max_y;

    MyNode *nw, *ne, *sw, *se;
};

#endif