import matplotlib.pyplot as plt
import numpy as np



def main():
    kpc = 3.086e19
    Msun= 1.989e30
    Myear = 1e6*365*24*3600.

    m = 1./kpc
    kg=1./Msun
    s = 1./Myear

    G = ((6.674*10e-11)*m**3)/(kg*s*s)
    print(G)
    
    print("\n")

    y = [0, 0, 0, 0, 0, 0]
    t_s = [0.005, 0.005, 0.007, 0.006, 0.006, 0.0056]
    x = [1, 2, 4, 6, 8, 9]
    t_t = [52.48+t_s[0], 74.34+t_s[1], 47.07+t_s[2], 27.87+t_s[3], 23.84+t_s[4], 19.811+t_s[5] ]
    for i in x:
        y = t_t/t_t[0]

    pt.plot(x, y)


if __name__== '__main__':
    main()