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
   

if __name__== '__main__':
    main()