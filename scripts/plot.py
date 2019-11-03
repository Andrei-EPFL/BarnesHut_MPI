import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import math
from itertools import product, combinations
import sys

def plot_earth_moon():
    x, y = np.loadtxt("./output.dat", usecols=(0,1), unpack=True)
    plt.plot(x, y)
    plt.show()


def draw_cuboid(x, y, z, ax, c='b'):
    # draw cuboid
    max_diff = [x[1] - x[0], y[1] - y[0]]
    plots = []
    for s, e in combinations(np.array(list(product(x, y))), 2):
        if np.sum(np.abs(s-e) == max_diff) == 1 + len(max_diff) - np.count_nonzero(max_diff):
                plots.append(ax.plot(*zip(s, e), color=c))
    return plots

def draw_cuboids(cuboids, ax):
    plots = []
    for i, c in enumerate(cuboids):
        plots += draw_cuboid(*c, ax)
    return plots




def main():

    file_object = open(sys.argv[1],"r")
    a = file_object.readlines()
    positions = []
    complete_positions = []
    count = 0
    min_pos = math.inf
    max_pos = -math.inf
    for i in range(len(a)):
        #print(len(a[i].split()))
        if len(a[i].split()) < 2:
            positions_in_time_step = np.array(positions)
            if np.max(positions_in_time_step) > max_pos:
                max_pos = np.max(positions_in_time_step)
            if np.min(positions_in_time_step) < min_pos:
                min_pos = np.min(positions_in_time_step)
            complete_positions.append(np.array(positions_in_time_step))
            count += 1
            positions = []
        else:
            positions.append([float(x) for x in a[i].split()]) 


    fig = plt.figure()
    ax = fig.gca()

    for i in range(0, len(complete_positions)):
        ax.set_xlim(200,300)
        ax.set_ylim(200,300)

        p1 = ax.scatter(complete_positions[i][:,0], complete_positions[i][:,1], c='r', marker='o')
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        #ax.view_init(30, 0.5*i)
      
        plt.pause(0.0001)
        plt.cla()

if __name__== '__main__':
    main()