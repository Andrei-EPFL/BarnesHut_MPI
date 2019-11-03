import numpy as np

def disk(orig, R, eps, N, omega):
    bodies = []
    x = R*np.random.rand(N)+orig[0]
    y = R*np.random.rand(N)+orig[1]
    for n in range(N):
        ##x = random.uniform(-R + orig[0], R + orig[0])
        ##y = random.uniform(-R + orig[1], R + orig[1])

        r = np.sqrt((x[n]-orig[0])**2 + (y[n]-orig[1])**2 +eps**2)

        if (r <=4*R):
            v = r * omega
            vx = -(y[n]-orig[1]) * v / r
            vy = (x[n]-orig[0]) * v / r 
            bodies.append((x[n],y[n], vx, vy, 1e11)) 
    return bodies


#bodies = disk([-1000, -1000, 0], 1000, 1, 1000, 0.000025)
#bodies += disk([1000, 1000, 0], 1000, 1, 1000, 0.000025)

bodies = disk([250, 250], 50, 0.01, 1000, 0)
#assivebody = [250, 250, 0, 0, 10e14]
#bodies += disk([1750, 1750], 100, 1, 100, 0)

#bodies += disk([400, 400], 200, 1, 50, 0.000025)

#bodies = disk([-1000, -1000, 1000], 1000, 1, 1000, 0.000025)
#bodies += disk([1000, 1000, -1000], 1000, 1, 1000, 0.000025)

#bodies = disk([-5000, -5000, 2000], 5000, 1, 25000, 0.0000025)
#bodies += disk([5000, 5000, -5000], 5000, 1, 25000, -0.0000025)

f = open("disk.txt", "w")
for b in bodies:
    for c in b:
        f.write(str(c) + " ")
    f.write("\n")
