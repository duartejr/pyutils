import numpy as np

def idw(x, y, v, grid, power):
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            distance = np.sqrt((x-i)**2+(y-j)**2)
            if (distance**power).min() == 0:
                grid[i,j] = v[(distance**power).argmin()]
            else:
                total = np.sum(1/(distance**power))
                grid[i,j] = np.sum(v/(distance**power)/total)
    return grid