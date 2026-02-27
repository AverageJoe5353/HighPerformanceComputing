import numpy as np
import matplotlib.pyplot as plt
import os

os.system('./runSpaceHash')

N_grid = 100




# Open position file 
dataunsorted = np.loadtxt('spatialHashGridDemo_unsorted.txt')
datasorted = np.loadtxt('spatialHashGridDemo_sorted.txt')

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

for j, data in enumerate([dataunsorted, datasorted]):
    # Make a 10x10 grid of squares
    axs[j].vlines(np.arange(0, 11), 0, 10)
    axs[j].hlines(np.arange(0, 11), 0, 10)

    id = data[:,0]
    x = data[:,1]
    y = data[:,2]
    cellX = data[:,3]
    cellY = data[:,4]
    cellKey = data[:,5]

    # Plot the particles
    axs[j].scatter(x, y, c=cellKey, s=25)

    for i in range(len(x)):
        # Write the cell key in top corner, colour coded by number
        axs[j].text(x[i]+0.075, y[i]+0.25, str(int(cellKey[i])), fontsize=12, color='black')

        # Write the cell coordinates
        axs[j].text(x[i]-0.1, y[i]-0.4, str(int(cellX[i])) + ', ' + str(int(cellY[i])), fontsize=12)


plt.show()