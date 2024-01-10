import numpy as np

radii = [0.5, 0.61, 0.8, 4.55]
Et = [2., 0.5, 1.5, 1.5]
vols = []
for i in range(len(radii)):
    if i == 0:
        vols.append(radii[i]*radii[i]*np.pi)
    else:
        vols.append((radii[i]*radii[i] - radii[i-1]*radii[i-1])*np.pi)

p = np.array([[0.5951, 0.0336, 0.1201, 0.2535],
              [0.2749, 0.0896, 0.2536, 0.3863],
              [0.1495, 0.0385, 0.3455, 0.4781],
              [0.0042, 0.00079, 0.00638, 0.9886]])

for i in range(4):
    p[i,:] *= vols[i] * Et[i]

print(p)