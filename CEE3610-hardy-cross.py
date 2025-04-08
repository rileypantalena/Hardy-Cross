import numpy as np

# calculate n and k for every pipe

n = 1.852

def kfunc(L, C, D):
    return (( 4.727 * L )) / ( C ** (1.852) * D ** (4.8704))

# kvalues = np.array([kfunc(1500, 100, 10/12), kfunc(400, 100, 10/12), kfunc(600, 120, 10/12), kfunc(800, 110, 8/12), kfunc(700, 110, 8/12), kfunc(1000, 90, 10/12), kfunc(750, 100, 10/12)])
# Qvalues = np.array([0.77574, 2.0, 0.784386, 0.2156135, 0.9632135, 1.4312265, 0.2156135])

# calculate head loss per pipe in loop 1 using the Hazen Williams method
kvalues1 = np.array([kfunc(1500, 100, 10/12), kfunc(400, 100, 10/12), kfunc(600, 120, 10/12), kfunc(800, 110, 8/12)])
Qvalues1 = np.array([0.77574*1.547, 2.0*1.547, 0.784386*1.547, 0.2156135*1.547])
headloss1 = kvalues1 * Qvalues1 ** n

sign1 = np.array([1, -1, -1, 1])
headlosssign1 = headloss1 * sign1

# calculate delta Q for loop 1
hloverQ1 = headlosssign1 / Qvalues1
deltaQ1 = (-headlosssign1.sum()) / (n * np.abs(hloverQ1.sum()))


# adjust the flows for all pipes in loop 1 using delta Q 1
updatedQ1 = Qvalues1 + deltaQ1


# repeat for loop 2:
# calculate head loss per pipe in loop 2 using the Hazen Williams method
kvalues2 = np.array([ kfunc(800, 110, 8/12), kfunc(700, 110, 8/12), kfunc(1000, 90, 10/12), kfunc(750, 100, 10/12)])
Qvalues2 = np.array([0.2156135*1.547, 0.9632135*1.547, 1.4312265*1.547, 0.2156135*1.547])
headloss2 = kvalues2 * Qvalues2 ** n

sign2 = np.array([-1, -1, 1, 1])
headlosssign2 = headloss2 * sign2

# calculate delta Q for loop 2
hloverQ2 = headlosssign2 / Qvalues2
deltaQ2 = (-headlosssign2.sum()) / (n * np.abs(hloverQ2.sum()))


# adjust the flows for all pipes in loop 1 using delta Q 2
updatedQ2 = Qvalues2 + deltaQ2


iterations = 100
tolerance = 0.01
iteration = 0

for i in range(iterations):

    headloss1 = kvalues1 * Qvalues1 ** n
    headlosssign1 = headloss1 * sign1

    # Calculate delta Q for loop 1
    hloverQ1 = headlosssign1 / Qvalues1
    deltaQ1 = (-headlosssign1.sum()) / (n * np.abs(hloverQ1.sum()))

    # Update Qvalues1
    Qvalues1_new = Qvalues1 + deltaQ1

    # Calculate head loss for loop 2
    headloss2 = kvalues2 * Qvalues2 ** n
    headlosssign2 = headloss2 * sign2

    # Calculate delta Q for loop 2
    hloverQ2 = headlosssign2 / (Qvalues2)
    deltaQ2 = (-headlosssign2.sum()) / (n * np.abs(hloverQ2.sum()))

    # Update Qvalues2
    Qvalues2_new = Qvalues2 + deltaQ2

    if abs(deltaQ1) < tolerance and abs(deltaQ2) < tolerance:
        break
    Qvalues1 = Qvalues1_new
    Qvalues2 = Qvalues2_new

    iteration += 1

print(f"Converged after {iteration} iterations")
print(deltaQ1)
print(deltaQ2)
print("Final Qvalues1:", Qvalues1)
print("Final Qvalues2:", Qvalues2)


# =^.^=
#  |-|_/



