import matplotlib.pyplot as plt
import math
import random

NUMBER_OF_PARTICLES = 700
SIGMA = 1
CUTOFF = 2.5 * SIGMA
MAX_DISPLACEMENT = 0.005 * SIGMA
TEMPERATURE = 1
DENSITY = 0.7
X, Y, Z = [], [], []
L = 10 * SIGMA

LJ_ENERGY, MC_STEPS = [], []
current_energy = 0


def initcontainer():
    while len(X) < NUMBER_OF_PARTICLES:
        x1 = random.uniform(0, L)
        y1 = random.uniform(0, L)
        z1 = random.uniform(0, L)

        X.append(x1)
        Y.append(y1)
        Z.append(z1)

        for i in range(len(X) - 1):
            if(distance(len(X) - 1, i) <= SIGMA): # otherwise potential will be positive
                X.pop()
                Y.pop()
                Z.pop()
                break


def distance(i, j):
    return math.sqrt((X[i] - X[j]) * (X[i] - X[j]) + (Y[i] - Y[j]) * (Y[i] - Y[j]) + (Z[i] - Z[j]) * (Z[i] - Z[j]))


def LJ_energy():
    energy = 0
    # Find distance between every pair of molecules
    for i in range(0, NUMBER_OF_PARTICLES):
        for j in range(i + 1, NUMBER_OF_PARTICLES):
            # Reduced form of LJ Energy = 4*[(sigma/r)^12 - (sigma/r)^6]
            dist = distance(i, j)
            if dist > CUTOFF:
                 # potential energy negligible, to avoid precision error
                continue
            invd = pow(SIGMA / dist, 6)
            energy += invd * (invd - 1)

    energy *= 4
    return energy


def findPotential(i):
    potential = 0
    for j in range(0, NUMBER_OF_PARTICLES):
        if i == j:
            continue
        dist = distance(i, j)
        if dist > CUTOFF:
            continue
        invd = pow(SIGMA / dist, 6)
        potential += invd * (invd - 1)
    return potential


def trial():
    global current_energy
    # random no from 0 to 699
    rn = random.randint(0, NUMBER_OF_PARTICLES - 1)
    # random real no from 0 to max_distance
    dist = random.uniform(0, MAX_DISPLACEMENT)

    initialPotential = findPotential(rn)
    initial = [X[rn], Y[rn], Z[rn]]
    X[rn] += random.uniform(0, L)
    Y[rn] += random.uniform(0, L)
    Z[rn] += random.uniform(0, L)

    # periodic boundary conditions
    if X[rn] > L:
        X[rn] -= L
    if Y[rn] > L:
        Y[rn] -= L
    if Z[rn] > L:
        Z[rn] -= L

    if X[rn] < 0:
        X[rn] += L
    if Y[rn] < 0:
        Y[rn] += L
    if Z[rn] < 0:
        Z[rn] += L

    finalPotential = findPotential(rn)
    delta = finalPotential - initialPotential
    if delta > 75:
        # Rejected, e^(-delta) > rand is neccesary
        X[rn], Y[rn], Z[rn] = initial
    else:
        current_energy += delta
        rand = random.uniform(0, 1)
        if delta > 0 and math.exp(-delta/TEMPERATURE) < rand:
            # Rejected
            X[rn], Y[rn], Z[rn] = initial
            current_energy -= delta


def mcs_cycle():
    # do 700 trials
    for i in range(NUMBER_OF_PARTICLES):
        trial()


if __name__ == "__main__":
    initcontainer()
    current_energy = LJ_energy()
    for i in range(1000):
        mcs_cycle()
        if i % 10 == 0:
            LJ_ENERGY.append(current_energy)
            MC_STEPS.append(i)

    # plotting the points
    plt.plot(MC_STEPS, LJ_ENERGY)

    # naming the axis
    plt.xlabel('MC Steps')
    plt.ylabel('LJ Energy')

    plt.title('LJ Energy VS MC Steps')

    # function to show the plot
    plt.show()