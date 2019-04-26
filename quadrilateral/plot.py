import numpy as np
from matplotlib import pyplot as plt

def main():
    node = np.loadtxt("result/node.txt")
    elem = np.loadtxt("result/element.txt", dtype=int)
    disp = np.loadtxt("result/displacement.txt")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for e in elem:
        t = np.array([node[e[i%4] - 1] for i in range(len(e) + 1)])
        ax.plot(t[:, 0], t[:, 1], color="k")

    for e in elem:
        t = np.array([node[e[i%4] - 1] for i in range(len(e) + 1)])
        s = np.array([disp[e[i%4] - 1] for i in range(len(e) + 1)])
        t += s
        ax.plot(t[:, 0], t[:, 1], color="r")
    plt.savefig("result/displace.png")

if __name__ == "__main__":
    main()