import os
import sys
import argparse
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from itertools import product
from statistics import mean
import pickle
import random
#----------------------------------making directories---------------------------------#


os.makedirs('data', exist_ok=True)
os.makedirs('graphs', exist_ok=True)
os.makedirs('params', exist_ok=True)

#----------------------------------functions---------------------------------#


def load_network(matrix, hourtime):
    matrix = matrix
    matrix = ((np.loadtxt("params/matrix_" + matrix +
                          '.txt', delimiter=" ")).astype(int)).tolist()
    init = np.loadtxt('params/init_{}.txt'.format(filename), delimiter=" ")
    unpe = np.loadtxt('params/perturb_{}.txt'.format(filename), delimiter=" ")
    post_exps = pickle.load(open("params/postexps_{}.p".format(filename), "rb"))
    return matrix, init, unpe, post_exps


def generate_network(graph, init, hourtime, postfix):
    #----------------------------------input setup---------------------------------#
    df = init  # initial values
    timecost = 3600*hourtime  # time in seconds
    n = 500  # time step
    t = np.linspace(0, timecost, n)  # time vector for initial simulation
    #---------------------------------ODE computing--------------------------------#
    df = odeint(tree_model, df, t)
    np.savetxt('data/solution_{}.txt'.format(new_filename), df, fmt='%1.4f')

    return df


def perturb(df, hourtime, perturbation):
    #----------------------perturbing system----------------------#
    timecost = 3600*hourtime
    n = 500
    t2 = np.linspace(0, timecost/4, n)
    dfo = df[-1]
    dfo = odeint(tree_model, dfo, t2)
    dfp = df[-1] + perturbation[::]
    dfp = odeint(tree_model, dfp, t2)

    #----------------------saving output for reusability----------------------#
    np.savetxt('data/steady_{}.txt'.format(new_filename), dfo, fmt='%1.4f')
    np.savetxt('data/perturbed_{}.txt'.format(new_filename), dfp, fmt='%1.4f')

    return dfo, dfp


def tree_model(df, t):
    '''
    this function reads in the initial conditions and solves the ordinary differntial equations.
    params:
    :df: (array) initial condition on y.
    :t: (array) the sequence of time points to solve for y.
    output:
    :ode_list: (array) values of y for each desired time in t, with y0 in the first row.
    '''
    #----------------------------------parameters---------------------------------#
    kmax = 1.2 * 10**(-2)  # mM/s
    kbas = 1.2 * 10**(-8)  # mM/s
    kdec = 1.2 * 10**(-5)  # 1/s
    alph = 2  # Hill coefficient
    kprot = 1  # mM
    #----------------------------------solving---------------------------------#
    ode_list = []
    for i in range(len(post_exps)):
        text = post_exps[i]
        if len(text) == 0:
            diffeq = 0
        else:
            stack = []
            for token in text:
                if token == 'and' or token == 'or':
                    a = stack.pop()
                    if a.isLogic == True:
                        if a.value > 0:
                            a = (df[a.index]**alph) / \
                                (kprot**alph + df[a.index]**alph)
                        elif a.value < 0:
                            a = (kprot**alph)/(kprot**alph + df[a.index]**alph)
                    else:
                        a = a.value
                    if len(stack) == 0:
                        if token == 'and':
                            b = Token(1)
                        if token == 'or':
                            b = Token(0)
                    else:
                        b = stack.pop()
                    if b.isLogic == True:
                        if b.value > 0:
                            b = (df[b.index]**alph) / \
                                (kprot**alph + df[b.index]**alph)
                        elif b.value < 0:
                            b = (kprot**alph)/(kprot**alph + df[b.index]**alph)
                    else:
                        b = b.value
                    if token == 'and':
                        stack.append(Token(a * b))
                    if token == 'or':
                        stack.append(Token(1 - (1 - a) * (1 - b)))
                elif (token.value == 1 or token.value == -1):
                    stack.append(token)
                else:
                    raise ValueError("unknown token {0}".format(token))
            diffeq = stack.pop().value
        diffeq *= kmax
        diffeq += kbas - (kdec * df[i])
        ode_list.append(diffeq)
    return ode_list


def graph_network(df, hourtime, dfp=None):
    columns = df.shape[1]
    if dfp is not None:
        t2 = np.linspace(0, (hourtime/4)/24, 500)
        #------------------graphing protein concentration of steady vs perturbed------------------#
        fig1 = plt.figure()
        fig1.subplots_adjust(hspace=0.4, wspace=0.4)
        fig1.set_size_inches(10, 4.8)
        scatterplot = fig1.add_subplot(211)
        scatterplotp = fig1.add_subplot(212)
        colors = ['r', 'y', 'g', 'c', 'b', 'm']
        lines = ['--', ':', '-.']
        graphs = []
        for r in product(lines, colors):
            graphs.append(r[1] + r[0])
        for i in range(columns):
            scatterplot.plot(t2, df[:, i], graphs[i %
                                                  18], label=str(chr(i+65))+"(t)")
            scatterplotp.plot(t2, dfp[:, i], graphs[i % 18],
                              label=str(chr(i+65))+"(t)")
        scatterplotp.set_xlabel('time (days)', fontsize=14)
        scatterplot.set_ylabel('protein concentration (' + u"\u03bcM)")
        scatterplotp.set_ylabel('protein concentration (' + u"\u03bcM)")
        scatterplot.set_title("Random Gene Regulatory Network of " +
                              str(columns) + " Genes", fontsize=16)
        scatterplotp.set_title(
            "Random Gene Regulatory Network with Perturbation", fontsize=16)
        if columns < 12:
            scatterplot.legend(loc=1)
            scatterplotp.legend(loc=1)
        scatterplot.set_xlim(0, (hourtime/4)/24)
        scatterplotp.set_xlim(0, (hourtime/4)/24)
        fig1.savefig(
            "graphs/peturbtimeseries_after_{}h_{}.png".format(int(hourtime/2), new_filename), dpi=100)
        plt.close()

        #----------------graphing distance between original and perturbed over time---------------#
        dist = abs(dfp[::] - dfo[::])
        rows = dist.shape[0]
        distances = []
        lle_distances = []
        for i in range(rows):
            distances.append(mean(dist[i]))
            lle_distances.append(np.linalg.norm(dist[i]))
        fig2 = plt.figure()
        displot = fig2.add_subplot(111)
        displot.plot(t2, distances)
        displot.set_title(
            "Distance Between Perturbed and Steady States", fontsize=16)
        displot.set_xlabel('time (days)', fontsize=14)
        displot.set_ylabel('distance', fontsize=14)
        fig2.savefig(
            "graphs/distafter_{}h_{}.png".format(int(hourtime/2), new_filename), dpi=100)
        np.savetxt(
            'data/distafter_{}h_{}.txt'.format(int(hourtime/2), new_filename), dist, fmt='%s')
        plt.close()

        #------------------------calculating largest lyapunov exponent [EXPERIMENTAL! not working]----------------------#
        # difference = mean(lle_distances)
        # perturb = np.linalg.norm(unpe)
        # lle = 1/n * np.log(difference/perturb)
        # np.savetxt('output/lle.txt', [lle], fmt='%s')
    else:
        t = np.linspace(0, hourtime, 500)
        #------------------------graphing protein concentration-----------------------#
        fig = plt.figure()
        scatterplot = fig.add_subplot(111)
        colors = ['r', 'y', 'g', 'c', 'b', 'm']
        lines = ['--', ':', '-.']
        graphs = []
        for r in product(lines, colors):
            graphs.append(r[1] + r[0])
        for i in range(columns):
            scatterplot.plot(t, df[:, i], graphs[i % 18],
                             label=str(chr(i+65))+"(t)")
        plt.xlabel('time (h)', fontsize=15)
        plt.ylabel('protein product concentration (µM)', fontsize=15)
        plt.title("Random Gene Regulatory Network of " +
                  str(columns) + " Genes", fontsize=16)
        if columns < 12:
            plt.legend(loc=2)
        plt.xlim(0, hourtime)
        fig.savefig("graphs/timeseries_{}.png".format(new_filename))
        plt.close()


class Token:
    def __init__(self, value, isLogic=False, index=None):
        self.value = value
        self.isLogic = isLogic
        self.index = index


if __name__ == '__main__':
    """
    the user's input will define the following:
    1. filename of matrix (without .txt)
    2. time of simulation in hours
    """
    parser = argparse.ArgumentParser(description='load_rgrn.py simulates a previously generated random gene regulatory network given'
                                     'the gene network matrix file and the simulation time of rerun')
    parser.add_argument('-g', '--numberOfGenes', required=True, dest='genes',
                        help='number of genes in the network [positive int]')
    parser.add_argument('-i', '--numberOfInteractions', required=True, dest='interactions',
                        help='number of interactions in network or connected or repressilator [nonnegative int, \'C\', or \'R\']')
    parser.add_argument('-t', '--hourtime', required=True, dest='hourtime',
                        help='simulation time in hours [positive int]')
    parser.add_argument('-s', '--seed', required=True, dest='seed',
                        help='seed for random number generator [positive int]')
    parser.add_argument('-n', '--newHourtime', required=True, dest='newhourtime',
                        help='simulation time of rerun in hours [positive int]')

    args = parser.parse_args()
    genes = int(args.genes)
    interactions = str(args.interactions)
    hourtime = int(args.hourtime)
    seed = int(args.seed)
    new_hourtime = int(args.newhourtime)

    filename = 'genes_{}_interactions_{}_hours_{}_seed_{}'.format(genes, interactions, hourtime, seed)
    new_filename = 'genes_{}_interactions_{}_hours_{}_seed_{}'.format(genes, interactions, new_hourtime, seed)
    random.seed(seed)

    graph, init, perturbation, post_exps = load_network(filename, hourtime)
    df = generate_network(graph, init, hourtime, post_exps)
    graph_network(df, hourtime)
    dfo, dfp = perturb(df, hourtime, perturbation)
    graph_network(dfo, hourtime, dfp)
