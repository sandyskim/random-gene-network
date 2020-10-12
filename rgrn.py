#-----------------------------------packages-----------------------------------#
from __future__ import division
import numpy as np
import scipy
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import random
import sys
from scipy.integrate import odeint
from itertools import product
from statistics import mean
import pickle
import math
import warnings
import time

#----------------------------------making directories---------------------------------#


os.makedirs('output', exist_ok=True)
os.makedirs('params', exist_ok=True)

#----------------------------------functions---------------------------------#


def initialize(genes, interactions):
    #----------------------------generate random matrix----------------------------#
    graph = np.zeros((genes, genes))
    if interactions == 'R':  # generate N-represillator network
        interactions = genes
        for aa in range(genes):
            graph[aa, (aa+1) % genes] = -1
    if interactions == 'C':  # generated connected network with 10% edge density
        interactions = int(genes**2/20)
        probability = 0.7
        rand = random.sample(range(0, (genes*genes)-1), interactions)
        for aa in range(0, interactions):
            graph[rand[aa]//genes, rand[aa] %
                  genes] = np.random.choice([-1, 1], 1, p=[1-probability, probability])
    else:  # generate random network with N genes
        interactions = int(interactions)
        probability = 0.7
        rand = random.sample(range(0, (genes*genes)-1), interactions)
        for aa in range(0, interactions):
            graph[rand[aa]//genes, rand[aa] %
                  genes] = np.random.choice([-1, 1], 1, p=[1-probability, probability])
    #----------------------generate random initial parameters----------------------#
    init = []
    for aa in range(genes):
        rand = 10**random.uniform(-2.0, 2.0)
        init.append(rand)

    #----------------------generate random perturbations----------------------#
    unpe = np.random.normal(size=len(init))
    unpe = unpe/(scipy.linalg.norm(unpe))
    unpe = unpe/10

    #----------------------saving parameters of system for reusability----------------------#
    np.savetxt('params/g_{}.txt'.format(interactions), graph, fmt=('%i'))
    np.savetxt('params/init.txt', init, fmt=('%1.2f'))
    np.savetxt('params/perturb.txt', unpe, fmt='%1.5f')

    return graph, init, unpe


def generate_network(graph, init, hourtime, postfix):
    #----------------------------------input setup---------------------------------#
    df = init  # initial values
    timecost = 3600*hourtime  # time in seconds
    n = 500  # time step
    t = np.linspace(0, timecost, n)  # time vector for initial simulation
    #---------------------------------ODE computing--------------------------------#
    df = odeint(tree_model, df, t)
    np.savetxt('output/ft_{}.txt'.format(hourtime), df, fmt='%1.4f')

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
    np.savetxt('output/steady.txt', dfo, fmt='%1.4f')
    np.savetxt('output/perturbed.txt', dfp, fmt='%1.4f')

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


def expressions(opds, oprs):
    '''
    params:
    :opds: (float) the operands
    :oprs: (string) the operators, either 'and' or 'or'.
    output:
    :exp: (array) the postfix expression
    :display_exp: (array) the postfix expression, readable to write in an output .txt
    '''
    exp = opds[:2]
    opds = opds[2:]
    opd_group = 0
    while len(opds) > 0:
        p = random.random()
        if p < 0.5:
            exp.append(opds[0])
            opds = opds[1:]
            opd_group += 1
        else:
            for i in range(opd_group):
                if len(oprs) > 1:
                    exp.append(oprs[0])
                    oprs = oprs[1:]
            opd_group = 0
    while len(oprs) > 0:
        exp.append(oprs[0])
        oprs = oprs[1:]
    return exp


def generate_postfix(matrix):
    tree_matrix = matrix.copy()
    post_exps = []
    display_post_exps = []
    for i in range(len(matrix)):
        # if the matrix is full of 0s, append an empty array
        if tree_matrix[i].count(0) == len(matrix):
            post_exps.append([])
        else:
            # replace each value in the tree_matrix with an object that contains more information
            for j in range(len(tree_matrix[i])):
                if tree_matrix[i][j] != 0:
                    tree_matrix[i][j] = Token(tree_matrix[i][j], True, j)
            # remove all 0s from the matrix
            for x in range(tree_matrix[i].count(0)):
                tree_matrix[i].remove(0)
            # take the operands and mix them
            opds = tree_matrix[i]
            random.shuffle(opds)
            # generate a list of operators, ordered randomly
            if len(opds) == 1:
                oprs = np.random.choice(['and', 'or'], len(opds), replace=True)
            else:
                oprs = np.random.choice(
                    ['and', 'or'], len(opds)-1, replace=True)
            # generate postfix expressions
            exps = expressions(opds, oprs)
            post_exps.append(exps)
    # this .txt is EMPTY, need to figure out how to do index [and/or] index ...
    np.savetxt('output/post_exp.txt', display_post_exps, fmt='%s')
    pickle.dump(post_exps, open("output/post_exps.p", "wb"))
    return post_exps


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
                              str(columns) + " Genes", fontsize=18)
        scatterplotp.set_title(
            "Random Gene Regulatory Network with Perturbation", fontsize=18)
        if columns < 12:
            scatterplot.legend(loc=1)
            scatterplotp.legend(loc=1)
        scatterplot.set_xlim(0, (hourtime/4)/24)
        scatterplotp.set_xlim(0, (hourtime/4)/24)
        fig1.savefig(
            "output/distODE_after{}h.png".format(int(hourtime/2)), dpi=100)
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
        fig2.savefig("output/dist.png", dpi=100)
        np.savetxt('output/dist.txt', dist, fmt='%s')
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
                  str(columns) + " Genes")
        if columns < 12:
            plt.legend(loc=2)
        plt.xlim(0, hourtime)
        fig.savefig("output/ODE_{}h.png".format(hourtime))
        plt.close()


#----------------------------------token class for postfix---------------------------------#


class Token:
    '''
    this class stores information about a given operand:
    params:
    :self.value: (int) protein product concentration
    :self.isLogic: (bool) checks if value is a logic operand (to distinguish 0 or 1 between logic and protein concentration)
    :self.index: (int) index of gene

    '''

    def __init__(self, value, isLogic=False, index=None):
        self.value = value
        self.isLogic = isLogic
        self.index = index


if __name__ == '__main__':
    """
    the user's input will define the following:
    1. number of genes we want to simulate
    2. repressilator OR the number of interactions between the genes
    3. time of simulation in hours
    """
    parser = argparse.ArgumentParser(description='rgrn.py generates and simulates a random gene regulatory network given'
                                     'the number of genes in the network, number of interactions in or characteristic of the network, and the simulation time in hours.')
    parser.add_argument('-g', '--numberOfGenes', required=True, dest='genes',
                        help='number of genes in the network [positive int]')
    parser.add_argument('-i', '--numberOfInteractions', required=True, dest='interactions',
                        help='number of interactions in network or connected or repressilator [nonnegative int, \'C\', or \'R\']')
    parser.add_argument('-t', '--hourtime', required=True, dest='hourtime',
                        help='simulation time in hours [positive int]')
    args = parser.parse_args()
    genes = int(args.genes)
    interactions = args.interactions
    hourtime = int(args.hourtime)

    graph, init, perturbation = initialize(genes, interactions)
    post_exps = generate_postfix(graph.tolist())
    df = generate_network(graph, init, hourtime, post_exps)
    graph_network(df, hourtime)
    dfo, dfp = perturb(df, hourtime, perturbation)
    graph_network(dfo, hourtime, dfp)
