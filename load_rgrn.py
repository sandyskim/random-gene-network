import rgrn
import os
import sys
import numpy as np
from scipy.integrate import odeint
import pickle
#----------------------------------making directories---------------------------------#
os.makedirs('output', exist_ok=True)
os.makedirs('params', exist_ok=True)

#----------------------------------functions---------------------------------#


def load_network(filename, hourtime):
    filename = filename
    matrix = ((np.loadtxt("params/" + filename +
                          '.txt', delimiter=" ")).astype(int)).tolist()
    init = np.loadtxt('params/init.txt', delimiter=" ")
    unpe = np.loadtxt("params/perturb.txt", delimiter=" ")
    post_exps = pickle.load(open("output/post_exps.p", "rb"))
    return matrix, init, unpe, post_exps


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
    filename = sys.argv[1]
    hourtime = int(sys.argv[2])
    graph, init, perturbation, post_exps = load_network(filename, hourtime)
    df = generate_network(graph, init, hourtime, post_exps)
    rgrn.graph_network(df, hourtime)
