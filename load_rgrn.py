#-----------------------------------packages-----------------------------------#
from __future__ import division
import numpy as np
import scipy
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from scipy.optimize import brentq
from scipy.optimize import fsolve
import random
import sys
from scipy.integrate import odeint
from itertools import combinations, permutations, product, islice
from statistics import mean
import math
import warnings
import time

#----------------------------------input setup---------------------------------#
""" 
the user's input will define the following:
1. filename of matrix (without .txt)
2. time of simulation in hours 
"""

filename = sys.argv[1]
hourtime = int(sys.argv[2])
togetnumber = filename.split("_")
interactions = int(togetnumber[1])
matrix = ((np.loadtxt("params/" + filename + '.txt', delimiter=" ")).astype(int))
size = matrix.shape[0]
matrix = matrix.tolist()
df0 = np.loadtxt('params/init.txt', delimiter=" ")
timecost = 3600*hourtime
n = 500
t = np.linspace(0, timecost, n)
t2 = np.linspace(0, timecost/2, n)
kmax = 1.2 * 10**(-2)  # mM/s
kbas = 1.2 * 10**(-8)  # mM/s
kdec = 1.2 * 10**(-5)  # 1/s
alph = 2 # Hill coefficient
kprot = 1  # mM
df = df0
unpe = np.loadtxt("params/perturb.txt", delimiter = " ")

#-------------------------------define functions-------------------------------#
def tree_model(df, t):
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
                            a = (df[a.index]**alph)/(kprot**alph + df[a.index]**alph)
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
                            b = (df[b.index]**alph)/(kprot**alph + df[b.index]**alph)
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

tree_matrix = matrix.copy()
post_exps = []
display_post_exps = []

class Token:
    def __init__(self, value, isLogic = False, index = None):
        self.value = value
        self.isLogic = isLogic
        self.index = index
 
def expressions(opds, oprs):
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
    display_exp = []
    for i in range(len(exp)):
        if exp[i] == 'and' or exp[i] == 'or':
            display_exp.append(exp[i])
        else:
            display_exp.append(exp[i].value)
    return exp, display_exp

for i in range(len(matrix)):
    if tree_matrix[i].count(0) == size:
        post_exps.append([])
    else:
        idx = 0
        for j in range(len(tree_matrix[i])):
            if tree_matrix[i][j] != 0:
                tree_matrix[i][j] = Token(tree_matrix[i][j], True, idx)
            idx += 1
        for x in range(tree_matrix[i].count(0)):
            tree_matrix[i].remove(0)
        opds = tree_matrix[i]  
        random.shuffle(opds)
        if len(opds) == 1:
                oprs = np.random.choice(['and', 'or'], len(opds), replace=True)
        else:
            oprs = np.random.choice(['and', 'or'], len(opds)-1, replace=True)
        exps = expressions(opds, oprs)
        post_exps.append(exps[0])
        display_post_exps.append(exps[1])
        
#---------------------------------ODE computing--------------------------------#
df = odeint(tree_model, df, t)

#----------------------saving output for reusability----------------------#
os.makedirs('output', exist_ok=True)
#np.savetxt('output/timing.txt', timefile, fmt='%i %i %1.1f %1.1f %1.5f')
np.savetxt('output/ft_{}.txt'.format(hourtime), df, fmt='%1.4f')
np.savetxt('output/post_exp.txt', display_post_exps, fmt='%s')

#----------------------perturbing the system----------------------#
dfo = df[-1]
dfo = odeint(tree_model, dfo, t2)
dfp = df[-1] + unpe
dfp = odeint(tree_model, dfp, t2)

#----------------------saving output for reusability----------------------#
np.savetxt('output/steady.txt',df,fmt='%1.4f')
np.savetxt('output/perturbed.txt',dfp,fmt='%1.4f')

#------------------------------input from terminal-----------------------------#
t1 = np.linspace(0,hourtime,500)
t2 = np.linspace(0,(hourtime/2)/24,500)
columns = df.shape[1]

#------------------------graphing protein concentration-----------------------#
fig = plt.figure()
scatterplot = fig.add_subplot(111)
colors = ['r', 'y', 'g', 'c', 'b', 'm']
lines = ['--', ':', '-.']
graphs = []
for r in product(lines, colors): 
    graphs.append(r[1] + r[0])
for i in range(columns):
    scatterplot.plot(t1, df[:,i], graphs[i%18] ,label=str(chr(i+65))+"(t)")
plt.xlabel('time (h)', fontsize= 15)
plt.ylabel('protein product concentration (µM)', fontsize=15)
plt.title("Random Gene Regulatory Network of " + str(columns) + " Genes")
if columns < 12:
    plt.legend(loc=2)
plt.xlim(0,hourtime)
fig.savefig("output/ODE_{}h.png".format(hourtime))
plt.close()

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
    scatterplot.plot(t2, dfo[:,i], graphs[i%18], label=str(chr(i+65))+"(t)")
    scatterplotp.plot(t2, dfp[:,i], graphs[i%18], label=str(chr(i+65))+"(t)")
scatterplotp.set_xlabel('time (days)',fontsize=14)
scatterplot.set_ylabel('protein concentration (' + u"\u03bcM)")
scatterplotp.set_ylabel('protein concentration (' + u"\u03bcM)")
scatterplot.set_title("Random Gene Regulatory Network of " + str(columns) + " Genes",fontsize=18)
scatterplotp.set_title("Random Gene Regulatory Network with Perturbation",fontsize=18)
if columns < 12:
    scatterplot.legend(loc=1)
    scatterplotp.legend(loc=1)
scatterplot.set_xlim(0,(hourtime/2)/24)
scatterplotp.set_xlim(0,(hourtime/2)/24)
fig1.savefig("output/distODE_after{}h.png".format(int(hourtime/2)), dpi=100)
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
displot.set_title("Distance Between Perturbed and Steady States",fontsize=16)
displot.set_xlabel('time (days)',fontsize=14)
displot.set_ylabel('distance',fontsize=14)
fig2.savefig("output/dist.png", dpi=100)
plt.close()

#------------------------calculating largest lyapunov exponent----------------------#
difference = mean(lle_distances)
perturb = np.linalg.norm(unpe)
lle = 1/n * np.log(difference/perturb)
np.savetxt('output/lle.txt', [lle], fmt='%s')