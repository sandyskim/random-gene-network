import numpy as np
import pandas as pd
import os
import argparse
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics.cluster import normalized_mutual_info_score
import sys

#----------------------------functions----------------------------#


def run_pca(pca_dim, hourtime):
    #----------------------------embed data using PCA----------------------------#
    data = np.loadtxt("output/ft_{}.txt".format(hourtime), delimiter=" ")
    data = pd.DataFrame(data)
    data = data[0:500]

    sc = StandardScaler()
    x = sc.fit_transform(data)

    if pca_dim == 2:
        pca = PCA(n_components=2)
        pca.fit(x)
        x_pca = pca.fit_transform(x)
        principal_data = pd.DataFrame(data=x_pca, columns=['PC 1', 'PC 2'])

    if pca_dim == 3:
        pca = PCA(n_components=3)
        pca.fit(x)
        x_pca = pca.fit_transform(x)
        principal_data = pd.DataFrame(
            data=x_pca, columns=['PC 1', 'PC 2', 'PC 3'])
    np.savetxt("output/principal_data_pca{}dim.txt".format(pca_dim),
               principal_data, fmt='%1.4f')

    fig = plt.figure(figsize=(8, 8))
    if pca_dim == 2:
        ax = fig.add_subplot(1, 1, 1)
    if pca_dim == 3:
        ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('PC 1', fontsize=15)
    ax.set_ylabel('PC 2', fontsize=15)
    if pca_dim == 2:
        ax.set_title('2 component PCA', fontsize=20)
        ax.scatter(x_pca[:, 0], x_pca[:, 1], alpha=0.2)
    if pca_dim == 3:
        ax.set_zlabel('PC 3', fontsize=15)
        ax.set_title('3 component PCA', fontsize=20)
        ax.scatter(x_pca[:, 0], x_pca[:, 1], x_pca[:, 2], alpha=0.2)
    fig.savefig("output/pca_{}dim.png".format(pca_dim))
    plt.close()

    return principal_data


def takens_embedding(data, delay, dimension):
    #----------------------------takens' embedding----------------------------#
    if delay*dimension > len(data):
        raise NameError('delay times dimension exceed length of data')
    embeddedData = np.array([data[0:len(data)-delay*dimension]])
    for i in range(1, dimension):
        embeddedData = np.append(
            embeddedData, [data[i*delay:len(data) - delay*(dimension - i)]], axis=0)
    return embeddedData


if __name__ == '__main__':

    """ 
    the user's input will define the following:
    1. hourtime of the simulation we want to embed
    2. PC we want to apply taken's
    3. delay in seconds
    """
    parser = argparse.ArgumentParser(description='takens.py takes high-dimensional time series data generated and reduces the dimensions using PCA'
                                     'and embeds a given principal components using the Takens\'s embedding given the time of the simulation in hours,'
                                     'the principal component to be embedded, and the time delay for the embedding.')
    parser.add_argument('-t', '--hourtime', required=True, dest='hourtime',
                        help='hourtime of the simulation we want to embed [positive int]')
    parser.add_argument('-p', '--principalComponent', required=True, dest='principalcomponent',
                        help='principal component we want to apply the embedding [1, 2, or 3]')
    parser.add_argument('-d', '--delay', required=True, dest='delay',
                        help='delay [positive int]')
    args = parser.parse_args()
    hourtime = int(args.hourtime)
    PC = int(args.principalcomponent)
    delay = int(args.delay)

    run_pca(2, hourtime)
    data = run_pca(3, hourtime)
    n = len(data)
    data = data['PC {}'.format(PC)]
    t = np.linspace(0, hourtime, n)

    embedded_data = takens_embedding(data, delay, PC-1)
    np.savetxt("output/takens_data.txt", embedded_data, fmt='%1.4f')

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(15, 14))
    ax[0].plot(t, data)
    ax[0].set_title('Time Series with delay = {}'.format(delay))
    ax[1].plot(embedded_data[0, :], embedded_data[1, :])
    ax[1].set_title('Takens Embedding Time Series')
    plt.savefig('output/takens2d.png')
    plt.clf()
    plt.close()

    fig = plt.figure(figsize=(8, 8))
    embedded_data = takens_embedding(data, delay, 3)
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.plot(embedded_data[0, :], embedded_data[1, :], embedded_data[2, :])
    ax1.set_title('3D Takens Embedding on 3D PCA')
    plt.savefig('output/takens3d.png')
    plt.close()
