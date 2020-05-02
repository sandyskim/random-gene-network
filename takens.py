import numpy as np 
import pandas as pd
import os
import math
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sys

hourtime = int(sys.argv[1])
PC = sys.argv[2]

#----------------------------embed data using PCA----------------------------#
def run_pca(pca_dim, hourtime):
    data = np.loadtxt("output/ft_{}.txt".format(hourtime), delimiter=" ")
    data = pd.DataFrame(data)
    data = data[0:500]
    
    sc = StandardScaler()
    x = sc.fit_transform(data)

    if pca_dim == 2:
        pca = PCA(n_components=2)
        pca.fit(x)
        x_pca = pca.fit_transform(x)
        principal_data = pd.DataFrame(data = x_pca, columns = ['PC 1', 'PC 2'])

    if pca_dim == 3:
        pca = PCA(n_components=3)
        pca.fit(x)
        x_pca = pca.fit_transform(x)
        principal_data = pd.DataFrame(data = x_pca, columns = ['PC 1', 'PC 2', 'PC 3'])
    np.savetxt("output/principal_data_pca{}dim.txt".format(pca_dim), principal_data, fmt='%1.4f')

    fig = plt.figure(figsize = (8,8))
    if pca_dim == 2:
        ax = fig.add_subplot(1,1,1) 
    if pca_dim == 3:
        ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('PC 1', fontsize = 15)
    ax.set_ylabel('PC 2', fontsize = 15)
    if pca_dim == 2:
        ax.set_title('2 component PCA', fontsize = 20)
        ax.scatter(x_pca[:,0], x_pca[:,1], alpha = 0.2)
    if pca_dim == 3:
        ax.set_zlabel('PC 3', fontsize = 15)
        ax.set_title('3 component PCA', fontsize = 20)
        ax.scatter(x_pca[:,0], x_pca[:,1], x_pca[:,2], alpha=0.2)
    fig.savefig("output/pca_{}dim.png".format(pca_dim))
    plt.close()

    return principal_data

data = run_pca(3, hourtime)

#----------------------------takens' embedding----------------------------#
def takens_embedding (data, delay, dimension):
    if delay*dimension > len(data):
        raise NameError('delay times dimension exceed length of data')    
    embeddedData = np.array([data[0:len(data)-delay*dimension]])
    for i in range(1, dimension):
        embeddedData = np.append(embeddedData, [data[i*delay:len(data) - delay*(dimension - i)]], axis=0)
    return embeddedData

data = data['PC {}'.format(PC)]
delay = 1
t = np.linspace(0, hourtime, 500)

embedded_data = takens_embedding(data, delay, 2)
np.savetxt("output/takens_data.txt", embedded_data, fmt='%1.4f')

fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(15,14))
ax[0].plot(t, data)
ax[0].set_title('Time Series with delay = {}'.format(delay))
ax[1].plot(embedded_data[0,:],embedded_data[1,:])
ax[1].set_title('Takens Embedding Time Series')
plt.savefig('output/takens2d.png')
plt.clf()
plt.close()

fig = plt.figure(figsize = (8,8))
embedded_data = takens_embedding(data, delay, 3)
ax1 = fig.add_subplot(111, projection='3d')
ax1.plot(embedded_data[0,:],embedded_data[1,:],embedded_data[2,:])
ax1.set_title('3D Takens Embedding on 3D PCA')
plt.savefig('output/takens3d.png')
plt.close()


'''
# cumulative variance

plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance')
plt.show()

# PCA information
print('Original number of features:', x.shape[1], file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
print('Reduced number of features:', x_pca.shape[1], file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
print('Variance attributed to principal component 1:', pca.explained_variance_ratio_[0], file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
print('Variance attributed to principal component 2:', pca.explained_variance_ratio_[1], file=open("output//pca_info_{}dim.txt".format(pca_dim), "a"))
if pca_dim == 3:
print('Variance attributed to principal component 3:', pca.explained_variance_ratio_[2], file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
print(pca.components_)
cov_mat = np.cov(x_pca.T)
eigen_vals, eigen_vecs = np.linalg.eig(cov_mat)
print('Covariance matrix:', cov_mat, file=open("{}/pca_info_{}dim.txt".format(i, pca_dim), "a"))
print('Eigenvalues of covariance matrix:', eigen_vals, file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
print('Eigenvectors of covariance matrix:', eigen_vecs, file=open("output/pca_info_{}dim.txt".format(pca_dim), "a"))
'''