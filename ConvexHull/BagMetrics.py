import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
import os
import plotly.graph_objects as go
import plotly.express as px

#NOTE: maybe preprocess to remove outliers (distance over e.g. 0.5m from mean or median point)? Only needed if OptiTrack streams false positives

def BagMetrics(filename, max_area, width, plot = False):
    path = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', filename)
    points3d = np.loadtxt(path, delimiter=",")


    #filter out outliers
    x_med = np.median(points3d[:,0])
    z_med = np.median(points3d[:,2])
    i_outliers = []
    for i, p in enumerate(points3d):
        if (np.abs(p[0]-x_med) > width) or (np.abs(p[2]-z_med) > width):
            i_outliers.append(i)

    outlier_points = points3d[i_outliers]
    points3d = np.delete(points3d, i_outliers, axis = 0)


    #filter out points that are too close
    pairdist = cdist(points3d, points3d, 'euclidean')
    pairdist += np.eye(pairdist.shape[0])#add large value to diagonal, so that distance from point to itself is ignored

    #remove rows in points3d that correspond to duplicate points so that only one is perserved
    #NOTE: define a good threshold - bag opening is flexible so markers can move closer and further appart!
    threshold = 0.001 #markers within this distance (in meters) are treated as the same
    
    filter = np.where(pairdist < threshold, 1, 0)
    filter = np.tril(filter) #take lower triangular matrix to remove double references to same pairs
    i_duplicate = np.where(np.any(filter, axis = 1)) #indexes of points too close to some other point
    duplicate_points = points3d[i_duplicate]
    points3d = np.delete(points3d, i_duplicate, axis=0) 

    #require that bag is not tilted so that projection directly onto xz plane gives area.
    #this should correspond to using the area as viewed by a top-down camera, and it would be unnecesarily complicated to deal 
    #with tilted planes

    rim_point_mask = np.zeros(points3d.shape[0], dtype=bool)

    for i, point  in enumerate(points3d):
        if np.abs(point[1] - np.max(points3d[:,1])) < np.abs(point[1] - np.min(points3d[:,1])):
            rim_point_mask[i] = 1

    rim_points = points3d[rim_point_mask, :]

    points2d = rim_points[:, [0,2]]


    hull2d = ConvexHull(points2d)

    hull3d = ConvexHull(points3d)
    volume3d = hull3d.volume

    #"volume" for 2D hull is actual area, and "area" is actual perimeter
    area_ratio = hull2d.volume #/ max_area

    #Use PCA major and minor axes for elongation measure (like in the AutoBag paper by Chen et al. 2023 https://doi.org/10.48550/arXiv.2210.17217)
    #If bounding box was used instead there would be problems if the bag is e.g. slim but diagonal wrt. the coordinate axes like so: / , 
    # because then the bounding box would be much larger than the real area, and it also could have equally long components 
    pca = PCA(n_components=2)
    pca.fit(points2d[hull2d.vertices])
    pca_axes = pca.components_
    pca_magnitude = pca.explained_variance_
    elongation = np.sqrt(pca_magnitude[1]) / np.sqrt(pca_magnitude[0]) 
    #short over long axis in my implementation (maybe other in AutoBag) > this way elongation is 1 if opening is perfectly round and a smaller value for worse cases

    #Plotting
    if plot:
        #Print here so that metrics are visible without having to close pyplot figures first
        print("Convex Hull area: ", area_ratio)
        print("Convex Hull elongation: ", elongation)

        print("3d hull volume: ", volume3d)

        convex_hull_plot_2d(hull2d) #plot 2d convex hull
        ax1 = plt.gca()
        ax1.set_aspect('equal') #Set equal axes so that orthogonality of components can be verified
        data_mean = np.mean(points2d, axis=0) #get mean of data for plotting ellipoid axes at right point

        plt.plot((data_mean[0], data_mean[0]+pca_axes[0,0]*np.sqrt(pca_magnitude[0])), (data_mean[1], data_mean[1]+pca_axes[0,1]*np.sqrt(pca_magnitude[0]))) #Plot major axis - note sqrt() to get standard deviation from variance for plotting
        plt.plot((data_mean[0], data_mean[0]+pca_axes[1,0]*np.sqrt(pca_magnitude[1])), (data_mean[1], data_mean[1]+pca_axes[1,1]*np.sqrt(pca_magnitude[1]))) #Plot minor axis - note sqrt() to get standard deviation from variance for plotting

        #Also plot negative directions to get better idea of ellipsoid fit
        plt.plot((data_mean[0], data_mean[0]-pca_axes[0,0]*np.sqrt(pca_magnitude[0])), (data_mean[1], data_mean[1]-pca_axes[0,1]*np.sqrt(pca_magnitude[0]))) #Plot major axis - note sqrt() to get standard deviation from variance for plotting
        plt.plot((data_mean[0], data_mean[0]-pca_axes[1,0]*np.sqrt(pca_magnitude[1])), (data_mean[1], data_mean[1]-pca_axes[1,1]*np.sqrt(pca_magnitude[1]))) #Plot minor axis - note sqrt() to get standard deviation from variance for plotting
        
        plt.show()

        #plot median point
        ax2 = plt.gca()
        ax2.set_aspect('equal') #Set equal axes so that orthogonality of components can be verified

        plt.plot(x_med, z_med, 'rx')
        #show where outlier points were removed
        for p in outlier_points:
            plt.plot(p[0], p[2], 'ro')

        #show where duplicate points were removed
        for p in duplicate_points:
            plt.plot(p[0], p[2], 'o', color='tab:orange')

        #show remaining points
        for p in points3d:
            plt.plot(p[0], p[2], 'go')

        plt.show()

        ch_points = points3d[hull3d.vertices]
        bottom_points = points3d[np.logical_not(rim_point_mask)]

        print("b:", bottom_points)

        fig = go.Figure(data=[go.Scatter3d(x=rim_points[:, 0], y=rim_points[:, 2], z=rim_points[:, 1],
                                   mode='markers')]).update_traces(marker=dict(color='red'))
        fig.add_trace(go.Mesh3d(x=ch_points[:, 0], 
                         y=ch_points[:, 2], #change y and z so that z points upwards
                         z=ch_points[:, 1], #change y and z so that z points upwards
                         color="blue", 
                         opacity=0.2,
                         alphahull=0))
        
        fig.add_trace(go.Scatter3d(x=bottom_points[:, 0], y=bottom_points[:, 2], z=bottom_points[:, 1],
                                   mode='markers', marker=dict(
        #size=12,
        color="blue",                # set color to an array/list of desired values
        opacity=0.8
    )))

        fig.show()

    return area_ratio, elongation

if __name__ == '__main__':
    A_CH, E_CH = BagMetrics("tracked_markers.csv", 0.085, 0.4, plot=True)
    #print("Convex Hull area: ", A_CH)
    #print("Convex Hull elongation: ", E_CH)
