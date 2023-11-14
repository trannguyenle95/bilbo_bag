import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
import os
import plotly.graph_objects as go
import plotly.express as px
import rospy
from sensor_msgs.msg import PointCloud
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import time
import alphashape


def calculate_metrics(bag, displayPlot = False):
    #Set widths and parameters for setting alpha.
    #Parameters k and b in "alpha = kx + b" are fitted from measured convex hull values "x" and manually selected appropriate alpha values.
    if bag == "A":
        width = 0.38
        k = -0.0443
        b = 24.281
    elif bag == "B":
        width = 0.41
        k = -0.0263
        b = 20.083
    elif bag == "C":
        width = 0.42
        k = -0.0226
        b = 18.917
    elif bag == "D":
        width = 0.49
        k = -0.0290
        b = 24.371
    elif bag == "E":
        width = 0.53
        k = -0.0131
        b = 15.092

    #Get points from NatNet rostopic
    cloud = rospy.wait_for_message("/natnet_ros/pointcloud", PointCloud, timeout=None)
    points3d = np.empty((len(cloud.points), 3)) #Note: y-axis denotes height, axes are only flipped when plotting the 3D plot

    for i, point in enumerate(cloud.points):
        #NOTE: this ix xyz order from the point cloud. As Motive exports with y-axis pointing upwards this is handled in BagMetrics()
        points3d[i][0] = point.x
        points3d[i][1] = point.y
        points3d[i][2] = point.z

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


    #get ponints from the 3D cloud belonging to the rim
    if np.max(points3d[:,1]) - np.min(points3d[:,1]) < 0.10: #if spread is less than 10cm then it is assumed only rim point are visible in initial state
        rim_point_mask = np.ones(points3d.shape[0], dtype=bool)
    else: #there are both rim and non-rim points
        rim_point_mask = np.zeros(points3d.shape[0], dtype=bool)
        for i, point  in enumerate(points3d):
            if np.abs(point[1] - np.max(points3d[:,1])) < np.abs(point[1] - np.min(points3d[:,1])):
                rim_point_mask[i] = 1
                
    rim_points = points3d[rim_point_mask, :]
    points2d = rim_points[:, [2,0]] #use x and z axis and reorder so plot is oriented same way as bag in lab (robot x is optitrack z, robot y is optitrack x)

    hull2d = ConvexHull(points2d)

    hull3d = ConvexHull(points3d)
    volume3d = hull3d.volume
    volume3d /= 0.001 #convert from m3 to liters
    rim_area_CH = hull2d.volume #Convex Hull area, volume" for 2D hull is actual area, and "area" is actual perimeter
    rim_area_CH /= 0.0001 #convert from m2 to cm2

    #area when approximating rim with potentially non-convex polygon using alpha shapes https://pypi.org/project/alphashape/
    alpha = k*rim_area_CH + b
    print("alpha value: ", alpha)
    alpha_shape = alphashape.alphashape(points2d, alpha)
    if isinstance(alpha_shape, Polygon):
        rim_area_poly = alpha_shape.area /  0.0001 #convert from m2 to cm2
    elif isinstance(alpha_shape, MultiPolygon):
        Polygons = list(alpha_shape.geoms)
        rim_area_poly = 0
        for poly in Polygons:
            rim_area_poly += poly.area /  0.0001 #convert from m2 to cm2


    #Use PCA major and minor axes for elongation measure (like in the AutoBag paper by Chen et al. 2023 https://doi.org/10.48550/arXiv.2210.17217)
    #If bounding box was used instead there would be problems if the bag is e.g. slim but diagonal wrt. the coordinate axes like so: / , 
    # because then the bounding box would be much larger than the real area, and it also could have equally long components 
    pca = PCA(n_components=2)
    pca.fit(points2d[hull2d.vertices])
    pca_axes = pca.components_

    pca_magnitude = pca.explained_variance_

    if abs(pca_axes[0][0]) > abs(pca_axes[1][0]): #major axis is mostly aligned with x-axis
        elongation = np.sqrt(pca_magnitude[1]) / np.sqrt(pca_magnitude[0]) 
    else: #major axis is mostly aligned with y-axis
        elongation = np.sqrt(pca_magnitude[0]) / np.sqrt(pca_magnitude[1])
    #elongation is now <1 if ellipse is too long along x axis, >1 if ellipse is too long along y axis, and 1 if ellipse is perfectly round


    #Plotting
    if displayPlot:
        #Print here so that metrics are visible without having to close pyplot figures first
        print("Convex Hull elongation: ", elongation)
        print("Rim area (cm2): ", rim_area_CH)
        print("3d hull volume (liters): ", volume3d)
        print("Non-convex rim area (cm2): ",  rim_area_poly)


        convex_hull_plot_2d(hull2d) #plot 2d convex hull
        ax1 = plt.gca()
        ax1.set_aspect('equal') #Set equal axes so that orthogonality of components can be verified
        data_mean = np.mean(points2d, axis=0) #get mean of data for plotting ellipoid axes at right point

        plt.plot((data_mean[0], data_mean[0]+pca_axes[0,0]*np.sqrt(pca_magnitude[0])), (data_mean[1], data_mean[1]+pca_axes[0,1]*np.sqrt(pca_magnitude[0]))) #Plot major axis - note sqrt() to get standard deviation from variance for plotting
        plt.plot((data_mean[0], data_mean[0]+pca_axes[1,0]*np.sqrt(pca_magnitude[1])), (data_mean[1], data_mean[1]+pca_axes[1,1]*np.sqrt(pca_magnitude[1]))) #Plot minor axis - note sqrt() to get standard deviation from variance for plotting

        #Also plot negative directions to get better idea of ellipsoid fit
        plt.plot((data_mean[0], data_mean[0]-pca_axes[0,0]*np.sqrt(pca_magnitude[0])), (data_mean[1], data_mean[1]-pca_axes[0,1]*np.sqrt(pca_magnitude[0]))) #Plot major axis - note sqrt() to get standard deviation from variance for plotting
        plt.plot((data_mean[0], data_mean[0]-pca_axes[1,0]*np.sqrt(pca_magnitude[1])), (data_mean[1], data_mean[1]-pca_axes[1,1]*np.sqrt(pca_magnitude[1]))) #Plot minor axis - note sqrt() to get standard deviation from variance for plotting
        
        #plot potentially non-convex polygon too
        if isinstance(alpha_shape, Polygon):
            plt.plot(*alpha_shape.exterior.xy, '--r')
        elif isinstance(alpha_shape, MultiPolygon):
            Polygons = list(alpha_shape.geoms)
            for poly in Polygons:
                plt.plot(*poly.exterior.xy, '--r')
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

        #plt.show()

        ch_points = points3d[hull3d.vertices]
        bottom_points = points3d[np.logical_not(rim_point_mask)]


        #In robot coordinte system x is z from optitrack, y is x from optitrack, and z is y from optitrack 

        fig = go.Figure(data=[go.Scatter3d(x=rim_points[:, 2],
                                y=rim_points[:, 0],
                                z=rim_points[:, 1],
                                mode='markers')]).update_traces(marker=dict(color='red'))
        fig.add_trace(go.Mesh3d(x=ch_points[:, 2], 
                         y=ch_points[:, 0],
                         z=ch_points[:, 1],
                         color="blue", 
                         opacity=0.2,
                         alphahull=0))
        
        fig.add_trace(go.Scatter3d(x=bottom_points[:, 2], y=bottom_points[:, 0], z=bottom_points[:, 1],
                                   mode='markers', marker=dict(
        #size=12,
        color="blue",                # set color to an array/list of desired values
        opacity=0.8
        )))

        fig.update_layout(scene_aspectmode='data')

        fig.show()

    return rim_area_CH, rim_area_poly, volume3d, elongation

if __name__ == '__main__':
    rospy.init_node('listener', anonymous=True)
    time.sleep(5.0) #wait 5s for time to adjust bag
    A_CH_rim, A_poly_rim, Vol, E_rim = calculate_metrics('E', displayPlot=True)
