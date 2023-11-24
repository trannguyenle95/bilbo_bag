import numpy as np
import pandas as pd
import os
import matplotlib
from matplotlib import pyplot as plt

datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs')

#df = pd.read_csv(datafolder+'/flip/Opt_DMP/C/Easy/C_Opt_DMP_10l_bag_flip_Easy1.csv')


Bag = "C"

if Bag == "A":
    V_max = 6.0
    A_max = 220
elif Bag == "B":
    V_max = 7.5
    A_max = 360
elif Bag == "C":
    V_max = 12.5
    A_max = 370
elif Bag == "D":
    V_max = 14.0
    A_max = 530
elif Bag == "E":
    V_max = 22.0
    A_max = 550


# for i in range(1,11):
#     df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(i)+'.csv')
#     print("i: ", i)
#     print(df.A_alpha_rim)
#     plt.plot(df.A_alpha_rim, label='run'+str(i))

# plt.axhline(y = A_max*0.6, color = 'r', linestyle = '--', xmin = 0, xmax = 10) 
# plt.legend()
# plt.show()

#print(df)

# df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(1)+'.csv')
# bagC_A1 = df.A_alpha_rim.to_numpy() #.reshape(-1,1)

# print("bagA_C1: ", bagC_A1)

# df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(2)+'.csv')
# bagC_A2 = df.A_alpha_rim.to_numpy() #.reshape(-1,1)

# print("bagA_C2: ", bagC_A2)

# stacked = np.vstack((bagC_A1, bagC_A2))

# print("stacked: ", stacked)
# print("stacked shape: ", stacked.shape)

#TODO: need to pad arrays!
# df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(1)+'.csv')
# bagC_A = df.A_alpha_rim.to_numpy().reshape(-1,1)
# for i in range(2,11):
#      df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(i)+'.csv')
#      A_np = df.A_alpha_rim.to_numpy().reshape(-1,1)
#      bagC_A  = np.hstack((bagC_A, A_np))
#      print("iter: ", i)


# #TODO: create mean/std plots over dmp+difficulty, one plot for each bag (consistent baseline!) and each metric!
#TODO: DO PADDING FOR PLOT, BUT NOT WHEN CALCULATING VALUES TO TABLE!!

#NOTE: have length either 10 or 20 depending on set of runs!

# max_actions = 10

# A_metrics = np.zeros((max_actions, 10))
# df = pd.read_csv(datafolder+'/flip/Opt_DMP/'+Bag+'/Easy/C_Opt_DMP_10l_bag_flip_Easy'+str(1)+'.csv')
# #A_i = df.A_alpha_rim.to_numpy().reshape(-1,1)
# A_i = df.A_alpha_rim.to_numpy() / A_max
# print("A_i: ", A_i)

# #A_padded = np.pad(A_i, ((0,max_actions-len(A_i)),(0,0)), 'constant', constant_values=(A_i[-1]))
# A_padded = np.pad(A_i, ((0,max_actions-len(A_i))), 'constant', constant_values=(A_i[-1]))

# print("A padded: ", A_padded)

# A_metrics[:,1] = A_padded

# print("A metrics: \n", A_metrics)

# print("bagA_C: ", bagC_A)

def extract_arrays(mode, dmp, bag, difficulty):
    if mode == "flip":
        max_actions = 10
    else: #both "full" and "refine" have max 20 actions
        max_actions = 20

    datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data', 'runs', mode, dmp, bag, difficulty)
    if bag == "A":
        V_max = 6.0
        A_max = 220
    elif bag == "B":
        V_max = 7.5
        A_max = 360
    elif bag == "C":
        V_max = 12.5
        A_max = 370
    elif bag == "D":
        V_max = 14.0
        A_max = 530
    elif bag == "E":
        V_max = 22.0
        A_max = 500

    rows = max_actions + 1 #+1 because first row is for initial state

    Actions = np.empty((10)) #number of actions until successful termination, NaN if failed to reach goal in max actions
    A_metrics = np.empty((rows, 10))
    V_metrics = np.empty((rows, 10))
    E_metrics = np.empty((rows, 10))

    print("A lim: ", A_max*0.6)
    print("V lim: ", V_max*0.7)
    print("A max: ", A_max)

    for i in range(1,11):
        df = pd.read_csv(datafolder+'/'+bag+'_'+dmp+'_10l_bag_flip_'+difficulty+str(i)+'.csv')
        A = df.A_alpha_rim.to_numpy()
        V = df.Vol.to_numpy()
        E = df.E_rim.to_numpy()

        if ((A[-1] >= 0.6*A_max) and (V[-1] >= 0.7*V_max)): #NOTE: do not REQUIRE that elongation is within limits? Change this if needed for "full" mode!!
            Actions[i-1] = len(A)-1 #minus one because the first row is for initial state
        else:
            Actions[i-1] = np.NaN

        #pad with final values
        A_pad = np.pad(A, ((0,rows-len(A))), 'constant', constant_values=(A[-1]))
        V_pad = np.pad(V, ((0,rows-len(V))), 'constant', constant_values=(V[-1]))
        E_pad = np.pad(E, ((0,rows-len(E))), 'constant', constant_values=(E[-1]))

        A_metrics[:,i-1] = A_pad
        V_metrics[:,i-1] = V_pad
        E_metrics[:,i-1] = E_pad

    return A_metrics, V_metrics, E_metrics, Actions #NOTE: return as absolute values for easy analysis, can divide A_metrics and V_metrics by max value to convert to percentages later if needed!



Am, Vm, Em, Act = extract_arrays('flip', 'Opt_DMP', 'C', 'Hard')

#TODO: TODO: measure DISTANCE FROM 1.0 IN ELONGATION (as different runs could rarely overshoot 1.0)

#TODO: process "Actions" and get failure rate and avg actions to complete from it!!


np.set_printoptions(linewidth=np.inf)
print("A metrics:\n", Am) #TODO: study in CSV if result is correct!!

print("V metrics:\n", Vm) #TODO: study in CSV if result is correct!!

print("E metrics:\n", Em) #TODO: study in CSV if result is correct!!

print("Actions :\n", Act) #TODO: study in CSV if result is correct!!


# NOTE
#IF I WANT METRICS FROM FINAL VALUES THEN I CAN JUST CALCULATE FROM INDEX 10 IN PADDED ARRAY
#I STILL NEED WAY TO EXTRACT FAIL / PASS RATE AND AVG RUNS TO COMPLETE!