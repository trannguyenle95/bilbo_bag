import numpy as np
import pandas as pd
import sys, os
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
from pathlib import Path


def extract_arrays(mode, dmp, bag, difficulty):
    src_dir = Path(os.path.abspath(os.path.dirname(sys.argv[0]))).parent
    if mode == "flip":
        max_actions = 10
    else: #both "full_pipeline" and "refine" have max 20 actions
        max_actions = 20
    if mode == "flip" or mode == "full_pipeline":
        datafolder = os.path.join(src_dir, 'Data', 'runs', mode, dmp, bag, difficulty)
    else: #refine has other structure
        datafolder = os.path.join(src_dir, 'Data', 'runs', mode)
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
    Num_actions = np.empty((10)) #number of actions until successful termination, NaN if failed to reach goal in max actions
    Full_success = np.empty((10))
    A_metrics = np.empty((rows, 10))
    V_metrics = np.empty((rows, 10))
    E_metrics = np.empty((rows, 10))

    flip_over_outliers = [] #runs that are terminated early because of flip-overs
    for i in range(1,11):
        if mode == "flip" or mode == "full_pipeline":
            df = pd.read_csv(datafolder+'/'+bag+'_'+dmp+'_10l_bag_flip_'+difficulty+str(i)+'.csv')
        else: #refine has other structure
            df = pd.read_csv(datafolder+'/'+bag+'_'+difficulty+str(i)+'.csv')


        A = df.A_alpha_rim.to_numpy()
        V = df.Vol.to_numpy()
        E = df.E_rim.to_numpy()
        Actions = df.Action.to_numpy()

        if ((A[-1] >= 0.6*A_max) and (V[-1] >= 0.7*V_max)):
            Num_actions[i-1] = len(A)-1 #minus one because the first row is for initial state
            if(abs(1-E[-1]) <= 0.2):
                Full_success[i-1] = 1
            else:
                Full_success[i-1] = 0
        else:
            Num_actions[i-1] = np.NaN
            if( (A.shape[0] < rows) and (Actions[-1] == 'F')):
                flip_over_outliers.append(i-1)
            Full_success[i-1] = 0

        #pad with final values
        A_pad = np.pad(A, ((0,rows-len(A))), 'constant', constant_values=(A[-1]))
        V_pad = np.pad(V, ((0,rows-len(V))), 'constant', constant_values=(V[-1]))
        E_pad = np.pad(E, ((0,rows-len(E))), 'constant', constant_values=(E[-1]))

        A_metrics[:,i-1] = A_pad
        V_metrics[:,i-1] = V_pad
        E_metrics[:,i-1] = E_pad

        failure_rate = np.isnan(Num_actions).sum() / 10
        if failure_rate < 1.0:
            avg_actions_in_successful = np.nanmean(Num_actions)
        else:
            avg_actions_in_successful = np.NaN #no successful runs

    
        full_success_rate = np.mean(Full_success)

    A_metrics = np.delete(A_metrics, flip_over_outliers, axis=1) #delete flip-over outliers from plot
    V_metrics = np.delete(V_metrics, flip_over_outliers, axis=1) #delete flip-over outliers from plot
    E_metrics = np.delete(E_metrics, flip_over_outliers, axis=1) #delete flip-over outliers from plot

    return A_metrics/A_max, V_metrics/V_max, abs(1-E_metrics), failure_rate, avg_actions_in_successful, full_success_rate


def plot_line(data, used_color, used_label, max_actions, ax):
    ax.plot(np.mean(data, axis = 1), color=used_color, label=used_label)
    ax.fill_between(np.linspace(0,max_actions, max_actions+1), np.mean(data, axis = 1)-np.std(data, axis = 1), np.mean(data, axis = 1)+np.std(data, axis = 1), color=used_color, alpha=0.3)

def write_table(A, V, E, fail_rate, act_successful, full_success_rate, mode, bag, difficulty, dmp, df):
    A_mean = np.mean(A, axis = 1)
    V_mean = np.mean(V, axis = 1)
    E_dist_mean = np.mean(E, axis = 1)
    df.loc[mode+'_'+bag+'_'+dmp+'_'+difficulty] = [A_mean[-1], V_mean[-1], E_dist_mean[-1], 1-fail_rate, act_successful, full_success_rate] 

def results(df, mode, bag, difficulty, metric, ax):
    if mode == 'flip':
        max_actions = 10
        A_tau_DMP, V_tau_DMP, E_tau_DMP, fail_tau_DMP, act_successful_tau_DMP, f_success_tau_DMP = extract_arrays(mode, 'tau_DMP', bag, difficulty)
        write_table(A_tau_DMP, V_tau_DMP, E_tau_DMP, fail_tau_DMP, act_successful_tau_DMP, f_success_tau_DMP, mode, bag, difficulty, 'tau_DMP', df)
        A_TC_DMP, V_TC_DMP, E_TC_DMP, fail_TC_DMP, act_successful_TC_DMP, f_success_TC_DMP = extract_arrays(mode, 'TC_DMP', bag, difficulty)
        write_table(A_TC_DMP, V_TC_DMP, E_TC_DMP, fail_TC_DMP, act_successful_TC_DMP, f_success_TC_DMP, mode, bag, difficulty, 'TC_DMP', df)
        A_Opt_DMP, V_Opt_DMP, E_Opt_DMP, fail_Opt_DMP, act_successful_Opt_DMP, f_success_Opt_DMP = extract_arrays(mode, 'Opt_DMP', bag, difficulty)
        write_table(A_Opt_DMP, V_Opt_DMP, E_Opt_DMP, fail_Opt_DMP, act_successful_Opt_DMP, f_success_Opt_DMP, mode, bag, difficulty, 'Opt_DMP', df)
    elif mode == 'refine':
            max_actions = 20
            A_refine, V_refine, E_refine, fail_refine, act_successful_refine, f_success_refine = extract_arrays(mode, 'none', bag, difficulty)
            write_table(A_refine, V_refine, E_refine, fail_refine, act_successful_refine, f_success_refine, mode, bag, difficulty, '', df)
    elif mode == 'full_pipeline':
            max_actions = 20
            if not(bag == 'D' or bag =='E'):
                flip_A_Opt_DMP, flip_V_Opt_DMP, flip_E_Opt_DMP, flip_fail_Opt_DMP, flip_act_successful_Opt_DMP, flip_f_success_tau_Opt_DMP = extract_arrays('flip', 'Opt_DMP', bag, difficulty) #for plotting comparison
            A_Opt_DMP, V_Opt_DMP, E_Opt_DMP, fail_Opt_DMP, act_successful_Opt_DMP, f_success_Opt_DMP = extract_arrays(mode, 'Opt_DMP', bag, difficulty)
            write_table(A_Opt_DMP, V_Opt_DMP, E_Opt_DMP, fail_Opt_DMP, act_successful_Opt_DMP, f_success_Opt_DMP, mode, bag, difficulty, 'Opt_DMP', df)
    if metric == 'A':
        if mode == 'flip':
            plot_line(A_tau_DMP, 'tab:blue', 'Area tau-DMP', max_actions, ax)
            plot_line(A_TC_DMP, 'tab:green', 'Area TC-DMP', max_actions, ax)
            plot_line(A_Opt_DMP, 'tab:orange', 'Area Opt-DMP', max_actions, ax)
        elif mode == 'refine':
            plot_line(A_refine, 'tab:orange', 'Area refine', max_actions, ax)
        elif mode == 'full_pipeline':
            if not(bag == 'D' or bag =='E'):
                plot_line(flip_A_Opt_DMP, 'tab:blue', 'Flip-only area Opt-DMP', 10, ax)
            plot_line(A_Opt_DMP, 'tab:orange', 'Area Opt-DMP', max_actions, ax)
        ax.axhline(y = 0.6, color = 'r', linestyle = '--', label='Area threshold') 
        ax.autoscale(enable=True, axis='x', tight=True)
        ax.set_ylim([0, 1.2])
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    elif metric == 'V':
        if mode == 'flip':
            plot_line(V_tau_DMP, 'tab:blue', 'Volume tau-DMP', max_actions, ax)
            plot_line(V_TC_DMP, 'tab:green', 'Volume TC-DMP', max_actions, ax)
            plot_line(V_Opt_DMP, 'tab:orange', 'Volume Opt-DMP', max_actions, ax)
        elif mode == 'refine':
            plot_line(V_refine, 'tab:orange', 'Volume refine', max_actions, ax)
        elif mode == 'full_pipeline':
            if not(bag == 'D' or bag =='E'):
                plot_line(flip_V_Opt_DMP, 'tab:blue', 'Flip-only volume Opt-DMP', 10, ax)
            plot_line(V_Opt_DMP, 'tab:orange', 'Volume Opt-DMP', max_actions, ax)
        ax.axhline(y = 0.7, color = 'r', linestyle = '--', label='Volume threshold') 
        ax.autoscale(enable=True, axis='x', tight=True)
        ax.set_ylim([0, 1.2])
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    elif metric == 'E':
        if mode == 'flip':
            plot_line(E_tau_DMP, 'tab:blue', 'Elongation tau-DMP', max_actions, ax)
            plot_line(E_TC_DMP, 'tab:green', 'Elongation TC-DMP', max_actions, ax)
            plot_line(E_Opt_DMP, 'tab:orange', 'Elongation Opt-DMP', max_actions, ax)
        elif mode == 'refine':
            plot_line(E_refine, 'tab:orange', 'Elongation refine', max_actions, ax)
        elif mode == 'full_pipeline':
            if not(bag == 'D' or bag =='E'):
                plot_line(flip_E_Opt_DMP, 'tab:blue', 'Flip-only elongation Opt-DMP', 10, ax)
            plot_line(E_Opt_DMP, 'tab:orange', 'Elongation Opt-DMP', max_actions, ax)
        ax.axhline(y = 0.2, color = 'b', linestyle = '--', label='Target elongation') 
        ax.autoscale(enable=True, axis='x', tight=True)
        ax.set_ylim([0, 1.2])
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    return


#------------------ initialize dataframe ------------------
df = pd.DataFrame(columns = ['A_mean_final', 'V_mean_final', 'E_mean_dist_final', 'Success_rate', 'Mean_act_successful', 'Full_success_rate'])
hspace_val = 0.25
# ------------------ flip-only runs ------------------

# ---------- Easy only -----------------
fig, axes = plt.subplots(nrows=3, ncols = 3)
fig.suptitle("Fling-only primitive", fontsize=26)
axes[0,0].set_title("Area")
axes[0,1].set_title("Volume")
axes[0,2].set_title("Distance to elongation 1.0")
axes[0,0].set_ylabel("BagA Easy", size='large')
axes[1,0].set_ylabel("BagB Easy", size='large')
axes[2,0].set_ylabel("BagC Easy", size='large')

results(df, 'flip', 'A', 'Easy', 'A', axes[0, 0]) # axes[row, col]
results(df, 'flip', 'A', 'Easy', 'V', axes[0, 1]) # axes[row, col]
results(df, 'flip', 'A', 'Easy', 'E', axes[0, 2]) # axes[row, col]

results(df, 'flip', 'B', 'Easy', 'A', axes[1, 0]) # axes[row, col]
results(df, 'flip', 'B', 'Easy', 'V', axes[1, 1]) # axes[row, col]
results(df, 'flip', 'B', 'Easy', 'E', axes[1, 2]) # axes[row, col]

results(df, 'flip', 'C', 'Easy', 'A', axes[2, 0]) # axes[row, col]
results(df, 'flip', 'C', 'Easy', 'V', axes[2, 1]) # axes[row, col]
results(df, 'flip', 'C', 'Easy', 'E', axes[2, 2]) # axes[row, col]

Ah, Al = axes[0,0].get_legend_handles_labels()
Eh, El = axes[0,2].get_legend_handles_labels()
Ah.append(Eh[-1]) #add graphic for target elongation
fig.legend(Ah, ['tau-DMP', 'TC-DMP', 'Opt-DMP', 'lower limit', 'upper limit'],  ncol = 5, loc='upper center',  bbox_to_anchor=(0.0, -0.06, 1.0, 1.0))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace_val)
fig.supxlabel('Actions', y = 0.06)

# ---------- Hard only -----------------
fig, axes = plt.subplots(nrows=3, ncols = 3)
fig.suptitle("Fling-only primitive", fontsize=26)
axes[0,0].set_title("Area")
axes[0,1].set_title("Volume")
axes[0,2].set_title("Distance to elongation 1.0")
axes[0,0].set_ylabel("BagA Hard", size='large')
axes[1,0].set_ylabel("BagB Hard", size='large')
axes[2,0].set_ylabel("BagC Hard", size='large')

results(df, 'flip', 'A', 'Hard', 'A', axes[0, 0]) # axes[row, col]
results(df, 'flip', 'A', 'Hard', 'V', axes[0, 1]) # axes[row, col]
results(df, 'flip', 'A', 'Hard', 'E', axes[0, 2]) # axes[row, col]

results(df, 'flip', 'B', 'Hard', 'A', axes[1, 0]) # axes[row, col]
results(df, 'flip', 'B', 'Hard', 'V', axes[1, 1]) # axes[row, col]
results(df, 'flip', 'B', 'Hard', 'E', axes[1, 2]) # axes[row, col]

results(df, 'flip', 'C', 'Hard', 'A', axes[2, 0]) # axes[row, col]
results(df, 'flip', 'C', 'Hard', 'V', axes[2, 1]) # axes[row, col]
results(df, 'flip', 'C', 'Hard', 'E', axes[2, 2]) # axes[row, col]

Ah, Al = axes[0,0].get_legend_handles_labels()
Eh, El = axes[0,2].get_legend_handles_labels()
Ah.append(Eh[-1]) #add graphic for target elongation
fig.legend(Ah, ['tau-DMP', 'TC-DMP', 'Opt-DMP', 'lower limit', 'upper limit'],  ncol = 5, loc='upper center',  bbox_to_anchor=(0.0, -0.06, 1.0, 1.0))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace_val)
fig.supxlabel('Actions', y = 0.06)



# ------------------ refinement only runs ------------------

fig, axes = plt.subplots(nrows=3, ncols = 3) #NOTE: Is 3 rows height sufficient, or should plot be higher/lower?? (other tests have 6 and 5 rows)
fig.suptitle("Refinement-only primitive", fontsize=26)
axes[0,0].set_title("Area")
axes[0,1].set_title("Volume")
axes[0,2].set_title("Distance to elongation 1.0")
axes[0,0].set_ylabel("BagA Hard", size='large')

results(df, 'refine', 'A', 'Hard', 'A', axes[0, 0]) # axes[row, col]
results(df, 'refine', 'A', 'Hard', 'V', axes[0, 1]) # axes[row, col]
results(df, 'refine', 'A', 'Hard', 'E', axes[0, 2]) # axes[row, col]

#deactivate rest of axes for shorter plot - cut figure to correct size when adding to thesis
axes[1, 0].axis('off')
axes[1, 1].axis('off')
axes[1, 2].axis('off')
axes[2, 0].axis('off')
axes[2, 1].axis('off')
axes[2, 2].axis('off')

Ah, Al = axes[0,0].get_legend_handles_labels()
Eh, El = axes[0,2].get_legend_handles_labels()
Ah.append(Eh[-1]) #add graphic for target elongation
fig.legend(Ah, ['refinement-only', 'lower limit', 'upper limit'],  ncol = 5, loc='upper center',  bbox_to_anchor=(0.0, -0.06, 1.0, 1.0))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace_val)
fig.supxlabel('Actions', y = 0.61)

# ------------------ full pipeline runs ------------------

#NOTE: ---------- Bags A-C only -----------------

fig, axes = plt.subplots(nrows=3, ncols = 3)
fig.suptitle("Full pipeline", fontsize=26)
axes[0,0].set_title("Area")
axes[0,1].set_title("Volume")
axes[0,2].set_title("Distance to elongation 1.0")
axes[0,0].set_ylabel("BagA Hard", size='large')
axes[1,0].set_ylabel("BagB Hard", size='large')
axes[2,0].set_ylabel("BagC Hard", size='large')

results(df, 'full_pipeline', 'A', 'Hard', 'A', axes[0, 0]) # axes[row, col]
results(df, 'full_pipeline', 'A', 'Hard', 'V', axes[0, 1]) # axes[row, col]
results(df, 'full_pipeline', 'A', 'Hard', 'E', axes[0, 2]) # axes[row, col]

results(df, 'full_pipeline', 'B', 'Hard', 'A', axes[1, 0]) # axes[row, col]
results(df, 'full_pipeline', 'B', 'Hard', 'V', axes[1, 1]) # axes[row, col]
results(df, 'full_pipeline', 'B', 'Hard', 'E', axes[1, 2]) # axes[row, col]

results(df, 'full_pipeline', 'C', 'Hard', 'A', axes[2, 0]) # axes[row, col]
results(df, 'full_pipeline', 'C', 'Hard', 'V', axes[2, 1]) # axes[row, col]
results(df, 'full_pipeline', 'C', 'Hard', 'E', axes[2, 2]) # axes[row, col]

Ah, Al = axes[0,0].get_legend_handles_labels()
Eh, El = axes[0,2].get_legend_handles_labels()
Ah.append(Eh[-1]) #add graphic for target elongation
Ah.append(Eh[1]) #add graphic for "flip-only" lines
fig.legend(Ah, ['fling-only Opt-DMP', 'full Opt-DMP', 'lower limit', 'upper limit'],  ncol = 5, loc='upper center',  bbox_to_anchor=(0.0, -0.06, 1.0, 1.0))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace_val)
fig.supxlabel('Actions', y = 0.06)


#NOTE: ---------- Bags D&E only -----------------

fig, axes = plt.subplots(nrows=4, ncols = 3)
fig.suptitle("Full pipeline", fontsize=26)
axes[0,0].set_title("Area")
axes[0,1].set_title("Volume")
axes[0,2].set_title("Distance to elongation 1.0")
axes[0,0].set_ylabel("BagD Hard", size='large')
axes[1,0].set_ylabel("BagE Hard", size='large')

results(df, 'full_pipeline', 'D', 'Hard', 'A', axes[0, 0]) # axes[row, col]
results(df, 'full_pipeline', 'D', 'Hard', 'V', axes[0, 1]) # axes[row, col]
results(df, 'full_pipeline', 'D', 'Hard', 'E', axes[0, 2]) # axes[row, col]

results(df, 'full_pipeline', 'E', 'Hard', 'A', axes[1, 0]) # axes[row, col]
results(df, 'full_pipeline', 'E', 'Hard', 'V', axes[1, 1]) # axes[row, col]
results(df, 'full_pipeline', 'E', 'Hard', 'E', axes[1, 2]) # axes[row, col]


#deactivate rest of axes for shorter plot - cut figure to correct size when adding to thesis
axes[2, 0].axis('off')
axes[2, 1].axis('off')
axes[2, 2].axis('off')
axes[3, 0].axis('off')
axes[3, 1].axis('off')
axes[3, 2].axis('off')

Ah, Al = axes[0,0].get_legend_handles_labels()
Eh, El = axes[0,2].get_legend_handles_labels()
Ah.append(Eh[-1]) #add graphic for target elongation
fig.legend(Ah, ['Opt-DMP', 'lower limit', 'upper limit'],  ncol = 5, loc='upper center',  bbox_to_anchor=(0.0, -0.06, 1.0, 1.0))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=hspace_val)
fig.supxlabel('Actions', y = 0.465)


# --------------------------------------------------------
plt.show() #show all plots generated above
np.set_printoptions(linewidth=np.inf)
print("df:\n", df) #show table filled out above
# --------------------------------------------------------