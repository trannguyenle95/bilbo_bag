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


def plot_three_results(df, mode, difficulty, bags_number, legend_list,
                       target_elongation= True, flip_only_lines=False,
                       dpi=300):
    hspace_val = 0.22
    wspace_val = 0.22

    fig, axes = plt.subplots(nrows=3, ncols = 3,dpi=dpi)
    fig.set_size_inches(10, 10, forward=True)

    # fig.suptitle("Fling-only primitive", fontsize=26)
    axes[0,0].set_title("Area")
    axes[0,1].set_title("Volume")
    axes[0,2].set_title("Distance to elongation")

    if bags_number == 1:
        axes[0,0].set_ylabel(f"BagA {difficulty}", size='large')
        bags = ["A"]
    elif bags_number == 2:
        axes[0,0].set_ylabel(f"BagD {difficulty}", size='large')
        axes[1,0].set_ylabel(f"BagE {difficulty}", size='large')
        bags = ["D", "E"]
    else:
        axes[0,0].set_ylabel(f"BagA {difficulty}", size='large')
        axes[1,0].set_ylabel(f"BagB {difficulty}", size='large')
        axes[2,0].set_ylabel(f"BagC {difficulty}", size='large')
        bags = ["A", "B", "C"]

    axes_list = ["A", "V", "E"]
    for i in range(bags_number):
        for j in range(3):
            results(df, mode, bags[i], difficulty, axes_list[j], axes[i, j]) # axes[row, col]

    Ah, Al = axes[0,0].get_legend_handles_labels()
    Eh, El = axes[0,2].get_legend_handles_labels()
    if target_elongation:
        Ah.append(Eh[-1]) #add graphic for target elongation
    if flip_only_lines:
        Ah.append(Eh[1]) #add graphic for "flip-only" lines

    fig.legend(Ah, legend_list,  ncol = len(legend_list), bbox_to_anchor=(0.4, 0.05, 0.5, 0.95))
    fig.subplots_adjust(left=0.1, bottom=0.05, right=0.98, top=0.88, 
                        wspace=wspace_val, hspace=hspace_val)
    # fig.supxlabel('Actions', y = 0.06)
    bags_names = ''.join(bags)
    plt.savefig(f'results_{mode}_{difficulty}_{bags_names}.pdf', dpi=300)


#------------------ initialize dataframe ------------------
df = pd.DataFrame(columns = ['A_mean_final', 'V_mean_final', 'E_mean_dist_final', 'Success_rate', 'Mean_act_successful', 'Full_success_rate'])
# ------------------ flip-only runs ------------------

# ---------- Easy only -----------------
dpi=300
plot_three_results(df, mode='flip', difficulty='Easy', dpi=300, bags_number=3,
                   legend_list=['tau-DMP', 'TC-DMP', 'Opt-DMP', 'lower limit', 'upper limit'])

# ---------- Hard only -----------------
plot_three_results(df, mode='flip', difficulty='Hard', dpi=300, bags_number=3,
                   legend_list=['tau-DMP', 'TC-DMP', 'Opt-DMP', 'lower limit', 'upper limit'])


# ------------------ refinement only runs ------------------
plot_three_results(df, mode='refine', difficulty='Hard', dpi=300, bags_number=1,
                   legend_list=['refinement-only', 'lower limit', 'upper limit'])

# ------------------ full pipeline runs ------------------

#NOTE: ---------- Bags A-C only -----------------
plot_three_results(df, mode='full_pipeline', difficulty='Hard', dpi=300, bags_number=3,
                   legend_list=['fling-only Opt-DMP', 'full Opt-DMP', 'lower limit', 'upper limit'],
                   flip_only_lines=True)


#NOTE: ---------- Bags D&E only -----------------
plot_three_results(df, mode='full_pipeline', difficulty='Hard', dpi=300, bags_number=2,
                   legend_list=['Opt-DMP', 'lower limit', 'upper limit'])


# --------------------------------------------------------
# plt.show() #show all plots generated above
# np.set_printoptions(linewidth=np.inf)
# print("df:\n", df) #show table filled out above
# --------------------------------------------------------