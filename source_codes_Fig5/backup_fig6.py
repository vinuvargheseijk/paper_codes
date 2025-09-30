import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
sns.set_context("poster")
#import pylustrator
#pylustrator.start()

directory = "./"

kcat_list = [20.0, 50.0, 200.0]
xticks = np.arange(0, len(kcat_list), 1)


def pValue_cat(p_value, significance_level = 0.05):
    if p_value > significance_level:
         text = 'n.s'
    if p_value < 0.05 and p_value > 0.01:
         text = '*'
    if p_value < 0.01 and p_value > 0.001:
         text = '* *'
    if p_value < 0.001:
         text = '* * *'
    return text

def sig_bars(ax, p_value, x1, x2, y, line_y, color, y_offset, sync_value, significance_level=0.05):
    """Add significance bars to a boxplot."""
    ax.plot([x1, x1, x2, x2], [line_y, line_y + y_offset, line_y + y_offset, line_y], color="k")
    
    # Annotate with asterisks based on p-value
   # if p_value > significance_level:
    value = pValue_cat(p_value, significance_level = 0.05)
    ax.text(( (x1 + x2) / 2 ), y, str(sync_value) + " " + value, color = "k", fontsize=16, ha='center')
    """
    if p_value < 0.05 and p_value > 0.01:
        ax.text(( (x1 + x2) / 2 ), y, str(sync_value) + '*', color = "k", fontsize=16, ha='center')
    if p_value < 0.01 and p_value > 0.001:
        ax.text(( (x1 + x2) / 2 ), y, str(sync_value) + '* *', color = "k", fontsize=16, ha='center')
    if p_value < 0.001:
        ax.text(( (x1 + x2) / 2 ), y, str(sync_value) + '* * *', color = "k", fontsize=16, ha='center')
    """
def sig_bar_param(ax, feature, kcats, y_pos, line_y, color, y_offset, sync_value, test_name):
    sig_test_kcats = kcats
    par_save = []
    for stest in sig_test_kcats:
       df = pd.read_csv("./" + str(stest) + ".csv")
       par_save.append(df[feature])
    if test_name == "t-test":
       print("Doing T-test")
       sig = stats.ttest_rel(par_save[0], par_save[1])
    if test_name == "wilcoxon":
       print("Doing Wilcoxon")
       sig = stats.wilcoxon(par_save[0], par_save[1])
    x_pos = []
    for i in sig_test_kcats:
        x_pos.append( kcat_list.index(i) )
    print(feature, kcats[0], kcats[1], sync_value, sig)
    sig_bars(ax, sig.pvalue, x_pos[0], x_pos[1], y_pos , line_y, color, y_offset,  sync_value, significance_level=0.05)

def just_the_tests(features, kcat, test_name):
       df = pd.read_csv("./" + str(kcat) + ".csv")
       if test_name == "t-test":
          print("Doing T-test")
          sig = stats.ttest_rel(df[features[0]], df[features[1]])
       if test_name == "wilcoxon":
          print("Doing Wilcoxon")
          sig = stats.wilcoxon(df[features[0]], df[features[1]])
       return sig.pvalue


def violin_split(list1, list2, var1, var2, group_name, parameter):
    arranged_list = []
    arranged_var = []
    arranged_group = []
    df = pd.DataFrame()
    for l1 in range(len(list1)):
        if list1[l1] > 0.01:
           arranged_list.append(list1[l1])
        else:
           arranged_list.append(np.nan)
        arranged_var.append(var1)
        arranged_group.append(group_name)
    for l2 in range(len(list2)):
        if list2[l2] > 0.01:
           arranged_list.append(list2[l2])
        else:
           arranged_list.append(np.nan)
        arranged_var.append(var2)
        arranged_group.append(group_name)
    df[parameter] = arranged_list
    df["Region"] = arranged_var
    df["kcat"] = arranged_group
    return df     



fig = plt.figure(figsize=(11.4, 14), layout="constrained")
spec = fig.add_gridspec(4, 3)
ax1 = fig.add_subplot(spec[0, 0])
ax2 = fig.add_subplot(spec[0, 1:2])
#ax3 = fig.add_subplot(spec[0, 2])
ax4 = fig.add_subplot(spec[1, :])
ax5 = fig.add_subplot(spec[2, :])
ax6 = fig.add_subplot(spec[3, :])

ax1.text(-0.2, 1.1, 'A',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax1.transAxes)
ax2.text(-0.2, 1.1, 'B',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax2.transAxes)
#ax3.text(-0.2, 1.1, 'C',
#        horizontalalignment='left',
#        verticalalignment='bottom',
#        transform=ax3.transAxes)
ax4.text(-0.2, 1.1, 'C',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax4.transAxes)
ax5.text(-0.2, 1.1, 'D',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax5.transAxes)
ax6.text(-0.2, 1.1, 'E',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax6.transAxes)

#Plotting spacingVs IBI
ax1_directory = "../chemistry_min_spacing_L50_no_delay/"
ax1_file = "spacingVsIBI_10.csv"
ax1Hnb = 10
df_ax1 = pd.read_csv(ax1_directory + ax1_file)
ax1_xticks = list(set(list(df_ax1["spacing"])))
ax1_xticks = sorted(ax1_xticks)
ax1_yticks = list(set(list(df_ax1["freq"])))
ax1_yticks = sorted(ax1_yticks)
print("FIRS PLOT DATA:", df_ax1)
print(ax1_xticks)
print(ax1_yticks)
#myColors = ["#0000ff", "#ffff00","#00ff00"]
myColors = ["#ffff00","#00ff00"]
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

ax1_yticks = np.round(1 / np.asarray(ax1_yticks), 1)
res = df_ax1.pivot(index = "freq", columns = "spacing", values = "N")
ax1 = sns.heatmap(res, ax = ax1, cmap = cmap, annot = False, cbar_kws={'ticks': [1, 2]}, xticklabels = ax1_xticks, yticklabels = ax1_yticks)
ax1.set_xticklabels(ax1.get_xticklabels(), rotation = 90)
ax1.set_yticklabels(ax1.get_yticklabels(), rotation = 0)
#ax1.invert_yaxis()
ax1.set_title(str(ax1Hnb) + " bursts")
ax1.set_xlabel("Spine Spacing $\mu m$")
ax1.set_ylabel("Inter burst \n interval (s)")

print("DATA: ", res.values)

"""
zm = np.ma.masked_less(res.values, 0.4)
x= np.arange(len(res.columns)+1)
y= np.arange(len(res.index)+1)
ax1.pcolor(x, y, zm, hatch='1', alpha=0.)
"""

files = ["onlyBG.csv", "50.0.csv"]
combined_df_life = pd.DataFrame()
combined_df_isd = pd.DataFrame()

onlyBG_dir = "./onlyBG_compare/"

sig_bg_sync = []
for file in files:
   df = pd.read_csv(onlyBG_dir + file)
   if file == "onlyBG.csv":
      group = "S1"
      print(group)
   if file == "50.0.csv":
      group = "S2"
   df_violin = violin_split(list(df["sync_life"]), list( df["out_sync_life"] ), "Rs", "Rns", group, "life")
   sig_bg_sync.append(stats.wilcoxon(df["sync_life"], df["out_sync_life"]).pvalue)
   combined_df_life = pd.concat([combined_df_life, df_violin], ignore_index = True)
   df_violin = violin_split(list( np.asarray(df["sync_isd"]) * 1e6 ), list( np.asarray(df["out_sync_isd"]) * 1e6 ), "Rs", "Rns", group, "ISD")
   combined_df_isd = pd.concat([combined_df_isd, df_violin], ignore_index = True)
print("pValues bgVssync: ", sig_bg_sync)   
ax2 = sns.violinplot(x="kcat", y="life", hue="Region", data=combined_df_life, split=True, ax = ax2, inner = "quart")
ax2.set_xlabel("")
ax2.set_ylabel("Life time (s)")
ax2.legend(frameon = False)
bg_sync_p_count = 0
for pV in sig_bg_sync:
   pValue = pValue_cat(pV, significance_level = 0.05)
   ax2.text(bg_sync_p_count, 65, pValue, color = "k", fontsize=16, ha='center')
   bg_sync_p_count += 1
#ax2 = sns.violinplot(x="kcat", y="ISD", hue="Region", data=combined_df_isd, split=True, ax = ax2, inner = "quart")

def make_violin_data(parameter1, parameter2, name, axis):
   sync_medians = []
   out_sync_medians = []
   df_20 = pd.read_csv("./" + str(20.0) + ".csv")
   if name == "life":
      df_violin_20 = violin_split(list(df_20[parameter1]), list( df_20[parameter2] ), "Rs", "Rns", str(20.0), name)
      sync_medians.append(np.median(df_20[parameter1]))
      out_sync_medians.append(np.median(df_20[parameter2]))
   if name == "ISD":
      df_violin_20 = violin_split(list(np.asarray(df_20[parameter1]) * 1e6), list( np.asarray(df_20[parameter2]) * 1e6 ), "Rs", "Rns", str(20.0), name)
      sync_medians.append(np.median(df_20[parameter1]) * 1e6)
      out_sync_medians.append(np.median(df_20[parameter2]) * 1e6)
   if name == "Area":
      df_violin_20 = violin_split(list(np.asarray(df_20[parameter1]) * 1e12), list( np.asarray(df_20[parameter2]) * 1e12 ), "Rs", "Rns", str(20.0), name)
      sync_medians.append(np.median(df_20[parameter1]) * 1e12)
      out_sync_medians.append(np.median(df_20[parameter2]) * 1e12)
   if name == "Spine density":
      df_violin_20 = violin_split(list(df_20[parameter1]), list( df_20[parameter2] ), "Rs", "Rns", str(20.0), name)
      print("Spine density")
      print(df_violin_20)
      print(df_20["out_sync_sp_density"])   
   """
   df_30 = pd.read_csv("./" + str(30.0) + ".csv")
   if name == "life":
      df_violin_30 = violin_split(list(df_30[parameter1]), list( df_30[parameter2] ), "+Sync", "-Sync", str(30.0), name)
      sync_medians.append(np.median(df_30[parameter1]))
      out_sync_medians.append(np.median(df_30[parameter2]))
   if name == "ISD":
      df_violin_30 = violin_split(list(np.asarray(df_30[parameter1]) * 1e6), list( np.asarray(df_30[parameter2]) * 1e6 ), "Rs", "-Rns", str(30.0), name)
      sync_medians.append(np.median(df_30[parameter1]) * 1e6)
      out_sync_medians.append(np.median(df_30[parameter2]) * 1e6)
   if name == "Area":
      df_violin_30 = violin_split(list(np.asarray(df_30[parameter1]) * 1e12), list( np.asarray(df_30[parameter2]) * 1e12 ), "Rs", "Rns", str(30.0), name)
      sync_medians.append(np.median(df_30[parameter1]) * 1e12)
      out_sync_medians.append(np.median(df_30[parameter2]) * 1e12)
   combined_df = pd.concat([df_violin_20, df_violin_30], ignore_index = True)
   """
   combined_df = df_violin_20


   #for kcat in [40.0, 50.0, 60.0, 80.0, 100.0, 200.0]:
   for kcat in [50.0, 200.0]:
      df = pd.read_csv("./" + str(kcat) + ".csv")
      if name == "life":
         df_violin = violin_split(list(df[parameter1]), list( df[parameter2] ), "Rs", "Rns", str(kcat), name)
         sync_medians.append(np.median(df[parameter1]))
         out_sync_medians.append(np.median(df[parameter2]))
      if name == "ISD":
         df_violin = violin_split(list(np.asarray(df[parameter1]) * 1e6), list( np.asarray(df[parameter2] ) * 1e6), "Rs", "Rns", str(kcat), name)
         sync_medians.append(np.median(df[parameter1]) * 1e6)
         out_sync_medians.append(np.median(df[parameter2]) * 1e6)
      if name == "Area":
         df_violin = violin_split(list(np.asarray(df[parameter1]) * 1e12), list( np.asarray(df[parameter2] ) * 1e12), "Rs", "Rns", str(kcat), name)
         sync_medians.append(np.median(df[parameter1]) * 1e12)
         out_sync_medians.append(np.median(df[parameter2]) * 1e12)
      if name == "Spine density":
         df_violin = violin_split(list(df[parameter1]), list( df[parameter2] ), "Rs", "Rns", str(kcat), name)
      combined_df = pd.concat([combined_df, df_violin], ignore_index = True)
   
   #axis.plot(xticks, sync_medians, linestyle = "--")
   #axis.plot(xticks, out_sync_medians, linestyle = "--")
   combined_df.to_csv("./combined.csv")
   axis = sns.violinplot(x="kcat", y=name, hue="Region", data=combined_df, split=True, ax = axis, inner = "quart")
   #axis.set_xticks(list(xticks), kcat_list)
   #axis.plot(xticks, sync_medians, linestyle = "--")
   #axis.plot(xticks, out_sync_medians, linestyle = "--")
    
make_violin_data("sync_life", "out_sync_life", "life", ax4)

make_violin_data("sync_sp_density", "out_sync_sp_density", "Spine density", ax5)

make_violin_data("sync_area", "out_sync_area", "Area", ax6)

#ax2 = sns.violinplot(x="group", y="ISD", hue="BG", data=combined_df, split=True)
#ax2.set_title("Healthy")
#ax2.set_ylim([0, 8])
"""
sync_medians = []
out_sync_medians = []
df_20 = pd.read_csv("./" + str(20.0) + ".csv")
sync_medians.append(np.median(df_20["sync_isd"]) * 1e6)
out_sync_medians.append(np.median(df_20["out_sync_isd"]) * 1e6)
df_violin_20 = violin_split(list(np.asarray(df_20["sync_isd"]) * 1e6), list( np.asarray(df_20["out_sync_isd"]) * 1e6 ), "+Sync", "-Sync", str(20.0), "ISD")
df_30 = pd.read_csv("./" + str(30.0) + ".csv")
sync_medians.append(np.median(df_30["sync_isd"]) * 1e6)
out_sync_medians.append(np.median(df_30["out_sync_isd"]) * 1e6)
df_violin_30 = violin_split(list(np.asarray(df_30["sync_isd"]) * 1e6), list( np.asarray(df_30["out_sync_isd"]) * 1e6 ), "+Sync", "-Sync", str(30.0), "ISD")
combined_df_isd = pd.concat([df_violin_20, df_violin_30], ignore_index = True)

for kcat in [40.0, 50.0, 60.0, 80.0, 100.0, 200.0]:
   df = pd.read_csv("./" + str(kcat) + ".csv")
   sync_medians.append(np.median(df["sync_isd"]) * 1e6)
   out_sync_medians.append(np.median(df["out_sync_isd"]) * 1e6)
   df_violin = violin_split(list(np.asarray(df["sync_isd"]) * 1e6), list( np.asarray(df["out_sync_isd"] ) * 1e6), "+Sync", "-Sync", str(kcat), "ISD")
   combined_df_isd = pd.concat([combined_df_isd, df_violin], ignore_index = True)
   
combined_df_isd.to_csv("./combined_isd.csv")
ax5 = sns.violinplot(x="kcat", y="ISD", hue="Region", data=combined_df_isd, split=True, ax = ax5, inner = "quart")
ax5.set_xticks(xticks, kcat_list)
ax5.plot(xticks, sync_medians, linestyle = "--")
ax5.plot(xticks, out_sync_medians, linestyle = "--")
"""

test_name = "wilcoxon"
ax4.set_title(test_name)
ax5.set_title(test_name)
kcat_selects = [[20.0, 50.0], [50.0, 200.0]]
y_offset = 5.0
y_pos = 70
line_y = 70
sig_bar_param(ax4, "sync_life", kcat_selects[0], y_pos, line_y,  "blue", y_offset, "Rs", test_name)
sig_bar_param(ax4, "out_sync_life", kcat_selects[0], y_pos + 5, line_y, "orange", y_offset, "Rns", test_name)

sig_bar_param(ax4, "sync_life", kcat_selects[1], y_pos, line_y, "blue", y_offset, "Rs", test_name)
sig_bar_param(ax4, "out_sync_life", kcat_selects[1], y_pos + 5, line_y, "orange", y_offset, "Rns", test_name)

kcat_selects = [[20.0, 200.0]]
y_offset = -5.0
y_pos = -1.0
line_y = -2.0
#sig_bar_param(ax4, "sync_life", kcat_selects[0], y_pos, line_y,  "blue", y_offset, "Rs", test_name)
#sig_bar_param(ax4, "out_sync_life", kcat_selects[0], y_pos - 5, line_y, "orange", y_offset, "Rns", test_name)

#sig_bar_param(ax1, "sync_life", kcat_selects[1], y_pos, line_y, "blue", y_offset, "+Sync", test_name)
#sig_bar_param(ax1, "out_sync_life", kcat_selects[1], y_pos - 5, line_y, "orange", y_offset, "-Sync", test_name)
kcat_selects = [[20.0, 50.0], [50.0, 200.0]]
y_pos = 0.8
line_y = 0.8
y_offset = 0.1
sig_bar_param(ax5, "sync_sp_density", kcat_selects[0], y_pos, line_y, "blue", y_offset, "Rs", test_name)
sig_bar_param(ax5, "out_sync_sp_density", kcat_selects[0], y_pos + 0.1, line_y, "orange", y_offset, "Rns", test_name)

sig_bar_param(ax5, "sync_sp_density", kcat_selects[1], y_pos, line_y, "blue", y_offset, "Rs", test_name)
sig_bar_param(ax5, "out_sync_sp_density", kcat_selects[1], y_pos + 0.1, line_y, "orange", y_offset, "Rns", test_name)

kcat_selects = [[20.0, 200.0]]
y_offset = -0.3
y_pos = -0.3
line_y = -0.1

#sig_bar_param(ax5, "sync_sp_density", kcat_selects[0], y_pos, line_y, "blue", y_offset, "Rs", test_name)
#sig_bar_param(ax5, "out_sync_sp_density", kcat_selects[0], y_pos + 0.1, line_y, "orange", y_offset, "Rns", test_name)

#sig_bar_param(ax2, "sync_isd", kcat_selects[1], y_pos, line_y,  "blue", y_offset, "+Sync", test_name)
#sig_bar_param(ax2, "out_sync_isd", kcat_selects[1], y_pos + 0.5, line_y, "orange", y_offset, "-Sync", test_name)
kcat_selects = [[20.0, 50.0], [50.0, 200.0]]
y_pos = 0.2
line_y = 0.2
y_offset = 0.01
sig_bar_param(ax6, "sync_area", kcat_selects[0], y_pos, line_y, "blue", y_offset, "Rs", test_name)
sig_bar_param(ax6, "out_sync_area", kcat_selects[0], y_pos + 0.1, line_y, "orange", y_offset, "Rns", test_name)

sig_bar_param(ax6, "sync_area", kcat_selects[1], y_pos, line_y, "blue", y_offset, "Rs", test_name)
sig_bar_param(ax6, "out_sync_area", kcat_selects[1], y_pos + 0.1, line_y, "orange", y_offset, "Rns", test_name)


kcat_selects = [[20.0, 200.0]]
y_offset = -0.01
y_pos = 0.1
line_y = 0.1

#sig_bar_param(ax6, "sync_area", kcat_selects[0], y_pos, line_y, "blue", y_offset, "Rs", test_name)
#sig_bar_param(ax6, "out_sync_area", kcat_selects[0], y_pos + 0.01, line_y, "orange", y_offset, "Rns", test_name)

ax5.set_ylim([0, 1.0])
ax4.legend(loc = 'center right', frameon = False)
ax4.set_ylabel("Life time (s)")
#ax5.set_yticks([0, 2, 4, 6, 8])
ax5.legend(loc = 'upper left', frameon = False)
ax2.spines[['right', 'top']].set_visible(False)
ax4.spines[['right', 'top']].set_visible(False)
ax5.spines[['right', 'top']].set_visible(False)
ax6.spines[['right', 'top']].set_visible(False)
ax6.legend(frameon = False)
#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
getattr(plt.figure(1), '_pylustrator_init', lambda: ...)()
plt.figure(1).ax_dict["<colorbar>"].set(position=[0.4248, 0.7975, 0.007714, 0.1254])
plt.figure(1).axes[0].set(position=[0.1691, 0.8037, 0.2444, 0.15])
plt.figure(1).axes[1].set(position=[0.5655, 0.8037, 0.3834, 0.15])
plt.figure(1).axes[1].texts[0].set(position=(-0.1966, 1.12))
plt.figure(1).axes[2].set(position=[0.1691, 0.542, 0.7798, 0.1254])
plt.figure(1).axes[2].texts[0].set(position=(-0.0732, 1.108))
plt.figure(1).axes[2].texts[1].set(position=(0.4718, 66.5))
plt.figure(1).axes[2].texts[2].set(position=(0.5, 79.25))
plt.figure(1).axes[2].texts[3].set(position=(1.463, 66.5))
plt.figure(1).axes[2].texts[4].set(position=(1.505, 79.25))
plt.figure(1).axes[3].set(position=[0.1691, 0.3157, 0.7798, 0.1254])
plt.figure(1).axes[3].get_legend().set(visible=False)
plt.figure(1).axes[3].texts[0].set(position=(-0.0732, 1.145))
plt.figure(1).axes[3].texts[1].set(position=(0.5, 0.7508))
plt.figure(1).axes[3].texts[2].set(position=(0.5, 0.9492))
plt.figure(1).axes[3].texts[3].set(position=(1.5, 0.7508))
plt.figure(1).axes[3].texts[4].set(position=(1.5, 0.9492))
plt.figure(1).axes[3].title.set(visible=False)
plt.figure(1).axes[4].set(position=[0.1691, 0.06701, 0.7798, 0.1254])
plt.figure(1).axes[4].get_legend().set(visible=False)
plt.figure(1).axes[4].texts[0].set(position=(-0.0732, 1.157))
plt.figure(1).axes[4].texts[1].set(position=(0.5, 0.1879))
plt.figure(1).axes[4].texts[2].set(position=(0.5, 0.2141))
plt.figure(1).axes[4].texts[3].set(position=(1.5, 0.1879))
plt.figure(1).axes[4].texts[4].set(position=(1.5, 0.2141))
#% end: automatic generated code from pylustrator
plt.show()






