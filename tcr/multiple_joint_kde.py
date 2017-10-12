import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set(style="darkgrid")
Pt_KDE = pd.read_csv('Pt31_KDE.csv')

# Subset the Pt dataset by Treatment
On = Pt_KDE.query("Treatment == 'On'")
Pre = Pt_KDE.query("Treatment == 'Pre'")

# Set up the figure
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect("auto")
axes = plt.gca()
#axes.set_xlim([0,150])
#axes.set_ylim([ymin,ymax])

#Draw the two density plots

#cmap = sns.cubehelix_palette(start=2.45, rot=0.05, light=0.85, as_cmap=True)
#cmap2 = sns.cubehelix_palette(start=0, rot=0.2, light=0.9, as_cmap=True)

ax = sns.kdeplot(On.Number_CDR3, On.H, n_levels=250,
                 cmap="Blues", shade=True, shade_lowest=False)
ax = sns.kdeplot(Pre.Number_CDR3, Pre.H, n_levels=250,
                 cmap="Reds", shade=True, shade_lowest=False)


# Add labels to the plot
#red = sns.color_palette("Reds")[-2]
#blue = sns.color_palette("Blues")[-2]
#ax.text(40, 1.1, "On", size=16, color=blue)
#ax.text(40, 0.7, "Pre", size=16, color=red)

plt.show()

fig.savefig('myimage.eps', format='eps', dpi=1200)
