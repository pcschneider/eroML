import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import math 
plt.rcParams.update({'font.size': 16})

PCTOCM = 3.0856776e18
#  Open color table 
colorsdata2 = ascii.read('../ero_data/EEM_dwarf_UBVIJHK_colors_Teff0319.txt')
temp = np.core.defchararray.find(colorsdata2.field('col4'), '...') != -1
colorsdata2.field('col4')[temp] = 'nan'
L_bol_theo = 10**(colorsdata2.field('col4').astype(np.float))*3.828e33
temp = np.core.defchararray.find(colorsdata2.field('col13'), '...') != -1
colorsdata2.field('col13')[temp] = 'nan'
M_G_theo = colorsdata2.field('col13').astype(np.float)
temp = np.core.defchararray.find(colorsdata2.field('col11'), '...') != -1
colorsdata2.field('col11')[temp] = 'nan'
BP_RP_theo = colorsdata2.field('col11').astype(np.float)
# create plot
L_X_sat = 10**(-3)*L_bol_theo
F_G_theo = 10**(-0.4*M_G_theo)*4052.97*2.5e-9
L_G_theo = 4*math.pi*10**2*PCTOCM**2*F_G_theo
activity_sat = np.log10(L_X_sat/L_G_theo)

fig, ax1 = plt.subplots()
ax1.plot(BP_RP_theo, activity_sat, 'k--')
plt.show()

Fx = np.ones(len(BP_RP_theo)) * 3.75e-15
Lx = 4*np.pi*(3.1e18 * 10)**2 * Fx # Lx at 10 pc
Lbol = 1e3 * Lx
Mg = -2.5*np.log10(Lbol / L_bol_theo) + M_G_theo
#Mg = M_G_theo
gi = np.where((BP_RP_theo > 0) & (BP_RP_theo < 4))[0]
pp = np.polyfit(BP_RP_theo[gi], Mg[gi], 4)

fig = plt.figure(figsize=(8,5))
fig.subplots_adjust(top=0.98, right=0.98, left=0.15, bottom=0.13)



#plt.plot(BP_RP_theo, Mg)
plt.plot(BP_RP_theo, np.polyval(pp, BP_RP_theo), lw=2)
ax1=plt.gca()
ypos = 19.7
ax1.axvline(x=-0.037, ymin=0, ymax=0.06, color='k', ls='-')
ax1.axvline(x=0.377, ymin=0, ymax=0.06, color='k', ls='-')
ax1.axvline(x=0.782, ymin=0, ymax=0.06, color='k', ls='-')
ax1.axvline(x=0.98, ymin=0, ymax=0.06, color='k', ls='-')
ax1.axvline(x=1.84, ymin=0, ymax=0.06, color='k', ls='-')
ax1.axvline(x=5, ymin=0, ymax=0.06, color='k', ls='-')
#ax1.text(-0.3, ypos, 'B', fontsize=12)
ax1.text(0.14, ypos, 'A', fontsize=14)
ax1.text(0.55, ypos, 'F', fontsize=14)
ax1.text(0.84, ypos, 'G', fontsize=14)
ax1.text(1.38, ypos, 'K', fontsize=14)
ax1.text(3.15, ypos, 'M', fontsize=14)

plt.xlim(-0.2, 4.2)
plt.ylim(19.8, 16.2)
plt.xlabel("BP - RP")
plt.ylabel("Gmag")
plt.show()
