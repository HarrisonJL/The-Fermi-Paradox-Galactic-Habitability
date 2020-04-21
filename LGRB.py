import numpy as np
import eagleSqlTools as sql
import matplotlib.pyplot as plt
import math
from matplotlib import pyplot
from scipy.integrate import quad
from numpy import exp, log10 as log
from scipy.integrate import quad


# Array of chosen simulations. Entries refer to the simulation name and comoving box length.
mySims = np.array([('RefL0100N1504', 100.)])
# This uses the easgleSQLTools module to connect to the database with your username and password.
# If the password is not given, the module will prompt for it.
con = sql.connect("jwx596", password="SPS530AQ")

#FromEAGLE
Lowfrac = np.array([0.01760045000900897, 0.02063130449373665, 0.022747567398593845, 0.025302491842799844, 0.027833066098888717, 0.032300326716953855, 0.03523550089243264, 0.03887629689478371, 0.042855860942132756, 0.04598091296152274, 0.05129097342649215])*5.68167
Highfrac = np.array([1,1,1,1,1,1,1,1,1,1,1]) - Lowfrac


for sim_name , sim_size in mySims:
    print (sim_name)
    myQuery = "SELECT \
                    SH.Redshift as z , \
                    SUM(SH.StarFormationRate) as SFR \
                FROM \
                    %s_SubHalo as SH \
                WHERE \
                    SH.Redshift >= 0 and SH.Redshift <= 21 \
                GROUP BY \
                    SH.Redshift \
                ORDER BY \
                    z"%(sim_name)


myData = sql.execute_query(con, myQuery)

SFR_density = myData['SFR'][:] / (float(sim_size)**3)



LowZfraction = (8e-5)*(myData['z'][:]**6) - 0.0021*(myData['z'][:]**5) + 0.0238*(myData['z'][:]**4) - 0.1322*(myData['z'][:]**3) + 0.3407*(myData['z'][:]**2) - 0.1066*myData['z'][:] + 0.1165
for i in range(len(LowZfraction)):
    if LowZfraction[i] >=1:
        LowZfraction[i] = 1
print(LowZfraction)

LGRBrateEAGLE = (0.007*0.22*0.69*0.21*0.13)*SFR_density[0:11]*1e9 * (Lowfrac +(Highfrac)/25)

#LGRBrate = (0.007*0.22*0.69*0.21*0.13)*SFR_density*1e9 * (LowZfraction +(1-LowZfraction)/25)
LGRBrate_high = LGRBrateEAGLE + LGRBrateEAGLE*(193/224)
LGRBrate_low =  LGRBrateEAGLE - LGRBrateEAGLE*(193/224)




plt.plot(myData['z'][0:11], LGRBrate_high, label = 'RefL0100N1504 uncertainty', color = 'k', linewidth = 1.5, linestyle = '--')
plt.plot(myData['z'][0:11], LGRBrate_low, color = 'k', linewidth = 1.5, linestyle = '--')
plt.plot(myData['z'][0:11], LGRBrateEAGLE, label = 'RefL0100N1504', color = 'k', linewidth = 1.5)

ax = plt.gca()
ax.fill_between(myData['z'][0:11], LGRBrate_low, LGRBrate_high, alpha =0.2, color = 'k')

#data from wanderman & piran

redshift = [0,0.016985138,0.050955414,0.08492569,0.118895966,0.152866242,0.186836518,0.220806794,0.25477707,0.288747346,0.33970276,0.390658174,0.441613588,0.475583864,0.526539278,0.560509554,0.59447983,0.628450106,0.662420382,0.713375796,0.76433121,0.798301486,0.8492569,0.900212314,0.93418259,0.985138004,1.036093418,1.087048832,1.138004246]
rate = [0.423, 0.452624132,0.488117738,0.513313296,0.539809393,0.58213984,0.627789731,0.660194848,0.694272646,0.7487157,0.787362748,0.849105685,0.915690343,0.987496396,1.038468775,1.092072236,1.148442589,1.207722655,1.270062627,1.33562045,1.404562223,1.477062618,1.553305323,1.633483508,1.717800313,1.806469365,1.899715318,1.99777442,2.100895117]
error_x = np.array([0.244003588,0.74362499])
error_y = np.array([0.00541612,0.031654601]) /0.0055
error_y_upper = np.array([0.015737266,0.045565131])  /0.0055
error_y_upper = np.array(error_y_upper) - np.array(error_y)
error_y_lower = np.array([0.000972542,0.02087349])  /0.0055
error_y_lower = np.array(error_y) - np.array(error_y_lower)
error_x_upper = np.array([0.013070373,0.498258449])
error_x_upper = np.array(error_x_upper) - np.array(error_x)
error_x_lower = np.array([0.489370129,0.988991531])
error_x_lower = np.array(error_x) - np.array(error_x_lower)

#rate is the aligned rate, so need to divide by beaming fraction 

rate = np.array(rate) / 0.0055

plt.plot(redshift, rate, label = 'Graff & Lien 2016', color = 'b', linewidth = 1.5)
#plt.errorbar(error_x, error_y, yerr = (error_y_lower,error_y_upper), xerr = (error_x_lower, error_x_upper), fmt = 'p', capsize = 2, elinewidth = 0.8, markersize = 5, ecolor ='b', markerfacecolor = 'b' )

ax = plt.gca()
ax.set_xlim([0,1])
ax.set_ylim([1,1e3])
ax.set_yscale('log')
ax2 = ax.twiny()
ax2.set_xlabel(r'Look back time [Gyrs]', fontsize=18)
ax2.set_xticks([0,0.08,0.17,0.26,0.37,0.49,0.64,0.82,1.0])
ax2.set_xticklabels([0,1,2,3,4,5,6,7,8])

ax.set_xlabel(r'Redshift', fontsize=18)
ax.set_ylabel(r'LGRB rate [yr$^{-1}$ cGpc$^{-3}$]', fontsize=18)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig('LGRBrate.png')
plt.show()


number_TypeIa = 1e6*1e9*1e-4*np.trapz(([0,0.07,0.15,0.25,0.36,0.48]),x = ([0.24,0.48,0.72,1.005,1.299,1.632]))

numberLGRBs = 1e-9*1e9*1e6*np.trapz(LGRBrateEAGLE[0:6], x = ([0,1.35,2.32,3.24,4.12,5.22]))
print(number_TypeIa)
print(numberLGRBs)