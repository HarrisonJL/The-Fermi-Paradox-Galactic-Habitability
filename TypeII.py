import numpy as np
import eagleSqlTools as sql
import matplotlib.pyplot as plt
import math
from matplotlib import pyplot
from scipy.integrate import quad
from numpy import exp, log10 as log
from scipy.integrate import quad



#Introducing our definitions of Chabrier's IMFs for binary and individual star systems
def chabrier03individual(m):
	k = 0.158 * exp(-(-log(0.08))**2/(2 * 0.69**2))
	return np.where(m <= 1,\
	        0.158*(1./m) * exp(-(log(m)-log(0.08))**2/(2 * 0.69**2)),\
	        k*m**-2.3)

def mchabrier03individual(m):
	k = 0.158 * exp(-(-log(0.08))**2/(2 * 0.69**2))
	return np.where(m <= 1,\
	        m*0.158*(1./m) * exp(-(log(m)-log(0.08))**2/(2 * 0.69**2)),\
	        k*m**-2.3)


# Array of chosen simulations. Entries refer to the simulation name and comoving box length.
mySims = np.array([('RefL0100N1504', 100.)])
# This uses the easgleSQLTools module to connect to the database with your username and password.
# If the password is not given, the module will prompt for it.
con = sql.connect("jwx596", password="SPS530AQ")

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
#Finding the star formation rate density at each snapshot
SFR_density = myData['SFR'][:] / (float(sim_size)**3)

#referring to the equation from the paper - finding the numerator and denominator of the equation
numerator = quad(chabrier03individual, 8, 40)[0] - quad(chabrier03individual, 8, 40)[1]
denominator = quad(mchabrier03individual, 0.1, 100)[0] - quad(mchabrier03individual, 0.1, 100)[1]

Rate_TypeII = SFR_density * (numerator/denominator)


#adding in Taylor et al data

redshift = [0.000966184,0.018357488,0.041545894,0.058937198,0.09178744,0.13236715,0.15942029,0.201932367,0.248309179,0.312077295,0.350724638,0.385507246,0.426086957,0.4647343,0.51884058,0.576811594,0.648309179,0.715942029,0.758454106,0.810628019,0.870531401,0.91884058,0.969082126,0.996135266,1.007729469]
rate = [6.55934E-05,7.06053E-05,7.65105E-05,8.1807E-05,9.4786E-05,0.000111304,0.000122239,0.000143541,0.000168556,0.000208818,0.00023873,0.000265717,0.00029974,0.000333623,0.000391763,0.000460035,0.000562341,0.00066034,0.000734986,0.000829096,0.000960634,0.001083637,0.001190094,0.001272481,0.00130701]

plt.plot(redshift, rate, label = 'Taylor et al.', color = 'r', linewidth = 1.5)

Li_x = [0.005]
Li_y = [6.96664E-05]
Li_y_error_upper = [8.98439E-05]
Li_y_error_upper = np.array(Li_y_error_upper) - np.array(Li_y)
Li_y_error_lower = [5.18939E-05]
Li_y_error_lower = np.array(Li_y) - np.array(Li_y_error_lower)

plt.errorbar(Li_x, Li_y, yerr = (Li_y_error_lower,Li_y_error_upper), fmt = 'd', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Li'  )

Cappellaro_1999_x = [0.005]
Cappellaro_1999_y = [4.18883E-05]
Cappellaro_1999_y_error_upper = [6.17586E-05]
Cappellaro_1999_y_error_upper = np.array(Cappellaro_1999_y_error_upper) - np.array(Cappellaro_1999_y)
Cappellaro_1999_y_error_lower = [2.23274E-05]
Cappellaro_1999_y_error_lower = np.array(Cappellaro_1999_y) - np.array(Cappellaro_1999_y_error_lower)

plt.errorbar(Cappellaro_1999_x, Cappellaro_1999_y, yerr = (Cappellaro_1999_y_error_lower,Cappellaro_1999_y_error_upper), fmt = 'p', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Cappellaro 1999'  )

Botticella_x = [0.209661836]
Botticella_y = [0.000128963]
Botticella_y_error_upper = [0.00019793]
Botticella_y_error_upper = np.array(Botticella_y_error_upper) - np.array(Botticella_y)
Botticella_y_error_lower = [6.2591E-05]
Botticella_y_error_lower = np.array(Botticella_y) - np.array(Botticella_y_error_lower)

plt.errorbar(Botticella_x, Botticella_y, yerr = (Botticella_y_error_lower,Botticella_y_error_upper), fmt = '*', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Botticella'  )

Cappellaro_2005_x = [0.257971014]
Cappellaro_2005_y = [0.000180225]
Cappellaro_2005_y_error_upper = [0.000269298]
Cappellaro_2005_y_error_upper = np.array(Cappellaro_2005_y_error_upper) - np.array(Cappellaro_2005_y)
Cappellaro_2005_y_error_lower = [8.86491E-05]
Cappellaro_2005_y_error_lower = np.array(Cappellaro_2005_y) - np.array(Cappellaro_2005_y_error_lower)

plt.errorbar(Cappellaro_2005_x, Cappellaro_2005_y, yerr = (Cappellaro_2005_y_error_lower,Cappellaro_2005_y_error_upper), fmt = '1', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Cappellaro 2005'  )

Bazin_x = [0.298550725]
Bazin_y = [0.000182653]
Bazin_y_error_upper = [0.000235555]
Bazin_y_error_upper = np.array(Bazin_y_error_upper) - np.array(Bazin_y)
Bazin_y_error_lower = [0.000136057]
Bazin_y_error_lower = np.array(Bazin_y) - np.array(Bazin_y_error_lower)

plt.errorbar(Bazin_x, Bazin_y, yerr = (Bazin_y_error_lower,Bazin_y_error_upper), fmt = '2', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Bazin', color ='k'  )

Dahlen_x = [0.30821256,0.696618357]
Dahlen_y = [0.000282038,0.000445105]
Dahlen_y_error_upper = [0.000512038,0.000765105]
Dahlen_y_error_upper = np.array(Dahlen_y_error_upper) - np.array(Dahlen_y)
Dahlen_y_error_lower = [5.25933E-05,0.000130701]
Dahlen_y_error_lower = np.array(Dahlen_y) - np.array(Dahlen_y_error_lower)
Dahlen_x_error_lower = [0.109178744,0.499516908]
Dahlen_x_error_lower = np.array(Dahlen_x) - np.array(Dahlen_x_error_lower)
Dahlen_x_error_upper = [0.509178744,0.899516908]
Dahlen_x_error_upper = np.array(Dahlen_x_error_upper) - np.array(Dahlen_x)

plt.errorbar(Dahlen_x, Dahlen_y, yerr = (Dahlen_y_error_lower,Dahlen_y_error_upper), xerr = (Dahlen_x_error_lower, Dahlen_x_error_upper), fmt = 'o', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Dahlen'  )


SVISS_x = [0.404830918,0.725603865]
SVISS_y = [0.000397043,0.000754931]
SVISS_y_error_upper = [0.000829096,0.001435413]
SVISS_y_error_upper = np.array(SVISS_y_error_upper) - np.array(SVISS_y)
SVISS_y_error_lower = [0.000120613,0.000361526]
SVISS_y_error_lower = np.array(SVISS_y) - np.array(SVISS_y_error_lower)

plt.errorbar(SVISS_x, SVISS_y, yerr = (SVISS_y_error_lower,SVISS_y_error_upper), fmt = '3', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'SVISS'  )


Graur_x = [0.652173913]
Graur_y = [0.000807192]
Graur_y_error_upper = [0.001901384]
Graur_y_error_upper = np.array(Graur_y_error_upper) - np.array(Graur_y)
Graur_y_error_lower = [0.000101384]
Graur_y_error_lower = np.array(Graur_y) - np.array(Graur_y_error_lower)

plt.errorbar(Graur_x, Graur_y, yerr = (Graur_y_error_lower,Graur_y_error_upper), fmt = '4', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'Graur'  )

plt.plot(myData['z'][:], Rate_TypeII, label = 'RefL0100N1504', color = 'k', linewidth = 1.5)






ax = plt.gca()
ax.set_xlim([0,1])
ax.set_ylim([1e-5,1e-2])
ax.set_yscale('log')
ax.legend(frameon=False, prop={'size': 8})

ax2 = ax.twiny()
ax2.set_xlabel(r'Look back time [Gyrs]', fontsize=18)
ax2.set_xticks([0,0.08,0.17,0.26,0.37,0.49,0.64,0.82,1.0])
ax2.set_xticklabels([0,1,2,3,4,5,6,7,8])

ax.set_xlabel(r'Redshift', fontsize=18)
ax.set_ylabel(r'Type II rate [yr$^{-1}$ cMpc$^{-3}$]', fontsize=18)


plt.tight_layout()
plt.savefig('TypeIISNR.png')
plt.show()

#The number of Type II's within the past 5 Gyrs
numberTypeII = np.trapz(Rate_TypeII[0:6], x = ([0,1.35,2.32,3.24,4.12,5.22]))*1e9*1e6

