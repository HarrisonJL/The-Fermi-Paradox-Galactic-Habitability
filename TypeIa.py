import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import pyplot

#Data for curve from Schaye et al. 2014
x = [1.010112124,1.027826038,1.054984083,1.08285972,1.117923605,1.150783368,1.181183522,1.212386754,1.248023126,1.288432677,1.322458997,1.349545109,1.381174603,1.417655567,1.446685934,1.476322056,1.502191612,1.546346383,1.587177903,1.640933726,1.696516675,1.779583779,1.855937938,1.949657774,2.018661157,2.102271193,2.170398352,2.22778466,2.273482798,2.340389039,2.402315864,2.458713236,2.523761285,2.590540151,2.651376582,2.697960274,2.745362424,2.793586739,2.850934639,2.909459799,2.977808515,3.038949804,3.110352368,3.165000092,3.201984672,3.267703848,3.334771879,3.393375534,3.483191055,3.617136283,3.702137945,3.789093694,3.889397903,4.015577572,4.086129764,4.182088729,4.292747032,4.419170979,4.536111157,4.710483173,4.849154121,4.948642836]
y = [0.223827565,0.237596864,0.255292223,0.272987582,0.302504814,0.330055895,0.353668432,0.377280969,0.40483205,0.436321675,0.467823783,0.491448804,0.522957154,0.550514476,0.578084282,0.597764518,0.625340565,0.652891646,0.68833854,0.727714615,0.763145904,0.806435556,0.84184812,0.873300294,0.889004535,0.900751507,0.904627633,0.904571457,0.896638194,0.888686205,0.868906102,0.860966597,0.84513128,0.825351176,0.8095221,0.789760723,0.769999345,0.754182752,0.730470347,0.706757942,0.683039295,0.655382104,0.627718672,0.607957294,0.5882084,0.56844078,0.54867316,0.524966997,0.497297323,0.457768325,0.426160107,0.406386246,0.374771786,0.343144843,0.323383465,0.307554389,0.287774286,0.258125977,0.236373481,0.210651233,0.192837281,0.180959232]




#errors 

#SDSS
SDSS_x = [1.027788735,1.095503004,1.130965487,1.188060293]
SDSS_y = [0.275072326,0.263100651,0.302479847,0.34773877]
SDSS_y_error_upper = [0.405250247,0.324247947,0.377424529,0.454247978]
SDSS_y_error_upper = np.array(SDSS_y_error_upper) - np.array(SDSS_y)
SDSS_y_error_lower = [0.204069308,0.221683525,0.270915321,0.30434613]
SDSS_y_error_lower = np.array(SDSS_y) - np.array(SDSS_y_error_lower)

#SDSS-DR7

SDSS_DR7_x = [1.09870099]
SDSS_DR7_y = [0.247315268]
SDSS_DR7_y_error_upper = [0.300566751]
SDSS_DR7_y_error_upper = np.array(SDSS_DR7_y_error_upper) - np.array(SDSS_DR7_y)
SDSS_DR7_y_error_lower = [0.194057543]
SDSS_DR7_y_error_lower = np.array(SDSS_DR7_y) - np.array(SDSS_DR7_y_error_lower)

#SNLS

SNLS_x = [1.273713615,1.361440006,1.461677347,1.562275039,1.662774716,1.764413,1.864369755,1.966673008,2.066122794]
SNLS_y = [0.323920255,0.412534369,0.412381446,0.550305378,0.548198787,0.664442003,0.652489054,0.883143532,0.84556196]
SNLS_y_error_upper = [0.47184659,0.544684682,0.53072501,0.662731764,0.650760089,0.780816296,0.782666975,1.090247891,1.096055839]
SNLS_y_error_upper = np.array(SNLS_y_error_upper) - np.array(SNLS_y)
SNLS_y_error_lower = [0.16218093,0.278417904,0.290093096,0.424069121,0.423938044,0.520460454,0.52033874,0.664207938,0.567454583]
SNLS_y_error_lower = np.array(SNLS_y) - np.array(SNLS_y_error_lower)

#GOODS

GOODS_x = [2.629416678]
GOODS_y = [0.415062277]
GOODS_y_error_upper = [0.994951985]
GOODS_y_error_upper = np.array(GOODS_y_error_upper) - np.array(GOODS_y)
GOODS_y_error_lower = [0.048197227]
GOODS_y_error_lower = np.array(GOODS_y) - np.array(GOODS_y_error_lower)

#SDF

SDF_x = [2.215071654,2.68179566]
SDF_y = [0.825688231,1.006736408]
SDF_y_error_upper = [1.089982617,1.547172019]
SDF_y_error_upper = np.array(SDF_y_error_upper) - np.array(SDF_y)
SDF_y_error_lower = [0.561381362,0.655650499]
SDF_y_error_lower = np.array(SDF_y) - np.array(SDF_y_error_lower)

#CLASH

CLASH_x = [1.413680415,1.933568918,2.580238306]
CLASH_y = [0.455845866,0.447282186,0.448633525]
CLASH_x_error_lower = [1,1.601415056,2.203053927]
CLASH_x_error_lower = np.array(CLASH_x) - np.array(CLASH_x_error_lower)
CLASH_x_error_upper = [1.587529631,2.196683464,2.778260488]
CLASH_x_error_upper = np.array(CLASH_x_error_upper) - np.array(CLASH_x)
CLASH_y_error_upper = [0.984453363,0.798374336,0.833253231]
CLASH_y_error_upper = np.array(CLASH_y_error_upper) - np.array(CLASH_y)
CLASH_y_error_lower = [0.014036133,0.198760701,0.133047565]
CLASH_y_error_lower = np.array(CLASH_y) - np.array(CLASH_y_error_lower)



plt.plot(x, y, label = 'RefL0100N1504', color = 'k', linewidth = 1.5)

plt.errorbar(SDSS_x, SDSS_y, yerr = (SDSS_y_error_lower,SDSS_y_error_upper), fmt = 'd', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'SDSS'  )

plt.errorbar(SDSS_DR7_x, SDSS_DR7_y, yerr = (SDSS_DR7_y_error_lower,SDSS_DR7_y_error_upper), fmt = 'o', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'SDSS-DR7')

plt.errorbar(SNLS_x, SNLS_y, yerr = (SNLS_y_error_lower,SNLS_y_error_upper), fmt = 'x', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'SNLS')

plt.errorbar(GOODS_x, GOODS_y, yerr = (GOODS_y_error_lower,GOODS_y_error_upper), fmt = 's', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'GOODS')

plt.errorbar(SDF_x, SDF_y, yerr = (SDF_y_error_lower,SDF_y_error_upper), fmt = 'p', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'SDF')

plt.errorbar(CLASH_x, CLASH_y, yerr = (CLASH_y_error_lower,CLASH_y_error_upper), xerr = (CLASH_x_error_lower,CLASH_x_error_upper), fmt = '*', capsize = 2, elinewidth = 0.8, markersize = 5, label = 'CLASH')


ax = plt.gca()
ax.set_xlim([1,4])
ax.set_ylim([0,1.7])
ax.set_xscale('log')
ax.legend(frameon=False, prop={'size': 9})


x_tick_locations = [1, 2, 3, 4]
x_tick_labels = ["0", "1.0", "2.0", "3.0"]

y_tick_locations = [0,0.5,1.0,1.5]
y_tick_labels = [0,0.5,1.0,1.5]

ax.set_xticks(x_tick_locations)
ax.set_xticklabels(x_tick_labels)
ax.set_yticks(y_tick_locations)
ax.set_yticklabels(y_tick_labels)

ax2 = ax.twiny()
ax2.set_xscale('log')
ax2.set_xlabel(r'Look back time [Gyrs]', fontsize=18)
ax2.set_xlim([1,4])
ax2.set_xticks([1.011295688,1.083943321,1.169228195,1.257249456,1.369070304,1.490796726,1.638831538,1.818751642,2.037515132,2.34107729,2.724108451,3.376139586])
ax2.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11])




ax.set_xlabel(r'Redshift', fontsize=18)
ax.set_ylabel(r'Type Ia rate [10$^{-4}$ yr$^{-1}$ cMpc$^{-3}$]', fontsize=18)
plt.tight_layout()
plt.savefig('TypeIaSNR.png')
plt.show()


