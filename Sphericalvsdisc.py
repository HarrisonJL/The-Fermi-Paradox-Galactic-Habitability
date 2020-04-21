from pandas import *
from mpl_toolkits.mplot3d import Axes3D
from sympy.solvers import solveset
from sympy import Symbol
import numpy as np
import matplotlib.pyplot as plt
import collections
import math
import random
from numpy import exp, log10 as log







#Defining the fucntion which will find the fraction of stars that are affected by the LGRBs

class Galaxy:
    def __init__(self, v, ratio, n_stars, n_LGRBs, seed =1):
        if seed == -1:
            random.seed(None)
        else:
            random.seed(seed)
        self.n_LGRBs = n_LGRBs
        self.n_stars = n_stars
        self.v = v
        self.ratio = ratio
        self.r = ((self.v / self.ratio)*(3/(4*np.pi)))**(1/3)
        

        self.h = self.r * self.ratio
        print(self.r, self.h)
        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())

        star_phi = np.array(ran)*2*np.pi

       
        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*np.pi

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        a  = ((np.array(ran))**(1/3))* self.r
     

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        b  = ((np.array(ran))**(1/3))*self.r

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        c  = ((np.array(ran))**(1/3))*self.h




        self.stars_x = a*np.sin(star_theta)*np.cos(star_phi)
        self.stars_y = b*np.sin(star_theta)*np.sin(star_phi)
        self.stars_z = c*np.cos(star_theta)

        


    def effected_by_LGRB(self):

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())

        LGRB_phi = np.array(ran)*2*np.pi

       
        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())
        
        LGRB_theta  = np.array(ran)*np.pi

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())
        
        a  = ((np.array(ran))**(1/3))* self.r
     

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())
        
        b  = ((np.array(ran))**(1/3))*self.r

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())
        
        c  = ((np.array(ran))**(1/3))*self.h




        self.LGRB_x = a*np.sin(LGRB_theta)*np.cos(LGRB_phi)
        self.LGRB_y = b*np.sin(LGRB_theta)*np.sin(LGRB_phi)
        self.LGRB_z = c*np.cos(LGRB_theta)




        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())
        LGRB_x_grad = -1 + 2*np.array(ran)

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())

        LGRB_y_grad = -1 + 2*np.array(ran)

        ran = []
        for i in range(self.n_LGRBs):
            ran.append(random.random())

        LGRB_z_grad = -1 + 2*np.array(ran)

        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)




        normaliser = (LGRB_x_grad**2 + LGRB_y_grad**2 + LGRB_z_grad**2)**0.5

        LGRB_x_grad = LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = LGRB_z_grad*(1/normaliser)




        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.n_LGRBs):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                theta = np.arccos(costheta)
                dist2 = (self.LGRB_x[j]-star_x_check)**2 +  (self.LGRB_y[j]-star_y_check)**2  +  (self.LGRB_z[j]-star_z_check)**2

                

                if (theta < theta_max ) and  ( dist2 < 2.16**2):
                    effected_star_x.append(star_x_check)
                    effected_star_y.append(star_y_check)
                    effected_star_z.append(star_z_check)

        

        effected_stars = list(zip(effected_star_x,effected_star_y,effected_star_z))
        effected_stars_filtered = []

        #removing duplicates in effected stars array
        for x in effected_stars:
            if x not in effected_stars_filtered:
                effected_stars_filtered.append(x)



        all_stars = list(zip(self.stars_x, self.stars_y, self.stars_z))

        joined_stars = all_stars + effected_stars_filtered

        d = collections.defaultdict(int)
        for x in joined_stars: d[x] += 1
        uneffected_stars = [x for x in joined_stars if d[x] == 1]

        uneffected_stars_x, uneffected_stars_y, uneffected_stars_z = zip(*uneffected_stars)

        if len(uneffected_stars_x) == len(all_stars):
            pass
        else:

            effected_star_x, effected_star_y, effected_star_z = zip(*effected_stars_filtered)



        return(effected_star_x, effected_star_y, effected_star_z, self.LGRB_x, self.LGRB_y, self.LGRB_z, uneffected_stars_x, uneffected_stars_y, uneffected_stars_z, LGRB_x_grad, LGRB_y_grad, LGRB_z_grad)



            

effected_galaxy_cumulative = np.zeros(99)
effected_galaxy_for_different_phi_1 = []
effected_galaxy_for_different_phi_2 = []
effected_galaxy_for_different_phi_3 = []
effected_galaxy_for_different_phi_4 = []
effected_galaxy_for_different_phi_5= []
effected_galaxy_for_different_phi_6 = []
effected_galaxy_for_different_phi_7 = []
effected_galaxy_for_different_phi_8 = []
effected_galaxy_for_different_phi_9 = []
effected_galaxy_for_different_phi_10 = []


for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_1.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_1 = np.array(effected_galaxy_for_different_phi_1)
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_2.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)  
effected_galaxy_for_different_phi_2 = np.array(effected_galaxy_for_different_phi_2) 
    
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_3.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_3 = np.array(effected_galaxy_for_different_phi_3)

for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_4.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100) 
    
effected_galaxy_for_different_phi_4 = np.array(effected_galaxy_for_different_phi_4)
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_5.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_5 = np.array(effected_galaxy_for_different_phi_5)
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_6.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)   
    
effected_galaxy_for_different_phi_6 = np.array(effected_galaxy_for_different_phi_6)
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_7.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_7 = np.array(effected_galaxy_for_different_phi_7)
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_8.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_8 = np.array(effected_galaxy_for_different_phi_8)

for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_9.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100)
    
effected_galaxy_for_different_phi_9 = np.array(effected_galaxy_for_different_phi_9)
radius = []
height = []
for i in range(1, 100):
    print(i)
    galaxy = Galaxy(7500, (i/100)*1, 1000, 1000, -1)
    radius.append(galaxy.r)
    height.append(galaxy.h)
    galaxy_effected = galaxy.effected_by_LGRB()
    effected_galaxy_for_different_phi_10.append((len(galaxy_effected[0])/(len(galaxy_effected[6])+len(galaxy_effected[0])))*100) 
    
effected_galaxy_for_different_phi_10 = np.array(effected_galaxy_for_different_phi_10)



effected_galaxy_cumulative =   effected_galaxy_for_different_phi_1 + effected_galaxy_for_different_phi_2 + effected_galaxy_for_different_phi_3 + effected_galaxy_for_different_phi_4 + effected_galaxy_for_different_phi_5 + effected_galaxy_for_different_phi_6 +effected_galaxy_for_different_phi_7+effected_galaxy_for_different_phi_8+effected_galaxy_for_different_phi_9+effected_galaxy_for_different_phi_10


average = effected_galaxy_cumulative/10


ratio = np.round(np.array(height)/np.array(radius), 4)





plt.plot(ratio, average, linewidth = 0.8, color='b')



plt.title('Spherical vs Disc Galaxy')
plt.xlabel(r'Height/Radius', fontsize=20)
plt.ylabel(r'Percentage of stars effected', fontsize=20)
plt.tight_layout()
plt.savefig('Spherical vs Disc.png')

plt.show()





