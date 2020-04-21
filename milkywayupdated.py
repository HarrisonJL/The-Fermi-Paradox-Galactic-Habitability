import eagleSqlTools as sql
import numpy as np
import matplotlib.pyplot as plt
import math
import random
from matplotlib import pyplot
from scipy.integrate import quad
from pandas import *
from mpl_toolkits.mplot3d import Axes3D
from sympy.solvers import solve
from sympy import Symbol
import collections
from sympy import I, re

class Galaxy:
    def __init__(self, h_a, h_b, h_c, n_stars, seed =1):
        if seed == -1:
            random.seed(None)
        else:
            random.seed(seed)

        self.n_stars = n_stars
        

        self.h_a = h_b
        self.h_b = h_b
        self.h_c = h_c

        print(h_b)

       
        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())

        c_max =  np.amax(-(self.h_c) * np.log((1 - np.array(ran))))

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        a_max =  np.amax(-(self.h_a) * np.log((1 - np.array(ran))))

        b_max =  np.amax(-(self.h_b) * np.log((1 - np.array(ran))) )

        print(a_max)
        print(b_max)
        print(c_max)


        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        
        rho  = (1 /( (np.sin(star_theta)/a_max)**2 + (np.cos(star_theta)/b_max)**2))**(1/2)
        self.h_r = self.h_a


        rho = ( -(    (self.h_r ) * np.log((1 - np.array(ran)**(1/2))) ))
        for i in range(len(rho)):
            while ( (b_max*rho[i]*np.sin(star_theta[i]))**2 + (a_max*rho[i]*np.cos(star_theta[i]))**2  - (a_max*b_max )**2 ) >= 0  :
                rho[i] = ( -(    (self.h_r ) * np.log((1 - random.random()**(1/2)))))
            else:
                pass
        









        self.stars_x = rho*np.sin(star_theta)
    
        self.stars_y = rho*np.cos(star_theta)



        
        self.stars_x_max = a_max

       


        self.stars_y_max = b_max

 
        z = ((1 - ((self.stars_x /a_max)**2)-((self.stars_y /b_max)**2))*(c_max**2))**(1/2)

        
        star_phi = np.arctan(rho/z)

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())
        


        self.stars_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran)))

        ran = []
        for i in range(self.n_stars):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        self.stars_z = self.stars_z * np.array(ran)

        self.stars_z_max = c_max


        self.volume = a_max * b_max* c_max *(4/3) *np.pi



    def effected_by_LGRB_step_1(self, Number_LGRBs_step_1, SFR_total, E, T):

        Number_LGRBs_step_1 = int(Number_LGRBs_step_1)
        self.Number_LGRBs_step_1 = Number_LGRBs_step_1

        ran = []
        for i in range(Number_LGRBs_step_1):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        ran = []
        for i in range(Number_LGRBs_step_1):
            ran.append(random.random())
        
        self.h_r = self.h_a


        rho = []

        for i in range(Number_LGRBs_step_1):
            r = Symbol("r")
            radius = solve((-(17000000*r**3-327000000*r**2+1044960000*r-921755351)/1200000000 - random.random()**(1/2)*5.306 ),cubic = True,complex=False)
            radius = np.array(radius[1])
            radius = re(radius+ (2*I))
            rho.append(radius)

        rho = np.array(rho)

        self.LGRB_x = rho*np.sin(star_theta)
    
        self.LGRB_y = rho*np.cos(star_theta)

        z = ((1 - ((self.LGRB_x /self.stars_x_max)**2)-((self.LGRB_y /self.stars_y_max)**2))*(self.stars_z_max**2))**(1/2)
 
        opp_adj = rho/z
        opp_adj = np.array(opp_adj, dtype = float)

        star_phi = np.arctan(opp_adj)

        ran = []
        for i in range(Number_LGRBs_step_1):
            ran.append(random.random())
        


        self.LGRB_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran)))
        ran = []
        for i in range(Number_LGRBs_step_1):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        self.LGRB_z = self.LGRB_z * np.array(ran)



        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)




       

        ran = []
        for i in range(self.Number_LGRBs_step_1):
            ran.append(random.random())
        self.LGRB_x_grad = -1 + 2*np.array(ran)

        ran = []
        for i in range(self.Number_LGRBs_step_1):
            ran.append(random.random())

        self.LGRB_y_grad = -1 + 2*np.array(ran)

        ran = []
        for i in range(self.Number_LGRBs_step_1):
            ran.append(random.random())

        self.LGRB_z_grad = -1 + 2*np.array(ran)


        normaliser = (self.LGRB_x_grad**2 + self.LGRB_y_grad**2 + self.LGRB_z_grad**2)**0.5

        LGRB_x_grad = self.LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = self.LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = self.LGRB_z_grad*(1/normaliser)


        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.Number_LGRBs_step_1):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                costheta = np.array(costheta, dtype = float)
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

    def effected_by_LGRB_step_2(self, Number_LGRBs_step_2, SFR_total, E, T):

        self.Number_LGRBs_step_2 = Number_LGRBs_step_2

        ran = []
        for i in range(Number_LGRBs_step_2):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())
        
        rho = []

        for i in range(Number_LGRBs_step_2):
            r = Symbol("r")
            radius = solve((-(17000000*r**3-327000000*r**2+1044960000*r-921755351)/1200000000 - random.random()**(1/2)*5.306 ),cubic = True,complex=False)
            radius = np.array(radius[1])
            radius = re(radius+ (2*I))
            rho.append(radius)

        rho = np.array(rho)


        LGRB_x = rho*np.sin(star_theta)
    
        LGRB_y = rho*np.cos(star_theta)

        z = ((1 - ((LGRB_x /self.stars_x_max)**2)-((LGRB_y /self.stars_y_max)**2))*(self.stars_z_max**2))**(1/2)
 
        opp_adj = rho/z
        opp_adj = np.array(opp_adj, dtype = float)

        star_phi = np.arctan(opp_adj)


        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())
        


        LGRB_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran)))

        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        LGRB_z = LGRB_z * np.array(ran)




        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)


        self.LGRB_x = np.append(self.LGRB_x, LGRB_x)
        self.LGRB_y = np.append(self.LGRB_y, LGRB_y)
        self.LGRB_z = np.append(self.LGRB_z, LGRB_z)



        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())
        self.LGRB_x_grad = np.append(self.LGRB_x_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())

        self.LGRB_y_grad = np.append(self.LGRB_y_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_2):
            ran.append(random.random())

        self.LGRB_z_grad = np.append(self.LGRB_z_grad, -1 + 2*np.array(ran))


        normaliser = (self.LGRB_x_grad**2 + self.LGRB_y_grad**2 + self.LGRB_z_grad**2)**0.5

        LGRB_x_grad = self.LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = self.LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = self.LGRB_z_grad*(1/normaliser)


        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.Number_LGRBs_step_1 + self.Number_LGRBs_step_2):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                costheta = np.array(costheta, dtype = float)
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

    def effected_by_LGRB_step_3(self, Number_LGRBs_step_3, SFR_total, E, T):

        self.Number_LGRBs_step_3 = Number_LGRBs_step_3

        ran = []
        for i in range(Number_LGRBs_step_3):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        rho = []

        for i in range(Number_LGRBs_step_3):
            r = Symbol("r")
            radius = solve((-(17000000*r**3-327000000*r**2+1044960000*r-921755351)/1200000000 -random.random()**(1/2)*5.306 ),cubic = True,complex=False)
            radius = np.array(radius[1])
            radius = re(radius+ (2*I))
            rho.append(radius)

        rho = np.array(rho)


        LGRB_x = rho*np.sin(star_theta)
    
        LGRB_y = rho*np.cos(star_theta)

        z = ((1 - ((LGRB_x /self.stars_x_max)**2)-((LGRB_y /self.stars_y_max)**2))*(self.stars_z_max**2))**(1/2)
 
        opp_adj = rho/z
        opp_adj = np.array(opp_adj, dtype = float)

        star_phi = np.arctan(opp_adj)


        ran = []
        for i in range(self.Number_LGRBs_step_3):
            ran.append(random.random())
        


        LGRB_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran))) 
        ran = []
        for i in range(self.Number_LGRBs_step_3):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        LGRB_z = LGRB_z * np.array(ran)




        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)


        self.LGRB_x = np.append(self.LGRB_x, LGRB_x)
        self.LGRB_y = np.append(self.LGRB_y, LGRB_y)
        self.LGRB_z = np.append(self.LGRB_z, LGRB_z)



        ran = []
        for i in range(self.Number_LGRBs_step_3):
            ran.append(random.random())
        self.LGRB_x_grad = np.append(self.LGRB_x_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_3):
            ran.append(random.random())

        self.LGRB_y_grad = np.append(self.LGRB_y_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_3):
            ran.append(random.random())

        self.LGRB_z_grad = np.append(self.LGRB_z_grad, -1 + 2*np.array(ran))


        normaliser = (self.LGRB_x_grad**2 + self.LGRB_y_grad**2 + self.LGRB_z_grad**2)**0.5

        LGRB_x_grad = self.LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = self.LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = self.LGRB_z_grad*(1/normaliser)


        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.Number_LGRBs_step_1 + self.Number_LGRBs_step_2 + self.Number_LGRBs_step_3):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                costheta = np.array(costheta, dtype = float)               
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


    def effected_by_LGRB_step_4(self, Number_LGRBs_step_4, SFR_total, E, T):

        self.Number_LGRBs_step_4 = Number_LGRBs_step_4

        ran = []
        for i in range(self.Number_LGRBs_step_4):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        rho = []

        for i in range(Number_LGRBs_step_4):
            r = Symbol("r")
            radius = solve((-(17000000*r**3-327000000*r**2+1044960000*r-921755351)/1200000000 - random.random()**(1/2)*5.306 ),cubic = True,complex=False)
            radius = np.array(radius[1])
            radius = re(radius+ (2*I))
            rho.append(radius)

        rho = np.array(rho)


        LGRB_x = rho*np.sin(star_theta)
    
        LGRB_y = rho*np.cos(star_theta)

        z = ((1 - ((LGRB_x /self.stars_x_max)**2)-((LGRB_y /self.stars_y_max)**2))*(self.stars_z_max**2))**(1/2)
 
        opp_adj = rho/z
        opp_adj = np.array(opp_adj, dtype = float)

        star_phi = np.arctan(opp_adj)


        ran = []
        for i in range(self.Number_LGRBs_step_4):
            ran.append(random.random())
        


        LGRB_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran)))

        ran = []
        for i in range(Number_LGRBs_step_4):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        LGRB_z = LGRB_z * np.array(ran)




        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)


        self.LGRB_x = np.append(self.LGRB_x, LGRB_x)
        self.LGRB_y = np.append(self.LGRB_y, LGRB_y)
        self.LGRB_z = np.append(self.LGRB_z, LGRB_z)



        ran = []
        for i in range(self.Number_LGRBs_step_4):
            ran.append(random.random())
        self.LGRB_x_grad = np.append(self.LGRB_x_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_4):
            ran.append(random.random())

        self.LGRB_y_grad = np.append(self.LGRB_y_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_4):
            ran.append(random.random())

        self.LGRB_z_grad = np.append(self.LGRB_z_grad, -1 + 2*np.array(ran))


        normaliser = (self.LGRB_x_grad**2 + self.LGRB_y_grad**2 + self.LGRB_z_grad**2)**0.5

        LGRB_x_grad = self.LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = self.LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = self.LGRB_z_grad*(1/normaliser)


        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.Number_LGRBs_step_1 + self.Number_LGRBs_step_2 + self.Number_LGRBs_step_3 + self.Number_LGRBs_step_4):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                costheta = np.array(costheta, dtype = float)               
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

    def effected_by_LGRB_step_5(self, Number_LGRBs_step_5, SFR_total, E, T):

        self.Number_LGRBs_step_5 = Number_LGRBs_step_5

        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())
        
        star_theta  = np.array(ran)*2*np.pi

        rho = []

        for i in range(Number_LGRBs_step_5):
            r = Symbol("r")
            radius = solve((-(17000000*r**3-327000000*r**2+1044960000*r-921755351)/1200000000 - random.random()**(1/2)*5.306 ),cubic = True,complex=False)
            radius = np.array(radius[1])
            radius = re(radius+ (2*I))
            rho.append(radius)

        rho = np.array(rho)


        LGRB_x = rho*np.sin(star_theta)
    
        LGRB_y = rho*np.cos(star_theta)

        z = ((1 - ((LGRB_x /self.stars_x_max)**2)-((LGRB_y /self.stars_y_max)**2))*(self.stars_z_max**2))**(1/2)
 
        opp_adj = rho/z
        opp_adj = np.array(opp_adj, dtype = float)

        star_phi = np.arctan(opp_adj)


        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())
        


        LGRB_z = -(self.h_c*abs(np.cos(star_phi))) * np.log((1 - np.array(ran))) 

        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())

        for i in range(len(ran)):
            if ran[i] > 0.5:
                ran[i] = -1
            else: 
                ran[i] = 1

        LGRB_z = LGRB_z * np.array(ran)




        effected_star_x = []
        effected_star_y = []
        effected_star_z = []

        theta_max = (6*np.pi/180)


        self.LGRB_x = np.append(self.LGRB_x, LGRB_x)
        self.LGRB_y = np.append(self.LGRB_y, LGRB_y)
        self.LGRB_z = np.append(self.LGRB_z, LGRB_z)



        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())
        self.LGRB_x_grad = np.append(self.LGRB_x_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())

        self.LGRB_y_grad = np.append(self.LGRB_y_grad, -1 + 2*np.array(ran))

        ran = []
        for i in range(self.Number_LGRBs_step_5 ):
            ran.append(random.random())

        self.LGRB_z_grad = np.append(self.LGRB_z_grad, -1 + 2*np.array(ran))


        normaliser = (self.LGRB_x_grad**2 + self.LGRB_y_grad**2 + self.LGRB_z_grad**2)**0.5

        LGRB_x_grad = self.LGRB_x_grad*(1/normaliser)
        LGRB_y_grad = self.LGRB_y_grad*(1/normaliser)
        LGRB_z_grad = self.LGRB_z_grad*(1/normaliser)


        for i in range(len(self.stars_x)):
            


            star_x_check = self.stars_x[i]
            star_y_check = self.stars_y[i]
            star_z_check = self.stars_z[i]

            for j in range(self.Number_LGRBs_step_1 + self.Number_LGRBs_step_2 + self.Number_LGRBs_step_3 + self.Number_LGRBs_step_4 + self.Number_LGRBs_step_5):



                LGRB_Star_vector_x_check   = star_x_check - self.LGRB_x[j]
                LGRB_Star_vector_y_check   = star_y_check - self.LGRB_y[j]
                LGRB_Star_vector_z_check   = star_z_check - self.LGRB_z[j]
                normaliser = (LGRB_Star_vector_x_check**2 + LGRB_Star_vector_y_check**2 + LGRB_Star_vector_z_check**2)**0.5
                    
                LGRB_Star_vector_x_check = (1/normaliser)* LGRB_Star_vector_x_check
                LGRB_Star_vector_y_check = (1/normaliser)* LGRB_Star_vector_y_check
                LGRB_Star_vector_z_check = (1/normaliser)* LGRB_Star_vector_z_check
                costheta = np.abs(LGRB_Star_vector_x_check*LGRB_x_grad[j] + LGRB_Star_vector_y_check*LGRB_y_grad[j] + LGRB_Star_vector_z_check*LGRB_z_grad[j])
                costheta = np.array(costheta, dtype = float)                
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

        self.uneffected_star_x = uneffected_stars_x
        self.uneffected_star_y = uneffected_stars_y
        self.uneffected_star_z = uneffected_stars_z

        if len(uneffected_stars_x) == len(all_stars):
            pass
        else:

            effected_star_x, effected_star_y, effected_star_z = zip(*effected_stars_filtered)

        self.effected_star_x = effected_star_x
        self.effected_star_y = effected_star_y
        self.effected_star_z = effected_star_z




        return(effected_star_x, effected_star_y, effected_star_z, self.LGRB_x, self.LGRB_y, self.LGRB_z, uneffected_stars_x, uneffected_stars_y, uneffected_stars_z, LGRB_x_grad, LGRB_y_grad, LGRB_z_grad)
    def radial_effect(self):

        a_100 = self.stars_x_max
        b_100 = self.stars_y_max
        c_100 = self.stars_z_max

        a_90 = ((3/(np.pi*4))*self.volume * ((90/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_90 = a_90 * (1 + (1-E)**2 - T) ** (1/2)
        c_90 = (1-E)*a_90

        a_80 = ((3/(np.pi*4))*self.volume * ((80/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_80 = a_80 * (1 + (1-E)**2 - T) ** (1/2)
        c_80 = (1-E)*a_80

        a_70 = ((3/(np.pi*4))*self.volume * ((70/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_70 = a_70 * (1 + (1-E)**2 - T) ** (1/2)
        c_70 = (1-E)*a_70

        a_60 = ((3/(np.pi*4))*self.volume * ((60/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_60 = a_60 * (1 + (1-E)**2 - T) ** (1/2)
        c_60 = (1-E)*a_60

        a_50 = ((3/(np.pi*4))*self.volume * ((50/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_50 = a_50 * (1 + (1-E)**2 - T) ** (1/2)
        c_50 = (1-E)*a_50

        a_40 = ((3/(np.pi*4))*self.volume * ((40/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_40 = a_40 * (1 + (1-E)**2 - T) ** (1/2)
        c_40 = (1-E)*a_40

        a_30 = ((3/(np.pi*4))*self.volume * ((30/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_30 = a_30 * (1 + (1-E)**2 - T) ** (1/2)
        c_30 = (1-E)*a_30

        a_20 = ((3/(np.pi*4))*self.volume * ((20/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_20 = a_20 * (1 + (1-E)**2 - T) ** (1/2)
        c_20 = (1-E)*a_20

        a_10 = ((3/(np.pi*4))*self.volume * ((10/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_10 = a_10 * (1 + (1-E)**2 - T) ** (1/2)
        c_10 = (1-E)*a_10

        uneffected_star_x_100_90 = []
        uneffected_star_y_100_90 = []
        uneffected_star_z_100_90 = []

        test_100 = ( (self.uneffected_star_x/a_100)**2 +  (self.uneffected_star_y/b_100)**2  )
        test_90 = ( (self.uneffected_star_x/a_90)**2 +  (self.uneffected_star_y/b_90)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_100[i] <= 1) and (test_90[i]) > 1:
                uneffected_star_x_100_90.append(self.uneffected_star_x[i])
                uneffected_star_y_100_90.append(self.uneffected_star_y[i])
                uneffected_star_z_100_90.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_100_90 = []
        effected_star_y_100_90 = []
        effected_star_z_100_90 = []

        test_100 = ( (self.effected_star_x/a_100)**2 +  (self.effected_star_y/b_100)**2  )
        test_90 = ( (self.effected_star_x/a_90)**2 +  (self.effected_star_y/b_90)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_100[i] <= 1) and (test_90[i]) > 1:
                effected_star_x_100_90.append(self.effected_star_x[i])
                effected_star_y_100_90.append(self.effected_star_y[i])
                effected_star_z_100_90.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_100_90)+ len(effected_star_x_100_90)) == 0:
            self.percentage_100_90 = 0
        else:
            self.percentage_100_90 = 100*(len(effected_star_x_100_90)/(len(uneffected_star_x_100_90)+ len(effected_star_x_100_90)))

        uneffected_star_x_90_80 = []
        uneffected_star_y_90_80 = []
        uneffected_star_z_90_80 = []

        test_90 = ( (self.uneffected_star_x/a_90)**2 +  (self.uneffected_star_y/b_90)**2  )
        test_80 = ( (self.uneffected_star_x/a_80)**2 +  (self.uneffected_star_y/b_80)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_90[i] <= 1) and (test_80[i]) > 1:
                uneffected_star_x_90_80.append(self.uneffected_star_x[i])
                uneffected_star_y_90_80.append(self.uneffected_star_y[i])
                uneffected_star_z_90_80.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_90_80 = []
        effected_star_y_90_80 = []
        effected_star_z_90_80 = []

        test_90 = ( (self.effected_star_x/a_90)**2 +  (self.effected_star_y/b_90)**2 )
        test_80 = ( (self.effected_star_x/a_80)**2 +  (self.effected_star_y/b_80)**2  )
        for i in range(len(self.effected_star_x)):
            if (test_90[i] <= 1) and (test_80[i]) > 1:
                effected_star_x_90_80.append(self.effected_star_x[i])
                effected_star_y_90_80.append(self.effected_star_y[i])
                effected_star_z_90_80.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_90_80)+ len(effected_star_x_90_80)) == 0:
            self.percentage_90_80 = 0
        else:
            self.percentage_90_80 = 100*(len(effected_star_x_90_80)/(len(uneffected_star_x_90_80)+ len(effected_star_x_90_80)))

        uneffected_star_x_80_70 = []
        uneffected_star_y_80_70 = []
        uneffected_star_z_80_70 = []

        test_80 = ( (self.uneffected_star_x/a_80)**2 +  (self.uneffected_star_y/b_80)**2 )
        test_70 = ( (self.uneffected_star_x/a_70)**2 +  (self.uneffected_star_y/b_70)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_80[i] <= 1) and (test_70[i]) > 1:
                uneffected_star_x_80_70.append(self.uneffected_star_x[i])
                uneffected_star_y_80_70.append(self.uneffected_star_y[i])
                uneffected_star_z_80_70.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_80_70 = []
        effected_star_y_80_70 = []
        effected_star_z_80_70 = []

        test_80 = ( (self.effected_star_x/a_80)**2 +  (self.effected_star_y/b_80)**2 )
        test_70 = ( (self.effected_star_x/a_70)**2 +  (self.effected_star_y/b_70)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_80[i] <= 1) and (test_70[i]) > 1:
                effected_star_x_80_70.append(self.effected_star_x[i])
                effected_star_y_80_70.append(self.effected_star_y[i])
                effected_star_z_80_70.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_80_70)+ len(effected_star_x_80_70)) == 0:
            self.percentage_80_70 = 0
        else:
            self.percentage_80_70 = 100*(len(effected_star_x_80_70)/(len(uneffected_star_x_80_70)+ len(effected_star_x_80_70)))

        uneffected_star_x_70_60 = []
        uneffected_star_y_70_60 = []
        uneffected_star_z_70_60 = []

        test_70 = ( (self.uneffected_star_x/a_70)**2 +  (self.uneffected_star_y/b_70)**2 )
        test_60 = ( (self.uneffected_star_x/a_60)**2 +  (self.uneffected_star_y/b_60)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_70[i] <= 1) and (test_60[i]) > 1:
                uneffected_star_x_70_60.append(self.uneffected_star_x[i])
                uneffected_star_y_70_60.append(self.uneffected_star_y[i])
                uneffected_star_z_70_60.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_70_60 = []
        effected_star_y_70_60 = []
        effected_star_z_70_60 = []

        test_70 = ( (self.effected_star_x/a_70)**2 +  (self.effected_star_y/b_70)**2  )
        test_60 = ( (self.effected_star_x/a_60)**2 +  (self.effected_star_y/b_60)**2  )
        for i in range(len(self.effected_star_x)):
            if (test_70[i] <= 1) and (test_60[i]) > 1:
                effected_star_x_70_60.append(self.effected_star_x[i])
                effected_star_y_70_60.append(self.effected_star_y[i])
                effected_star_z_70_60.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_70_60)+ len(effected_star_x_70_60)) == 0:
            self.percentage_70_60 = 0
        else:
            self.percentage_70_60 = 100*(len(effected_star_x_70_60)/(len(uneffected_star_x_70_60)+ len(effected_star_x_70_60)))


        uneffected_star_x_60_50 = []
        uneffected_star_y_60_50 = []
        uneffected_star_z_60_50 = []

        test_60 = ( (self.uneffected_star_x/a_60)**2 +  (self.uneffected_star_y/b_60)**2 )
        test_50 = ( (self.uneffected_star_x/a_50)**2 +  (self.uneffected_star_y/b_50)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_60[i] <= 1) and (test_50[i]) > 1:
                uneffected_star_x_60_50.append(self.uneffected_star_x[i])
                uneffected_star_y_60_50.append(self.uneffected_star_y[i])
                uneffected_star_z_60_50.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_60_50 = []
        effected_star_y_60_50 = []
        effected_star_z_60_50 = []

        test_60 = ( (self.effected_star_x/a_60)**2 +  (self.effected_star_y/b_60)**2  )
        test_50 = ( (self.effected_star_x/a_50)**2 +  (self.effected_star_y/b_50)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_60[i] <= 1) and (test_50[i]) > 1:
                effected_star_x_60_50.append(self.effected_star_x[i])
                effected_star_y_60_50.append(self.effected_star_y[i])
                effected_star_z_60_50.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_60_50)+ len(effected_star_x_60_50)) == 0:
            self.percentage_60_50 = 0
        else:
            self.percentage_60_50 = 100*(len(effected_star_x_60_50)/(len(uneffected_star_x_60_50)+ len(effected_star_x_60_50)))

        uneffected_star_x_50_40 = []
        uneffected_star_y_50_40 = []
        uneffected_star_z_50_40 = []

        test_50 = ( (self.uneffected_star_x/a_50)**2 +  (self.uneffected_star_y/b_50)**2  )
        test_40 = ( (self.uneffected_star_x/a_40)**2 +  (self.uneffected_star_y/b_40)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_50[i] <= 1) and (test_40[i]) > 1:
                uneffected_star_x_50_40.append(self.uneffected_star_x[i])
                uneffected_star_y_50_40.append(self.uneffected_star_y[i])
                uneffected_star_z_50_40.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_50_40 = []
        effected_star_y_50_40 = []
        effected_star_z_50_40 = []

        test_50 = ( (self.effected_star_x/a_50)**2 +  (self.effected_star_y/b_50)**2  )
        test_40 = ( (self.effected_star_x/a_40)**2 +  (self.effected_star_y/b_40)**2  )
        for i in range(len(self.effected_star_x)):
            if (test_50[i] <= 1) and (test_40[i]) > 1:
                effected_star_x_50_40.append(self.effected_star_x[i])
                effected_star_y_50_40.append(self.effected_star_y[i])
                effected_star_z_50_40.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_50_40)+ len(effected_star_x_50_40)) == 0:
            self.percentage_50_40 = 0
        else:
            self.percentage_50_40 = 100*(len(effected_star_x_50_40)/(len(uneffected_star_x_50_40)+ len(effected_star_x_50_40)))

        uneffected_star_x_40_30 = []
        uneffected_star_y_40_30 = []
        uneffected_star_z_40_30 = []

        test_40 = ( (self.uneffected_star_x/a_40)**2 +  (self.uneffected_star_y/b_40)**2  )
        test_30 = ( (self.uneffected_star_x/a_30)**2 +  (self.uneffected_star_y/b_30)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_40[i] <= 1) and (test_30[i]) > 1:
                uneffected_star_x_40_30.append(self.uneffected_star_x[i])
                uneffected_star_y_40_30.append(self.uneffected_star_y[i])
                uneffected_star_z_40_30.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_40_30 = []
        effected_star_y_40_30 = []
        effected_star_z_40_30 = []

        test_40 = ( (self.effected_star_x/a_40)**2 +  (self.effected_star_y/b_40)**2  )
        test_30 = ( (self.effected_star_x/a_30)**2 +  (self.effected_star_y/b_30)**2  )
        for i in range(len(self.effected_star_x)):
            if (test_40[i] <= 1) and (test_30[i]) > 1:
                effected_star_x_40_30.append(self.effected_star_x[i])
                effected_star_y_40_30.append(self.effected_star_y[i])
                effected_star_z_40_30.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_40_30)+ len(effected_star_x_40_30)) == 0:
            self.percentage_40_30 = 0
        else:
            self.percentage_40_30 = 100*(len(effected_star_x_40_30)/(len(uneffected_star_x_40_30)+ len(effected_star_x_40_30)))

        uneffected_star_x_30_20 = []
        uneffected_star_y_30_20 = []
        uneffected_star_z_30_20 = []

        test_30 = ( (self.uneffected_star_x/a_30)**2 +  (self.uneffected_star_y/b_30)**2  )
        test_20 = ( (self.uneffected_star_x/a_20)**2 +  (self.uneffected_star_y/b_20)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_30[i] <= 1) and (test_20[i]) > 1:
                uneffected_star_x_30_20.append(self.uneffected_star_x[i])
                uneffected_star_y_30_20.append(self.uneffected_star_y[i])
                uneffected_star_z_30_20.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_30_20 = []
        effected_star_y_30_20 = []
        effected_star_z_30_20 = []

        test_30 = ( (self.effected_star_x/a_30)**2 +  (self.effected_star_y/b_30)**2  )
        test_20 = ( (self.effected_star_x/a_20)**2 +  (self.effected_star_y/b_20)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_30[i] <= 1) and (test_20[i]) > 1:
                effected_star_x_30_20.append(self.effected_star_x[i])
                effected_star_y_30_20.append(self.effected_star_y[i])
                effected_star_z_30_20.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_30_20)+ len(effected_star_x_30_20)) == 0:
            self.percentage_30_20 = 0
        else:
            self.percentage_30_20 = 100*(len(effected_star_x_30_20)/(len(uneffected_star_x_30_20)+ len(effected_star_x_30_20)))
 
        uneffected_star_x_20_10 = []
        uneffected_star_y_20_10 = []
        uneffected_star_z_20_10 = []

        test_20 = ( (self.uneffected_star_x/a_20)**2 +  (self.uneffected_star_y/b_20)**2  )
        test_10 = ( (self.uneffected_star_x/a_10)**2 +  (self.uneffected_star_y/b_10)**2  )
        for i in range(len(self.uneffected_star_x)):
            if (test_20[i] <= 1) and (test_10[i]) > 1:
                uneffected_star_x_20_10.append(self.uneffected_star_x[i])
                uneffected_star_y_20_10.append(self.uneffected_star_y[i])
                uneffected_star_z_20_10.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_20_10 = []
        effected_star_y_20_10 = []
        effected_star_z_20_10 = []

        test_20 = ( (self.effected_star_x/a_20)**2 +  (self.effected_star_y/b_20)**2  )
        test_10 = ( (self.effected_star_x/a_10)**2 +  (self.effected_star_y/b_10)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_20[i] <= 1) and (test_10[i]) > 1:
                effected_star_x_20_10.append(self.effected_star_x[i])
                effected_star_y_20_10.append(self.effected_star_y[i])
                effected_star_z_20_10.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_20_10)+ len(effected_star_x_20_10)) == 0:
            self.percentage_20_10 = 0
        else:
            self.percentage_20_10 = 100*(len(effected_star_x_20_10)/(len(uneffected_star_x_20_10)+ len(effected_star_x_20_10)))

        uneffected_star_x_10_0 = []
        uneffected_star_y_10_0 = []
        uneffected_star_z_10_0 = []

        test_10 = ( (self.uneffected_star_x/a_10)**2 +  (self.uneffected_star_y/b_10)**2  )

        for i in range(len(self.uneffected_star_x)):
            if (test_10[i] <= 1):
                uneffected_star_x_10_0.append(self.uneffected_star_x[i])
                uneffected_star_y_10_0.append(self.uneffected_star_y[i])
                uneffected_star_z_10_0.append(self.uneffected_star_z[i])
            else:
                pass
    
        effected_star_x_10_0 = []
        effected_star_y_10_0 = []
        effected_star_z_10_0 = []

        test_10 = ( (self.effected_star_x/a_10)**2 +  (self.effected_star_y/b_10)**2 )
        for i in range(len(self.effected_star_x)):
            if (test_10[i] <= 1):
                effected_star_x_10_0.append(self.effected_star_x[i])
                effected_star_y_10_0.append(self.effected_star_y[i])
                effected_star_z_10_0.append(self.effected_star_z[i])
            else:
                pass
        
        if (len(uneffected_star_x_10_0)+ len(effected_star_x_10_0)) == 0:
            self.percentage_10_0 = 0
        else:
            self.percentage_10_0 = 100*(len(effected_star_x_10_0)/(len(uneffected_star_x_10_0)+ len(effected_star_x_10_0)))
    def star_dist(self):

        a_100 = self.stars_x_max
        b_100 = self.stars_y_max
        c_100 = self.stars_z_max

        a_90 = ((3/(np.pi*4))*self.volume * ((90/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_90 = a_90 * (1 + (1-E)**2 - T) ** (1/2)
        c_90 = (1-E)*a_90

        a_80 = ((3/(np.pi*4))*self.volume * ((80/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_80 = a_80 * (1 + (1-E)**2 - T) ** (1/2)
        c_80 = (1-E)*a_80

        a_70 = ((3/(np.pi*4))*self.volume * ((70/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_70 = a_70 * (1 + (1-E)**2 - T) ** (1/2)
        c_70 = (1-E)*a_70

        a_60 = ((3/(np.pi*4))*self.volume * ((60/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_60 = a_60 * (1 + (1-E)**2 - T) ** (1/2)
        c_60 = (1-E)*a_60

        a_50 = ((3/(np.pi*4))*self.volume * ((50/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_50 = a_50 * (1 + (1-E)**2 - T) ** (1/2)
        c_50 = (1-E)*a_50

        a_40 = ((3/(np.pi*4))*self.volume * ((40/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_40 = a_40 * (1 + (1-E)**2 - T) ** (1/2)
        c_40 = (1-E)*a_40

        a_30 = ((3/(np.pi*4))*self.volume * ((30/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_30 = a_30 * (1 + (1-E)**2 - T) ** (1/2)
        c_30 = (1-E)*a_30

        a_20 = ((3/(np.pi*4))*self.volume * ((20/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_20 = a_20 * (1 + (1-E)**2 - T) ** (1/2)
        c_20 = (1-E)*a_20

        a_10 = ((3/(np.pi*4))*self.volume * ((10/100)**3 )* (1/((1-E)*(1-T+T*(1-E)**2)**(1/2))) ) **(1/3)
        b_10 = a_10 * (1 + (1-E)**2 - T) ** (1/2)
        c_10 = (1-E)*a_10

        stars_x_100_90 = []
        stars_y_100_90 = []
        stars_z_100_90 = []

        test_100 = ( (self.stars_x/a_100)**2 +  (self.stars_y/b_100)**2  )
        test_90 = ( (self.stars_x/a_90)**2 +  (self.stars_y/b_90)**2  )
        for i in range(len(self.stars_x)):
            if (test_100[i] <= 1) and (test_90[i]) > 1:
                stars_x_100_90.append(self.stars_x[i])
                stars_y_100_90.append(self.stars_y[i])
                stars_z_100_90.append(self.stars_z[i])
            else:
                pass
    

        self.percentage_100_90 = 100*(len(stars_x_100_90)/(len(self.stars_x)))

        stars_x_90_80 = []
        stars_y_90_80 = []
        stars_z_90_80 = []

        test_90 = ( (self.stars_x/a_90)**2 +  (self.stars_y/b_90)**2  )
        test_80 = ( (self.stars_x/a_80)**2 +  (self.stars_y/b_80)**2  )
        for i in range(len(self.stars_x)):
            if (test_90[i] <= 1) and (test_80[i]) > 1:
                stars_x_90_80.append(self.stars_x[i])
                stars_y_90_80.append(self.stars_y[i])
                stars_z_90_80.append(self.stars_z[i])
            else:
                pass
    
    

        self.percentage_90_80 = 100*(len(stars_x_90_80)/(len(self.stars_x)))

        stars_x_80_70 = []
        stars_y_80_70 = []
        stars_z_80_70 = []

        test_80 = ( (self.stars_x/a_80)**2 +  (self.stars_y/b_80)**2 )
        test_70 = ( (self.stars_x/a_70)**2 +  (self.stars_y/b_70)**2  )
        for i in range(len(self.stars_x)):
            if (test_80[i] <= 1) and (test_70[i]) > 1:
                stars_x_80_70.append(self.stars_x[i])
                stars_y_80_70.append(self.stars_y[i])
                stars_z_80_70.append(self.stars_z[i])
            else:
                pass
 
        self.percentage_80_70 = 100*(len(stars_x_80_70)/(len(self.stars_x)))
        stars_x_70_60 = []
        stars_y_70_60 = []
        stars_z_70_60 = []

        test_70 = ( (self.stars_x/a_70)**2 +  (self.stars_y/b_70)**2 )
        test_60 = ( (self.stars_x/a_60)**2 +  (self.stars_y/b_60)**2  )
        for i in range(len(self.stars_x)):
            if (test_70[i] <= 1) and (test_60[i]) > 1:
                stars_x_70_60.append(self.stars_x[i])
                stars_y_70_60.append(self.stars_y[i])
                stars_z_70_60.append(self.stars_z[i])
            else:
                pass


        self.percentage_70_60 = 100*(len(stars_x_70_60)/(len(self.stars_x)))



        stars_x_60_50 = []
        stars_y_60_50 = []
        stars_z_60_50 = []

        test_60 = ( (self.stars_x/a_60)**2 +  (self.stars_y/b_60)**2 )
        test_50 = ( (self.stars_x/a_50)**2 +  (self.stars_y/b_50)**2  )
        for i in range(len(self.stars_x)):
            if (test_60[i] <= 1) and (test_50[i]) > 1:
                stars_x_60_50.append(self.stars_x[i])
                stars_y_60_50.append(self.stars_y[i])
                stars_z_60_50.append(self.stars_z[i])
            else:
                pass
    
        self.percentage_60_50 = 100*(len(stars_x_60_50)/(len(self.stars_x)))


        stars_x_50_40 = []
        stars_y_50_40 = []
        stars_z_50_40 = []

        test_50 = ( (self.stars_x/a_50)**2 +  (self.stars_y/b_50)**2  )
        test_40 = ( (self.stars_x/a_40)**2 +  (self.stars_y/b_40)**2  )
        for i in range(len(self.stars_x)):
            if (test_50[i] <= 1) and (test_40[i]) > 1:
                stars_x_50_40.append(self.stars_x[i])
                stars_y_50_40.append(self.stars_y[i])
                stars_z_50_40.append(self.stars_z[i])
            else:
                pass
    
        self.percentage_50_40 = 100*(len(stars_x_50_40)/(len(self.stars_x)))


        stars_x_40_30 = []
        stars_y_40_30 = []
        stars_z_40_30 = []

        test_40 = ( (self.stars_x/a_40)**2 +  (self.stars_y/b_40)**2  )
        test_30 = ( (self.stars_x/a_30)**2 +  (self.stars_y/b_30)**2  )
        for i in range(len(self.stars_x)):
            if (test_40[i] <= 1) and (test_30[i]) > 1:
                stars_x_40_30.append(self.stars_x[i])
                stars_y_40_30.append(self.stars_y[i])
                stars_z_40_30.append(self.stars_z[i])
            else:
                pass
    
        self.percentage_40_30 = 100*(len(stars_x_40_30)/(len(self.stars_x)))


        stars_x_30_20 = []
        stars_y_30_20 = []
        stars_z_30_20 = []

        test_30 = ( (self.stars_x/a_30)**2 +  (self.stars_y/b_30)**2  )
        test_20 = ( (self.stars_x/a_20)**2 +  (self.stars_y/b_20)**2  )
        for i in range(len(self.stars_x)):
            if (test_30[i] <= 1) and (test_20[i]) > 1:
                stars_x_30_20.append(self.stars_x[i])
                stars_y_30_20.append(self.stars_y[i])
                stars_z_30_20.append(self.stars_z[i])
            else:
                pass
    
        self.percentage_30_20 = 100*(len(stars_x_30_20)/(len(self.stars_x)))

 
        stars_x_20_10 = []
        stars_y_20_10 = []
        stars_z_20_10 = []

        test_20 = ( (self.stars_x/a_20)**2 +  (self.stars_y/b_20)**2  )
        test_10 = ( (self.stars_x/a_10)**2 +  (self.stars_y/b_10)**2  )
        for i in range(len(self.stars_x)):
            if (test_20[i] <= 1) and (test_10[i]) > 1:
                stars_x_20_10.append(self.stars_x[i])
                stars_y_20_10.append(self.stars_y[i])
                stars_z_20_10.append(self.stars_z[i])
            else:
                pass
    
        self.percentage_20_10 = 100*(len(stars_x_20_10)/(len(self.stars_x)))


        stars_x_10_0 = []
        stars_y_10_0 = []
        stars_z_10_0 = []

        test_10 = ( (self.stars_x/a_10)**2 +  (self.stars_y/b_10)**2  )

        for i in range(len(self.stars_x)):
            if (test_10[i] <= 1):
                stars_x_10_0.append(self.stars_x[i])
                stars_y_10_0.append(self.stars_y[i])
                stars_z_10_0.append(self.stars_z[i])
            else:
                pass
    

            self.percentage_10_0 = 100*(len(stars_x_10_0)/(len(self.stars_x)))
mySims = np.array([('RefL0100N1504', 100.)])
con = sql.connect("jwx596", password="SPS530AQ")













for sim_name , sim_size in mySims:
    print (sim_name)
    myQuery = "SELECT \
                    SH.HalfMassRad_star as r , \
                    MK.Ellipticity as E , \
                    MK.Triaxiality as T \
                FROM \
                    %s_SubHalo as SH, \
                    %s_MorphoKinem as MK, \
                    %s_Aperture as AP \
                WHERE \
                    MK.GalaxyID = AP.GalaxyID  and \
                    SH.GalaxyID = AP.GalaxyID and \
                    SH.GalaxyID = 17007154 and \
                    AP.ApertureSize = 100"%(sim_name, sim_name, sim_name)

myData = sql.execute_query(con, myQuery)

Half_mass_radius = 0.67
T = 0
E = 0.92

a_halfmass = (Half_mass_radius ** 3 / ((1-E)*(1 - T + T*(1 - E)**2)**(1/2)))**(1/3)
b_halfmass = a_halfmass * (1 + (1-E)**2 - T) ** (1/2)
c_halfmass = (1-E)*a_halfmass


h_a = 2.15
h_b = 2.15
h_c = 2.15*0.08

##Now we need the historic SFR to find the total number of LGRBs over the last 5Gyrs##

#But first we must find its main branch progenitor IDs###

for sim_name , sim_size in mySims:
    print (sim_name)
    myQuery = "SELECT  \
                    SH.GalaxyID as ID, \
                    SH.Redshift as z, \
                    SH.Stars_Metallicity as Z,  \
                    AP.SFR as SFR \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Subhalo as REF , \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.Aperturesize = 100 and \
                    REF.GalaxyID =17007154 and  \
                    (( SH.SnapNum > REF.SnapNum and REF.GalaxyID  \
                    between SH.GalaxyID and SH.TopLeafID ) or  \
                    (SH.SnapNum <= REF.SnapNum and SH.GalaxyID  \
                    between REF.GalaxyID and REF.TopLeafID )) \
                ORDER BY \
                    SH.Redshift"%(sim_name, sim_name, sim_name)

myData = sql.execute_query(con, myQuery)

#Number of LGRBs to have occured in galaxies we are interested in 
SFR = 2

LGRBrate = (0.007*0.22*0.69*0.21*0.046)*SFR*0.1



#low metallicty fraction for up to redshift = 4




f_z = 0.1



galaxy = Galaxy(h_a, h_b, h_c, 1000)


galaxy.star_dist()

cumulative = [galaxy.percentage_10_0, galaxy.percentage_10_0+galaxy.percentage_20_10, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30,  galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40,     galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40+galaxy.percentage_60_50, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40+galaxy.percentage_60_50+galaxy.percentage_70_60, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40+galaxy.percentage_60_50+galaxy.percentage_70_60+galaxy.percentage_80_70, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40+galaxy.percentage_60_50+galaxy.percentage_70_60+galaxy.percentage_80_70 +galaxy.percentage_90_80, galaxy.percentage_10_0+galaxy.percentage_20_10+galaxy.percentage_30_20+galaxy.percentage_40_30 +galaxy.percentage_50_40+galaxy.percentage_60_50+galaxy.percentage_70_60+galaxy.percentage_80_70 +galaxy.percentage_90_80+galaxy.percentage_100_90                                   ]
percentile = [5*(galaxy.stars_y_max/100), 15*(galaxy.stars_y_max/100), 25*(galaxy.stars_y_max/100), 35*(galaxy.stars_y_max/100), 45*(galaxy.stars_y_max/100), 55*(galaxy.stars_y_max/100), 65*(galaxy.stars_y_max/100), 75*(galaxy.stars_y_max/100), 85*(galaxy.stars_y_max/100), 95*(galaxy.stars_y_max/100)]

plt.plot( percentile,cumulative, linewidth = 0.8)
plt.show()


myQuery_17007159 = "SELECT  \
                    SH.GalaxyID as ID, \
                    AP.SFR as SFR , \
                    AP.ApertureSize as Size \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.ApertureSize between 10 and 100 and \
                    SH.GalaxyID =17007159 \
                ORDER BY \
                    AP.ApertureSize"%(sim_name, sim_name)

myData_17007159 = sql.execute_query(con, myQuery_17007159)

Number_LGRBs_17007159 = int (SFR * 1.1e9 * (0.007*0.22*0.69*0.21*0.046)*f_z) 
print(Number_LGRBs_17007159)


SFR_total =  myData_17007159['SFR'][6]


SFR_100_70 = myData_17007159['SFR'][6] - myData_17007159['SFR'][5]
SFR_70_50 = myData_17007159['SFR'][5] - myData_17007159['SFR'][4]
SFR_50_40 = myData_17007159['SFR'][4] - myData_17007159['SFR'][3]
SFR_40_30 = myData_17007159['SFR'][3] - myData_17007159['SFR'][2]
SFR_30_20 = myData_17007159['SFR'][2] - myData_17007159['SFR'][1]
SFR_20_10 = myData_17007159['SFR'][1] - myData_17007159['SFR'][0]
SFR_10_0 = myData_17007159['SFR'][0]



galaxy_effected_step_1 = galaxy.effected_by_LGRB_step_1(Number_LGRBs_17007159,SFR_total,E,T)
percentage_effected_step_1 = len(galaxy_effected_step_1[0])/(len(galaxy_effected_step_1[0])+len(galaxy_effected_step_1[6]))*100
print(percentage_effected_step_1)


fig = plt.figure()

ax1 = fig.add_subplot(111, projection = '3d')

ax1.set_zlim3d(-15, 15)
ax1.set_ylim3d(-15, 15)
ax1.set_xlim3d(-15, 15)
ax1.scatter(galaxy_effected_step_1[6], galaxy_effected_step_1[7], galaxy_effected_step_1[8], marker='o')
ax1.scatter(galaxy_effected_step_1[0], galaxy_effected_step_1[1], galaxy_effected_step_1[2], marker='o', color ='r')
ax1.scatter(galaxy_effected_step_1[3], galaxy_effected_step_1[4], galaxy_effected_step_1[5], marker='o', color ='g')



ax1.set_xlabel("$x - Radius$ [$KPc$]")
ax1.set_ylabel("$y - Radius$ [$KPc$]")
ax1.set_zlabel("$z - Width$ [$KPc$]")

plt.show()


myQuery_17007158 = "SELECT  \
                    SH.GalaxyID as ID, \
                    AP.SFR as SFR , \
                    AP.ApertureSize as Size \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.ApertureSize between 10 and 100 and \
                    SH.GalaxyID =17007158 \
                ORDER BY \
                    AP.ApertureSize"%(sim_name, sim_name)

myData_17007158 = sql.execute_query(con, myQuery_17007158)

Number_LGRBs_17007158 = int (SFR* 0.88e9 * (0.007*0.22*0.69*0.21*0.046)*f_z)

print(Number_LGRBs_17007158)


SFR_total =  myData_17007158['SFR'][6]


SFR_100_70 = myData_17007158['SFR'][6] - myData_17007158['SFR'][5]
SFR_70_50 = myData_17007158['SFR'][5] - myData_17007158['SFR'][4]
SFR_50_40 = myData_17007158['SFR'][4] - myData_17007158['SFR'][3]
SFR_40_30 = myData_17007158['SFR'][3] - myData_17007158['SFR'][2]
SFR_30_20 = myData_17007158['SFR'][2] - myData_17007158['SFR'][1]
SFR_20_10 = myData_17007158['SFR'][1] - myData_17007158['SFR'][0]
SFR_10_0 = myData_17007158['SFR'][0]


galaxy_effected_step_2 = galaxy.effected_by_LGRB_step_2(Number_LGRBs_17007158,SFR_total,E,T)
percentage_effected__step_2 = len(galaxy_effected_step_2[0])/(len(galaxy_effected_step_2[0])+len(galaxy_effected_step_2[6]))*100
print(percentage_effected__step_2)



myQuery_17007157 = "SELECT  \
                    SH.GalaxyID as ID, \
                    AP.SFR as SFR , \
                    AP.ApertureSize as Size \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.ApertureSize between 10 and 100 and \
                    SH.GalaxyID =17007157 \
                ORDER BY \
                    AP.ApertureSize"%(sim_name, sim_name)

myData_17007157 = sql.execute_query(con, myQuery_17007157)

Number_LGRBs_17007157  = int (SFR * 0.92e9 * (0.007*0.22*0.69*0.21*0.046)*f_z)

print(Number_LGRBs_17007157)



SFR_total =  myData_17007157['SFR'][6]


SFR_100_70 = myData_17007157['SFR'][6] - myData_17007157['SFR'][5]
SFR_70_50 = myData_17007157['SFR'][5] - myData_17007157['SFR'][4]
SFR_50_40 = myData_17007157['SFR'][4] - myData_17007157['SFR'][3]
SFR_40_30 = myData_17007157['SFR'][3] - myData_17007157['SFR'][2]
SFR_30_20 = myData_17007157['SFR'][2] - myData_17007157['SFR'][1]
SFR_20_10 = myData_17007157['SFR'][1] - myData_17007157['SFR'][0]
SFR_10_0 = myData_17007157['SFR'][0]


galaxy_effected_step_3 = galaxy.effected_by_LGRB_step_3(Number_LGRBs_17007157,SFR_total,E,T)
percentage_effected__step_3 = len(galaxy_effected_step_3[0])/(len(galaxy_effected_step_3[0])+len(galaxy_effected_step_3[6]))*100
print(percentage_effected__step_3)



myQuery_17007156 = "SELECT  \
                    SH.GalaxyID as ID, \
                    AP.SFR as SFR , \
                    AP.ApertureSize as Size \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.ApertureSize between 10 and 100 and \
                    SH.GalaxyID =17007156 \
                ORDER BY \
                    AP.ApertureSize"%(sim_name, sim_name)

myData_17007156 = sql.execute_query(con, myQuery_17007156)

Number_LGRBs_17007156 =  int (SFR * 0.97e9 * (0.007*0.22*0.69*0.21*0.046)*f_z)

print(Number_LGRBs_17007156)


SFR_total =  myData_17007156['SFR'][6] 


SFR_100_70 = myData_17007156['SFR'][6] - myData_17007156['SFR'][5]
SFR_70_50 = myData_17007156['SFR'][5] - myData_17007156['SFR'][4]
SFR_50_40 = myData_17007156['SFR'][4] - myData_17007156['SFR'][3]
SFR_40_30 = myData_17007156['SFR'][3] - myData_17007156['SFR'][2]
SFR_30_20 = myData_17007156['SFR'][2] - myData_17007156['SFR'][1]
SFR_20_10 = myData_17007156['SFR'][1] - myData_17007156['SFR'][0]
SFR_10_0 = myData_17007156['SFR'][0]


galaxy_effected_step_4 = galaxy.effected_by_LGRB_step_4(Number_LGRBs_17007156,SFR_total,E,T)
percentage_effected__step_4 = len(galaxy_effected_step_4[0])/(len(galaxy_effected_step_4[0])+len(galaxy_effected_step_4[6]))*100
print(percentage_effected__step_4)



myQuery_17007155 = "SELECT  \
                    SH.GalaxyID as ID, \
                    AP.SFR as SFR , \
                    AP.ApertureSize as Size \
                FROM \
                    %s_Subhalo as SH ,  \
                    %s_Aperture as AP  \
                WHERE \
                    AP.GalaxyID = SH.GalaxyID  and  \
                    AP.ApertureSize between 10 and 100 and \
                    SH.GalaxyID =17007155 \
                ORDER BY \
                    AP.ApertureSize"%(sim_name, sim_name)

myData_17007155 = sql.execute_query(con, myQuery_17007155)

Number_LGRBs_17007155 =  int (SFR * 1.35e9 * (0.007*0.22*0.69*0.21*0.046)*f_z)

print(Number_LGRBs_17007155)

SFR_total =  myData_17007155['SFR'][6]


SFR_100_70 = myData_17007155['SFR'][6] - myData_17007155['SFR'][5]
SFR_70_50 = myData_17007155['SFR'][5] - myData_17007155['SFR'][4]
SFR_50_40 = myData_17007155['SFR'][4] - myData_17007155['SFR'][3]
SFR_40_30 = myData_17007155['SFR'][3] - myData_17007155['SFR'][2]
SFR_30_20 = myData_17007155['SFR'][2] - myData_17007155['SFR'][1]
SFR_20_10 = myData_17007155['SFR'][1] - myData_17007155['SFR'][0]
SFR_10_0 = myData_17007155['SFR'][0]


galaxy_effected_step_5 = galaxy.effected_by_LGRB_step_5(Number_LGRBs_17007155,SFR_total,E,T)
percentage_effected__step_5 = len(galaxy_effected_step_5[0])/(len(galaxy_effected_step_5[0])+len(galaxy_effected_step_5[6]))*100
print(percentage_effected__step_5)





galaxy.radial_effect()

percentage_by_percentile = [galaxy.percentage_10_0, galaxy.percentage_20_10, galaxy.percentage_30_20, galaxy.percentage_40_30, galaxy.percentage_50_40, galaxy.percentage_60_50, galaxy.percentage_70_60, galaxy.percentage_80_70, galaxy.percentage_90_80, galaxy.percentage_100_90 ]
percentile = [5*(galaxy.stars_y_max/100), 15*(galaxy.stars_y_max/100), 25*(galaxy.stars_y_max/100), 35*(galaxy.stars_y_max/100), 45*(galaxy.stars_y_max/100), 55*(galaxy.stars_y_max/100), 65*(galaxy.stars_y_max/100), 75*(galaxy.stars_y_max/100), 85*(galaxy.stars_y_max/100), 95*(galaxy.stars_y_max/100)]





plt.plot(percentile, percentage_by_percentile, linewidth = 0.8, color='b')



plt.title('Radial probability of extinction')
plt.xlabel(r'Radius [KPc]', fontsize=20)
plt.ylabel(r'Percentage of stars effected', fontsize=20)
plt.tight_layout()
plt.savefig('MW Radial 2.png')

plt.show()



















