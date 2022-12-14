import sys
import os
sys.path.append("../../PsrPopPy2-master/lib/python")
import numpy as np, astropy, scipy
import os
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.interpolate as pol

import dosurvey as dosurvey, populate

#Create a number of populations and see how many 1757-like binaries we would detect

def combine_spline(beta,corr_func_MP,corr_func_IP):
    if beta < 10.:
        h = corr_func_MP(beta)
    
    else:
        h = max(corr_func_MP(beta),corr_func_IP(beta))
    
    return h

pulsar=sys.argv[1]

#default intensity profile correction
corr_func_MP = 1
corr_func_IP = 1
corr_norm = 1

if pulsar=="0509":
    log_p = np.log10(76.5)
    pulse_duty = 18.
    model_name = './0509_z330.model'
    print "On 0509"

if pulsar=="1102":
    log_p = np.log10(27.28500685)
    pulse_duty = 6.
    model_name = './1102_z330.model'
    print "On 1102"

if pulsar=="1946":
    log_p = 1.2278867046136734
    pulse_duty = 6.
    model_name = './1946_z330.model'
    print "On 1946'"

if pulsar=="0737A":
    log_p = 1.3560258571931227
    pulse_duty = 27.
    model_name = './0737A_z330.model'
    print "On 0737A"

if pulsar=="1757":
    log_p = 1.3324384599156054
    pulse_duty = 6.
    model_name = './1757_z330.model'
    print "On 1757"

if pulsar=="1534":
    log_p = 1.5786392099680724
    pulse_duty = 4.
    model_name = './1534_z330.model'
    print "On 1534"

if pulsar=="1913":
    log_p = 1.7709992051639407
    pulse_duty = 16.9
    model_name = './1913_z330.model'
    print "On 1913"

if pulsar=="0737B":
    log_p = 3.4471580313422194
    pulse_duty = 1.3
    model_name = './0737B_z330.model'
    print "On 0737B"

if pulsar=="1756":
    log_p = np.log10(28.462)
    pulse_duty = 3.
    model_name = './1756_2251_z330.model'
    print "On 1756"

if pulsar=="1906":

    log_p = np.log10(144.073)
    pulse_duty = 1.
    model_name = './1906_0746_spline.model'
    print "On 1906+0746"

    #new intensity profile correction
    profile_MP = np.loadtxt("mainpulse_mirror_log.txt").transpose()  #load determined profiles
    profile_IP = np.loadtxt("interpulse_mirror_log.txt").transpose()

    t_MP = np.append([-20.,-19,-18.,-17.,-16.,], np.append(profile_MP[0][36:58],[16.,17.,18.,19.,20.]))  #internal knots MP
    t_IP = np.append(profile_IP[0][2:6], np.append([-5.1,-5.,-4.,-3.,0.,3.,4.,5.,5.1], profile_IP[0][88:92])) + 18.82    #internal knots IP

    corr_func_MP = pol.LSQUnivariateSpline(profile_MP[0], profile_MP[1], t_MP)
    corr_func_IP = pol.LSQUnivariateSpline(profile_IP[0]+18.82, profile_IP[1], t_IP)
 
    corr_norm = max(corr_func_MP(np.arange(-12.,12.,0.001)))

if pulsar=="1208":
    log_p = np.log10(28.714)
    pulse_duty = 8.
    model_name = './1208_z330.model'
    print "On 1208"
    
#survey_dict = {'1946': 'PALFA_one_v_older', '0737A': 'PHSURV', '1757': 'HTRU_low_1757', '1534': '1534_survey', '1913': '1913_survey', '0737B': 'PHSURV', '1102': 'PALFA_one_v_older'}

surveys = ['PALFA_one_v_older', 'PHSURV', 'HTRU_low_1757', '1534_survey', 'PMSURV', 'GBNCC', 'MMGPS-L']
#surveys = ['MMGPS-L']
pops = np.arange(500, 4100, 100)   #default: 10,6000,100
runs_of_pop = 100  # type: int, default 100
runs_per_pop = 100

#lum_params = [-1.5, 0.94]     #From Bagchi, Lorimer, Jayanth, 2011
lum_params = [-1.1, 0.9]       #Fiducial parameters
    
detections = np.full((runs_of_pop, runs_per_pop), 0)
    
os.system("mkdir ./pop_result_"+pulsar+"_general/")

for xx, npop in enumerate(pops):

    filepath = "./pop_result_"+pulsar+"_general/npop={:5d}.npy".format(npop)  # check if the file already exists
    
    # print "On population:{}".format(npop)
    if os.path.exists(filepath):
        print "Pfad {}.npy existiert".format(npop)
        continue
    
    for ii in range(runs_of_pop):
    
        pop = populate.generate(npop, pulsar, corr_func=combine_spline, corr_func_MP=corr_func_MP, corr_func_IP=corr_func_IP, corr_norm=corr_norm, pDistPars = [log_p, 0.0], duty_percent = pulse_duty, orbits = False, nostdout = True, zscale = 0.33, siDistPars=[-1.4, 1.0], lumDistPars = lum_params)

        pop.write(outf = model_name)
    
        for jj in range(runs_per_pop):
                
            population = dosurvey.loadModel(model_name)
            
            output = dosurvey.run(population, surveys, pulsar, corr_func=combine_spline, corr_func_MP=corr_func_MP, corr_func_IP=corr_func_IP, corr_norm=corr_norm, nostdout = True)

            detections[ii][jj] = output[0][2].ndet + output[1][2].ndet + output[2][2].ndet + output[3][2].ndet + output[4][2].ndet + output[5][2].ndet + output[6][2].ndet
#            detections[ii][jj] = output[0][2].ndet

    #filepath = "./pop_result_z330_1906_0746/npop={:5d}.npy".format(npop)
    np.save(filepath, detections)
    
print pulsar+" all Done!"

