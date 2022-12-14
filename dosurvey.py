#!/usr/bin/python

import sys
import argparse
import math
import random
from numpy.random import rand
import numpy as np
from scipy.special import erf

import cPickle

from population import Population
from pulsar import Pulsar
from survey import Survey


class Detections:
    """Just a simple object to store survey detection summary"""
    def __init__(self,
                 ndet=None,
                 ndisc=None,
                 nsmear=None,
                 nout=None,
                 nbr=None,
                 ntf=None):
        self.ndet = ndet
        self.ndisc = ndisc
        self.nsmear = nsmear
        self.nout = nout
        self.nbr = nbr
        self.nfaint = ntf


def loadModel(popfile='populate.model', popmodel=None):
    """Loads in either a model from disk (popfile, cPickle),
       or pass in a model from memory (popmodel)"""
    if popmodel is None:
        with open(popfile, 'rb') as f:
            pop = cPickle.load(f)
    else:
        pop = popmodel

    return pop


def write(surveyPops,
          extension='.results',
          nores=False,
          asc=False,
          summary=False):
    """Write a survey results population to a binary file."""

    for surv, survpop, detected in surveyPops:
        # create an output file, if required
        if not nores:
            if surv is not None:
                s = ''.join([surv, extension])

                survpop.write(s)
            else:
                s = 'allsurveys.results'
                survpop.write(s)

        # Write ascii file if required
        if asc and surv is not None:
            survpop.write_asc(surv + '.det')

        if summary and surv is not None:
            # Write a summary file for the survey (if true)
            filename = ''.join([surv, '.summary'])
            s = 'Detected {0}'.format(detected.ndet)
            s = '\n'.join([s, 'Ndiscovered {0}'.format(detected.ndisc)])
            s = '\n'.join([s, 'Nsmear {0}'.format(detected.nsmear)])
            s = '\n'.join([s, 'Nfaint {0}'.format(detected.nfaint)])
            s = '\n'.join([s, 'Nout {0}'.format(detected.nout)])
            s = '\n'.join([s, 'Nbr {0}'.format(detected.nbr)])
            s = ''.join([s, '\n'])

            with open(filename, 'w') as output:
                output.write(s)

def bimodal_gaussian(x,sigma1,x1,A1,sigma2,x2,A2):

    gaussian1 = A1 * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2))
    gaussian2 = A2 * np.exp(-(x - x2) ** 2 / (2 * sigma2 ** 2))

    return gaussian1 + gaussian2

def skewed_gaussian(x,sigma,x0,A,sk):

	gaussian = A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
	skew = (1+erf(sk*(x-x0)/(sigma*np.sqrt(2))))

	return gaussian*skew

def bimodal_skewed_gaussian(x,sigma1,x1,A1,sk1,sigma2,x2,A2,sk2):

    gaussian1 = A1 * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2))
    gaussian2 = A2 * np.exp(-(x - x2) ** 2 / (2 * sigma2 ** 2))
    skew1 = (1+erf(sk1*(x-x1)/(sigma1*np.sqrt(2))))
    skew2 = (1+erf(sk2*(x-x2)/(sigma2*np.sqrt(2))))

    return gaussian1*skew1 + gaussian2*skew2

def random_beam_dimensions():

    got_major=False

    while got_major==False:
        major_x_rand=rand()*45+20
        major_y_rand=rand()*65
        major_y=skewed_gaussian(major_x_rand,9.004263903263585,30.185883926792673,35.92039069270681,4.221258645066683)
        if major_y_rand<major_y:
            major=major_y_rand
            got_major=True

    got_minor=False

    while got_minor==False:
        minor_x_rand=rand()*15+20
        minor_y_rand=rand()*252
        minor_y=bimodal_skewed_gaussian(minor_x_rand,1.5813541343296091,23.19434345354499,131.6197129070408,8.523764790070745,1.3592721993829742,26.502940358203706,91.68825508633608,5.631889059801915)
        if minor_y_rand<minor_y:
            minor=minor_y_rand
            got_minor=True

    return major,minor #arcsec

def random_tiling_dimensions():

    got_tiling=False

    while got_tiling==False:
        tiling_x_rand=rand()*17+4
        tiling_y_rand=rand()*108
        tiling_y=bimodal_gaussian(tiling_x_rand,0.8982772232058504,12.5459929791275,80.58422913080236,2.263747678383036,12.760781687193873,27.395069356103154)
        if tiling_y_rand<tiling_y:
            tiling=tiling_y_rand
            got_tiling=True

    return tiling #arcmin

def Gaussian2D_degradation(x, y, xSigma, ySigma):

    return np.exp(-x**2/(2*xSigma**2))*np.exp(-y**2/(2*ySigma**2))

def run(pop,
        surveyList,
        psr_name,
        corr_func=1,
        corr_func_MP=1,
        corr_func_IP=1,
        corr_norm=1,
        nostdout=False,
        allsurveyfile=False,
        scint=False,
        accelsearch=False,
        jerksearch=False,
        rratssearch=False):

    """ Run the surveys and detect the pulsars."""

    # print the population
    if not nostdout:
        print "Running doSurvey on population..."
        print pop

    # loop over the surveys we want to run on the pop file
    surveyPops = []

    for surv in surveyList:

        s = Survey(surv, psr_name, corr_func, corr_func_MP, corr_func_IP, corr_norm)
        s.discoveries = 0
        if not nostdout:
            print "\nRunning survey {0}".format(surv)

        # create a new population object to store discovered pulsars in
        survpop = Population()
        # HERE SHOULD INCLUDE THE PROPERTIES OF THE ORIGINAL POPULATION

        # counters
        nsmear = 0
        nout = 0
        ntf = 0
        ndet = 0
        nbr = 0
        # loop over the pulsars in the population list
        for psr in pop.population:
            # pulsar could be dead (evolve!) - continue if so
            if psr.dead:
                continue

            # is the pulsar over the detection threshold?
            (snr, beam_deg_factor, offset) = s.SNRcalc(psr,pop,accelsearch,jerksearch,rratssearch)

            #print snr
            ######################################################################################
            #This section is to add in degradation due to orbital motion of DNS pulsars.
            #Remove this if doing literally anything else.
            #Right now this is very ad-hoc and manual. Optimization possible, maybe not worth right now.
            deg_fac_0509 = {'PALFA_one_v_older': 0.9999, 'PHSURV': 0.9999, 'HTRU_low_1757': 0.9983, '1534_survey': 0.9999, 'PMSURV': 0.8678, 'GBNCC': 0.9999, 'MMGPS-L': 0.9983 }
            deg_fac_1102 = {'PALFA_one_v_older': 0.9997, 'PHSURV': 0.9997, 'HTRU_low_1757': 0.9912, '1534_survey': 0.9999, 'PMSURV': 0.4649, 'GBNCC': 0.9999, 'MMGPS-L': 0.9913 } #This is for 1913+1102
            deg_fac_1913 = {'PALFA_one_v_older': 0.9953, 'PHSURV': 0.9956, 'HTRU_low_1757': 0.9569, '1534_survey': 0.9999, 'PMSURV': 0.7915, 'GBNCC': 0.9999, 'MMGPS-L': 0.9936 }
            deg_fac_1946 = {'PALFA_one_v_older': 0.9514, 'PHSURV': 0.9543, 'HTRU_low_1757': 0.6513, '1534_survey': 0.9999, 'PMSURV': 0.2368, 'GBNCC': 0.9995, 'MMGPS-L': 0.6514 }
            deg_fac_1757 = {'PALFA_one_v_older': 0.9710, 'PHSURV': 0.9716, 'HTRU_low_1757': 0.9255, '1534_survey': 0.9999, 'PMSURV': 0.5080, 'GBNCC': 0.9995, 'MMGPS-L': 0.9255 }
            deg_fac_1534 = {'PALFA_one_v_older': 0.9999, 'PHSURV': 0.9999, 'HTRU_low_1757': 0.9994, '1534_survey': 0.9999, 'PMSURV': 0.7759, 'GBNCC': 0.9999, 'MMGPS-L': 0.9995 }
            deg_fac_0737A= {'PALFA_one_v_older': 0.9910, 'PHSURV': 0.9916, 'HTRU_low_1757': 0.8371, '1534_survey': 0.9999, 'PMSURV': 0.2991, 'GBNCC': 0.9999, 'MMGPS-L': 0.8371 }
            deg_fac_1756 = {'PALFA_one_v_older': 0.9999, 'PHSURV': 0.9999, 'HTRU_low_1757': 0.9982, '1534_survey': 0.9999, 'PMSURV': 0.5598, 'GBNCC': 0.9999, 'MMGPS-L': 0.9983 }
            deg_fac_1906 = {'PALFA_one_v_older': 0.9999, 'PHSURV': 0.9999, 'HTRU_low_1757': 0.9994, '1534_survey': 0.9999, 'PMSURV': 0.7337, 'GBNCC': 0.9999, 'MMGPS-L': 0.9995 }
            deg_fac_1208 = {'PALFA_one_v_older': 0.9999, 'PHSURV': 0.9999, 'HTRU_low_1757': 0.9998, '1534_survey': 0.9999, 'PMSURV': 0.8509, 'GBNCC': 0.9999, 'MMGPS-L': 0.9999 }
            
            if surv == "MMGPS-L":
                # Compute dimensions of beam.
                (major_axis,minor_axis)=random_beam_dimensions()

                #There are 2 possible configurations for this: the old tiling and the new tiling.
                #The new tiling covers the entire circle, but it has a smaller overlap.
                #The old tiling covers a central region of the beam with a tighter overlap.

              	tiling_decider=rand()

                if tiling_decider > 0.3915282003:

                    #Squares for 480 beams spread across an hexagon with inner diameter of 29.85.
                    #1.267911768=sqrt((6*(sin(30)/cos(30))*r**2)/480
                    (x,y)=(rand()*1.267911768-0.633955884,rand()*1.267911768-0.633955884)

                    beam_deg_factor=beam_deg_factor*Gaussian2D_degradation(x, y, 2*major_axis/60, 2*minor_axis/60)

                elif tiling_decider < 0.3915282003:

                    tiling_size=random_tiling_dimensions()

                    #Older tiling with 480 beams spread across a circle with an inner diameter fraction to be decided by the distribution.
                    if offset>tiling_size*1.05:

                        beam_deg_factor=0

                    else:

                        beam_region=np.sqrt((np.pi*(tiling_size*0.95)**2)/480)
                        (x,y)=(rand()*beam_region-beam_region/2,rand()*beam_region-beam_region/2)

                        beam_deg_factor=beam_deg_factor*Gaussian2D_degradation(x, y, 2*major_axis/60, 2*minor_axis/60)

            snr = snr*beam_deg_factor
                
            if snr==-1 or snr==-2 or snr==-3:
            	snr = snr
            elif psr_name=="0509":
                snr = snr * (deg_fac_0509[surv] ** 2)
            elif psr_name=="1102":
                snr = snr * (deg_fac_1102[surv] ** 2)
            elif psr_name=="1913":
                snr = snr * (deg_fac_1913[surv] ** 2)
            elif psr_name=="1946":
                snr = snr * (deg_fac_1946[surv] ** 2)
            elif psr_name=="1757":
                snr = snr * (deg_fac_1757[surv] ** 2)
            elif psr_name=="1534":
                snr = snr * (deg_fac_1534[surv] ** 2)
            elif psr_name=="0737A":
                snr = snr * (deg_fac_0737A[surv] ** 2)
            elif psr_name=="1756":
                snr = snr * (deg_fac_1756[surv] ** 2)
            elif psr_name=="1906":
                snr = snr * (deg_fac_1906[surv] ** 2)
            elif psr_name=="1208":
                snr = snr * (deg_fac_1208[surv] ** 2)

            ######################################################################################
            # add scintillation, if required
            # modifying S/N rather than flux is sensible because then
            # a pulsar can have same flux but change S/N in repeated surveys

            if scint:
                snr = s.scint(psr, snr)

            if snr > s.SNRlimit:
                ndet += 1
                psr.snr = snr
                survpop.population.append(psr)

                # check if the pulsar has been detected in other
                # surveys
                if not psr.detected:
                    # if not, set to detected and increment
                    # number of discoveries by the survey
                    psr.detected = True
                    s.discoveries += 1

            elif snr == -1.0:
                nsmear += 1
            elif snr == -2.0:
                nout += 1
            elif snr == -3.0:
                nbr += 1
            else:
                ntf += 1

        # report the results
        if not nostdout:
            print "Total pulsars in model = {0}".format(len(pop.population))
            print "Number detected by survey {0} = {1}".format(surv, ndet)
            print "Of which are discoveries = {0}".format(s.discoveries)
            print "Number too faint = {0}".format(ntf)
            print "Number smeared = {0}".format(nsmear)
            print "Number out = {0}".format(nout)
            if rratssearch:
                print "Number didn't burst = {0}".format(nbr)
            print "\n"


        d = Detections(ndet=ndet,
                       ntf=ntf,
                       nsmear=nsmear,
                       nout=nout,
                       nbr=nbr,
                       ndisc=s.discoveries)

        surveyPops.append([surv, survpop, d]) #MCiB: the only thing that find_alpha_general.py reads is.

    if allsurveyfile:
        allsurvpop = Population()
        allsurvpop.population = [psr for psr in pop.population if psr.detected]
        surveyPops.append([None, allsurvpop, None])

    return surveyPops


if __name__ == '__main__':
    """ 'Main' function; read in options, then survey the population"""
    # Parse command line arguments

    parser = argparse.ArgumentParser(
        description='Run a survey on your population model')
    parser.add_argument(
        '-f', metavar='fname', default='populate.model',
        help='file containing population model (def=populate.model')

    parser.add_argument(
        '-surveys', metavar='S', nargs='+', required=True,
        help='surveys to use to detect pulsars (required)')

    parser.add_argument(
        '--noresults', nargs='?', const=True, default=False,
        help='flag to switch off pickled .results file (def=False)')
    parser.add_argument(
        '--singlepulse', nargs='?', const=True, default=False,
        help='Rotating Radio Transients uses single pulse snr')

    parser.add_argument(
        '--asc', nargs='?', const=True, default=False,
        help='flag to create ascii population file (def=False)')

    parser.add_argument(
        '--summary', nargs='?', const=True, default=False,
        help='flag to create ascii summary file (def=False)')

    parser.add_argument(
        '--nostdout', nargs='?', const=True, default=False,
        help='flag to switch off std output (def=False)')

    parser.add_argument(
        '--allsurveys', nargs='?', const=True, default=False,
        help='write additional allsurv.results file (def=False)')

    parser.add_argument(
        '--scint', nargs='?', const=True, default=False,
        help='include model scintillation effects (def=False)')

    parser.add_argument(
        '--accel', nargs='?', const=True, default=False,
        help='use accel search for MSPs (def=False)')

    parser.add_argument(
        '--jerk', nargs='?', const=True, default=False,
        help='use accel & jerk search for MSPs (def=False)')

    args = parser.parse_args()

    # Load a model population
    population = loadModel(popfile=args.f)
    # run the population through the surveys
    surveyPopulations = run(population,
                            args.surveys,
                            nostdout=args.nostdout,
                            allsurveyfile=args.allsurveys,
                            scint=args.scint,
                            accelsearch=args.accel,
                            jerksearch=args.jerk,
                            rratssearch=args.singlepulse)

    # write the output files
    write(surveyPopulations,
          nores=args.noresults,
          asc=args.asc,
          summary=args.summary)
