import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import poisson
from numpy.random import rand

def displaced_poisson(k,lamb,a,A):
    '''poisson function, parameter lamb is the fit parameter'''
    return A*lamb**(k-a)*np.exp(-lamb)/factorial(k-a)

def gaussian(x,sigma,x0,A):

	gaussian = A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

	return gaussian

def bimodal_gaussian(x,sigma1,x1,A1,sigma2,x2,A2):

	gaussian1 = A1 * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2))
	gaussian2 = A2 * np.exp(-(x - x2) ** 2 / (2 * sigma2 ** 2))

	return gaussian1 +  gaussian2

def skewed_gaussian(x,sigma,x0,A,sk):

	gaussian = A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
	skew = (1+erf(sk*(x-x0)/(sigma*np.sqrt(2))))

	return gaussian*skew

def bimodal_skewed_gaussian(x,sigma1,x1,A1,sk1,sigma2,x2,A2,sk2):

	gaussian1 = A1 * np.exp(-(x - x1) ** 2 / (2 * sigma1 ** 2))
	gaussian2 = A2 * np.exp(-(x - x2) ** 2 / (2 * sigma2 ** 2))
	skew1 = (1+erf(sk1*(x-x1)/(sigma1*np.sqrt(2))))
	skew2 = (1+erf(sk2*(x-x2)/(sigma2*np.sqrt(2))))

	return gaussian1*skew1+gaussian2*skew2


collected_data=np.loadtxt("collected_data.txt",delimiter=",",dtype=str).T

major=np.maximum(collected_data[3].astype(float),collected_data[4].astype(float))
minor=np.minimum(collected_data[3].astype(float),collected_data[4].astype(float))

indices=[]

for element in collected_data[0]:
	if "MSGPS_L_" in element:
		indices.append(True)
	else:
		indices.append(False)

major=major[indices]
minor=minor[indices]
ra=collected_data[1,indices].astype(float)
dec=collected_data[2,indices].astype(float)
time=collected_data[5,indices].astype(float)
tiling_size=collected_data[6,indices].astype(float)

#major=major[(minor>20) & (minor<30)]
#ra=collected_data[1,(minor>20) & (minor<30)]
#dec=collected_data[2,(minor>20) & (minor<30)]
#time=collected_data[5,(minor>20) & (minor<30)]
#tiling_size=collected_data[6,(minor>20) & (minor<30)]
#minor=minor[(minor>20) & (minor<30)]

#major=major[(ra<260) & (ra>120)]
#minor=minor[(ra<260) & (ra>120)]
#dec=dec[(ra<260) & (ra>120)]
#time=time[(ra<260) & (ra>120)]
#tiling_size=tiling_size[(ra<260) & (ra>120)]
#ra=ra[(ra<260) & (ra>120)]

#major=major[dec<-30]
#minor=minor[dec<-30]
#ra=ra[dec<-30]
#time=time[dec<-30]
#tiling_size=tiling_size[dec<-30]
#dec=dec[dec<-30]

size_axis=np.arange(0,100,0.2)
major_histo=np.zeros(np.size(size_axis),int)

for size in major:
	i=0
	for bn in size_axis:
		if size>bn and size<size_axis[i+1]:
			major_histo[i]=major_histo[i]+1
			break
		i=i+1

minor_histo=np.zeros(np.size(size_axis),int)

for size in minor:
	i=0
	for bn in size_axis:
		if size>bn and size<size_axis[i+1]:
			minor_histo[i]=minor_histo[i]+1
			break
		i=i+1

(major_parameters, major_cov) = curve_fit(skewed_gaussian, size_axis+0.1, major_histo, p0=[9,30,500,4])
(minor_parameters, minor_cov) = curve_fit(bimodal_skewed_gaussian, size_axis+0.1, minor_histo, p0=[3.4,23,2000,4,3.4,27,2000,4])

print("Major axis sigma:",major_parameters[0],", center:",major_parameters[1],", amplitude:",major_parameters[2]," skewness:",major_parameters[3])
print("Minor axis 1 sigma:",minor_parameters[0],", center:",minor_parameters[1],", amplitude:",minor_parameters[2],", skewness:",minor_parameters[3])
print("Minor axis 2 sigma:",minor_parameters[4],", center:",minor_parameters[5],", amplitude:",minor_parameters[6],", skewness:",minor_parameters[7])
print("Minor axis 1/2 relative amplitudes:",minor_parameters[2]/minor_parameters[6])

print("Pointings before tiling change (October 5, 2022):",np.size(time[time<2021.847]))
print("pointings after tiling change (October 5, 2022):",np.size(time[time>2021.847]))

plt.scatter(ra,major,marker=".",label="Semi-major axis")
plt.scatter(ra,minor,marker=".",label="Semi-minor axis")
plt.xlabel("RA (degrees)")
plt.ylabel("Beam size (arcsec)")
plt.legend()
plt.show()

plt.scatter(dec,major,marker=".",label="Semi-major axis")
plt.scatter(dec,minor,marker=".",label="Semi-minor axis")
plt.xlabel("DEC (degrees)")
plt.ylabel("Beam size (arcsec)")
plt.legend()
plt.show()

plt.scatter(time,major,marker=".",label="Major axis")
plt.scatter(time,minor,marker=".",label="Minor axis")
plt.xlabel("time (years)")
plt.ylabel("Beam size (arcsec)")
plt.legend()
plt.show()

plt.scatter(time,tiling_size,marker=".")
plt.xlabel("time (years)")
plt.ylabel("Tiling size (arcsec)")
plt.show()

plt.scatter(ra,tiling_size,marker=".")
plt.xlabel("RA (degrees)")
plt.ylabel("Tiling size (arcsec)")
plt.show()

plt.scatter(dec,tiling_size,marker=".")
plt.xlabel("DEC (degrees)")
plt.ylabel("Tiling size (arcmin)")
plt.show()

plt.scatter(major,minor,marker=".")
plt.xlabel("Major axis (arcsec)")
plt.ylabel("Minor axis (arcsec)")
plt.show()

model=np.arange(0,100,0.01)

plt.plot(size_axis+0.1,major_histo,label="Major axis")
plt.plot(model,skewed_gaussian(model,major_parameters[0],major_parameters[1],major_parameters[2],major_parameters[3]),"k-",label="Major axis")
plt.plot(size_axis+0.1,minor_histo,label="Minor axis")
plt.plot(model,bimodal_skewed_gaussian(model,minor_parameters[0],minor_parameters[1],minor_parameters[2],minor_parameters[3],minor_parameters[4],minor_parameters[5],minor_parameters[6],minor_parameters[7]),"k-",label="Minor axis")
plt.xlabel("Beam size (arcsec)")
plt.ylabel("Pointings")
plt.xlim(20,65)
plt.show()

major_x_rand=rand(5000)*45+20
major_y_rand=rand(5000)*65
major_y=skewed_gaussian(major_x_rand,major_parameters[0],major_parameters[1],major_parameters[2],major_parameters[3])
major_x_take=major_x_rand[major_y_rand<major_y]

minor_x_rand=rand(5000)*15+20
minor_y_rand=rand(5000)*252
minor_y=bimodal_skewed_gaussian(minor_x_rand,minor_parameters[0],minor_parameters[1],minor_parameters[2],minor_parameters[3],minor_parameters[4],minor_parameters[5],minor_parameters[6],minor_parameters[7])
minor_x_take=minor_x_rand[minor_y_rand<minor_y]

size=min(np.size(major_x_take),np.size(minor_x_take))

new_major=major_x_take[:size]
new_minor=minor_x_take[:size]

plt.scatter(major,minor,marker=".",label="MMGPS-L")
plt.scatter(new_major,new_minor,marker=".",label="simulated")
plt.xlabel("Semi-major axis (arcsec)")
plt.ylabel("Semi-minor axis (arcsec)")
plt.legend()
plt.tight_layout()
plt.show()

size_axis=np.arange(5,25,0.2)
tiling_histo=np.zeros(np.size(size_axis),int)

for size in tiling_size[time<2021.846]:
	i=0
	for bn in size_axis:
		if size>bn and size<size_axis[i+1]:
			tiling_histo[i]=tiling_histo[i]+1
			break
		i=i+1

(tiling_parameters, major_cov) = curve_fit(bimodal_gaussian, size_axis+0.1, tiling_histo, p0=[1.193682939192035,12.582730128392928,500,3,12.582730128392928,500])
print("Tiling size 1 sigma:",tiling_parameters[0],", center1:",tiling_parameters[1],", amplitude:",tiling_parameters[2])
print("Tiling size 1 sigma:",tiling_parameters[3],", center2:",tiling_parameters[4],", amplitude:",tiling_parameters[5])

model=np.arange(2,25,0.01)

plt.plot(size_axis+0.1,tiling_histo)
plt.plot(model,bimodal_gaussian(model,tiling_parameters[0],tiling_parameters[1],tiling_parameters[2],tiling_parameters[3],tiling_parameters[4],tiling_parameters[5]),"k-")
plt.xlabel("Tiling size (arcmin)")
plt.ylabel("Pointings")
plt.show()
