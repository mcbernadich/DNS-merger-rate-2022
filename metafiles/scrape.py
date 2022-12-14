import numpy as np
#from ast import literal_eval
import json
import glob

folders=glob.glob("/beegfs/DATA/TRAPUM/SCI-20200703-MK-01/**")

collected_data=open("collected_data.txt","w")
debugging_info=open("debugging_info.txt","w")

for folder in folders:

	metafiles=glob.glob(folder+"/**/*.meta")

	for metafile in metafiles:

		print(metafile)

		with open(metafile, 'r') as f:
			try:
				data = json.load(f)  
			except:
				continue

		beam_data = data['beamshape']
		bore_data = data['boresight']
		time = data['utc_start'].split(" ")
		date=time[0].split("/")
		hour=time[1].split(":")
		time=float(date[0])+float(date[1])/12+float(date[2])/365+float(hour[0])/8760+float(hour[1])/525600+float(hour[2])/31536000

		ra=bore_data.split(",")[2].split(":")
		dec=bore_data.split(",")[3].split(":")

		ra=180*(float(ra[0])+float(ra[1])/60+float(ra[2])/3600)/12
		if float(dec[0])<0:
			dec=float(dec[0])-float(dec[1])/60-float(dec[2])/3600
		elif float(dec[0])>0:
			dec=float(dec[0])+float(dec[1])/60+float(dec[2])/3600


		beam_vals = list(data['beams'].values())
		beam_keys = list(data['beams'].keys())		
		largest=0
		i=0
		for val in beam_vals:

			if val.split(",")[0]=="unset":
				continue

			beam_ra=val.split(",")[2].split(":")
			beam_dec=val.split(",")[3].split(":")

			beam_ra=180*(float(beam_ra[0])+float(beam_ra[1])/60+float(beam_ra[2])/3600)/12
			if float(beam_dec[0])<0:
				beam_dec=float(beam_dec[0])-float(beam_dec[1])/60-float(beam_dec[2])/3600
			elif float(beam_dec[0])>0:
				beam_dec=float(beam_dec[0])+float(beam_dec[1])/60+float(beam_dec[2])/3600

			beam_dist=60*np.sqrt(((beam_ra-ra)*np.cos(np.pi*dec/180))**2+(beam_dec-dec)**2)

			if beam_dist>largest:
				largest=beam_dist
				furthest_ra=beam_ra
				furthest_dec=beam_dec
				largest_key=beam_keys[i]

			i=i+1

		collected_data.write(beam_vals[0].split(",")[0]+","+str(ra)+","+str(dec)+","+str(3600*float(beam_data.split(",")[1].split(":")[1]))+","+str(3600*float(beam_data.split(",")[0].split(":")[1]))+","+str(time)+","+str(largest)+"\n")
		debugging_info.write(largest_key+","+str(furthest_ra)+","+str(furthest_dec)+","+str(np.cos(np.pi*dec/180))+"\n")
