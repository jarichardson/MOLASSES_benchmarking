#!/usr/bin/python -u

import numpy as np
import sys, os
import warnings

if len(sys.argv) != 2:
	print("Lava Flow Code path required!")
	print("    Syntax: ./benchmarking_lvl1_test2.py <molasses-model-path>")
	sys.exit();
	
#DEFINE LAVA FLOW MODEL FROM USER INPUT##########
model = sys.argv[1]
print "          Benchmarking Trail Level 1, Test 2"
print "  Test for self-similarity of flow on different slopes\n"
print "  Lava Flow Model: "+model+"\n"

#DEFINE OTHER FILES##############################
tstdir     = "bmark_lvl1_t2"
configfile = tstdir+"/"+model+".bm_l1t2.cfg"
demfile    = tstdir+"/bm_l1t2_DEM.asc"
outputfilea = tstdir+"/"+model+".bm_l1t2_out.xyz"
outputfiler = tstdir+"/"+model+".bm_l1t2_out.tif"
resultsfile= tstdir+"/"+model+".bm_l1t2_results.dat"
ventfile   = tstdir+"/"+model+".bm_l1t2_vent.dat"

if not os.path.isdir(tstdir):
	os.mkdir(tstdir)

#DEFINE CONSTANT DEM PARAMETERS##################
xurcorner     = 600.0
yurcorner     = xurcorner
xllcorner     = 0.0
yllcorner     = 0.0
channel_slope = 20
off_channel_azimuths = 15
pin_elev = 500.0      #location where both sides of channel are the same elev.
pin_x    = (xurcorner+xllcorner)*0.5     #DEM pin location
pin_y    = (yurcorner+yllcorner)*0.8

#DEFINE FLOW PARAMETERS##########################
residual_thickness = 10.0
vent_x             = pin_x
vent_y             = pin_y
total_volume       = 75000

#write out vent coordinates
vf = open(ventfile,"w")
vf.write("%0.2f\t%0.2f" % (vent_x,vent_y))
vf.close()


resf = open(resultsfile,"a")
resf.write("CellSize,PulseVol,FlowLength\n")

#WRITE DEM#######################################
for c in range(11):
	cellsize = (c+4)
	
	xurcorner     = 600.0
	if ((xurcorner-xllcorner)%cellsize)!=0:
		xurcorner = xurcorner-((xurcorner-xllcorner)%cellsize)
	
	yurcorner = xurcorner
	
	print "\nGrid Cell Resolution: "+str(cellsize)
		
	ncols         = int((yurcorner-yllcorner)/cellsize)
	nrows         = int((xurcorner-xllcorner)/cellsize)

	#Create ASCII DEM header
	dem_header  = "ncols   "+str(ncols)+"\n"
	dem_header += "nrows   "+str(nrows)+"\n"
	dem_header += "xllcorner       "+str(xllcorner)+"\n"
	dem_header += "yllcorner       "+str(yllcorner)+"\n"
	dem_header += "cellsize        "+str(cellsize)+"\n"

	demf = open(demfile,"w")
	demf.write(dem_header)

	#define DEM coordinates and associated geometry
	side_slopes   = channel_slope/np.cos(np.radians(off_channel_azimuths))  #degrees to south
	left_azimuth  = 180 - off_channel_azimuths #direction left channel side dips
	right_azimuth = 180 + off_channel_azimuths #direction right channel side dips

	xs = np.tile(np.arange(xllcorner,xurcorner,cellsize), nrows)
	ys = np.repeat(np.arange(yurcorner,yllcorner,-1*cellsize), ncols)
	dxs       = xs - pin_x #negative indicates left of centroid
	dys       = ys - pin_y #negative indicates below centroid
	distances = ((dxs)**2+(dys)**2)**0.5
	warnings.filterwarnings("ignore") #if distance=0, the next line will cause a warning
	dirs      = np.arccos(dys/distances) #deviation from north, blind to cw or ccw.
	warnings.filterwarnings("default") #restart warning detection
	dirs[np.where(distances==0)] = 0 #correct for location at centroid
	dirs[np.where(dxs<0)] -= 2*np.pi #correct for ccw angles
	dirs[np.where(dirs<0)] *= -1
	ondip_dzs  = np.tan(np.radians(side_slopes))*distances

	left_dzs  = ondip_dzs * np.cos(np.radians(left_azimuth)-dirs)
	right_dzs = ondip_dzs * np.cos(np.radians(right_azimuth)-dirs)

	elevs = pin_elev - right_dzs
	elevs[np.where(dxs<=0)] = pin_elev - left_dzs[np.where(dxs<=0)]

	ect = 0
	for e in elevs:
		demf.write("%0.2f\t" % e)
		ect+=1
		if ect%ncols==0:
			demf.write("\n") #new line for each row
	#close DEM file
	demf.close()
	sys.stdout.write("\r DEM written.                \n")

	#CREATE FLOW CONFIGURATION FILE##################
	for pulse_factor in range(10):
		pulse_volume  = (cellsize**2)*residual_thickness #"residual volume"
		#if pulse_factor==0:
			#pulse_volume *= 0.15
		#else:
		pulse_volume = 10**(pulse_factor/5.0 + 1.0)
	
	
	
		cfgf = open(configfile,"w")
		cfgf.write("DEM_FILE = "+demfile+"\n")
		cfgf.write("OUTFILE_A_THICKNESS = "+outputfilea+"\n")
		cfgf.write("OUTFILE_R_NEW_ELEV = "+outputfiler+"\n")
		cfgf.write("MODAL_THICKNESS = "+str(residual_thickness)+"\n")
		cfgf.write("NEW_VENT\n")
		cfgf.write("VENT_EASTING = "+str(vent_x)+"\n")
		cfgf.write("VENT_NORTHING = "+str(vent_y)+"\n")
		cfgf.write("VENT_PULSE_VOLUME = "+str(pulse_volume)+"\n")
		cfgf.write("VENT_TOTAL_VOLUME = "+str(total_volume))
		cfgf.close()

		#RUN LAVA FLOW CODE############################
		sys.stdout.write(" Running Lava Flow Code (Pulse %0.3e)...\n" % pulse_volume)
		err_out = os.system("./"+model+" "+configfile+"| grep -i 'ERROR'")
		if err_out!=256:
			print "Critical error in Benchmark. Ending prematurely\n"
			break
	
		#ANALYZE FLOW OUTPUT###########################
		#  Read xyz output
		flowxs,flowys,flowzs = np.loadtxt(outputfilea,unpack=True)
	
		#  Locate farthest point from vent (runout)
		flow_dists = ((flowxs-vent_x)**2+(flowys-vent_y)**2)**0.5
		maxflowdist = np.max(flow_dists)
		resf.write("%0.3f\t%0.3f\t%0.3f\n" % (cellsize,np.log10(pulse_volume),np.max(flow_dists)))

	#MAP FLOW OUTPUT###############################
	sys.stdout.write("Mapping Flow...\n")

	cptfile = tstdir+"/"+model+(".%02d-m.cpt" % cellsize)
	gradfile = tstdir+"/"+model+(".%02d-m.gradient" % cellsize)
	mapfile = tstdir+"/"+model+(".%02d-m.eps" % cellsize)

	os.system("makecpt -Cgray -I -T"+str(np.min(elevs))+"/"+str(np.max(elevs))+"/2 >"+cptfile)
	os.system("grdgradient "+outputfiler+" -A0/270 -G"+gradfile)
	os.system("grdimage "+outputfiler+" -C"+cptfile+" -I"+gradfile+" `grdinfo -I1 "+demfile+"` -JX6i -Ba100WSen -K > "+mapfile)
	os.system("grdcontour "+outputfiler+" -A20 -R -J -O -K >> "+mapfile)
	#os.system("psxy " + outputfilea + " -R -J -Sc0.05c -O -K >> " + mapfile)
	os.system("psxy " + ventfile + " -R -J -O -St0.2c -G66 -Wwhite >> " + mapfile)
	os.system("ps2raster "+mapfile+" -A -Tg -P ")
	os.remove(mapfile)
	os.remove(cptfile)
	os.remove(gradfile)

resf.close()
os.remove(outputfiler)
os.remove(configfile)
os.remove(outputfilea)
os.remove(ventfile)
