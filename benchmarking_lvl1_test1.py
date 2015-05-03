#!/usr/bin/python -u

#Benchmarking Test for MOLASSES lava flow codes
#LEVEL 1:
#        Self-similarity tests
#TEST  1:
#        Slope direction test
#
#This trial tests the ability of a lava flow code to replicate the same flow
#morphology independent of slope direction. A lava flow is runout with the same
#parameters over a simple sloped (18 degrees, Osmond and Griffiths 2001, Fig 8)
#DEM that rotates 180 times at 2 degree increments.

#Created by Jacob Richardson, April 2015 (jarichardson@mail.usf.edu)

#ALGORITHM
#Define lava flow model from user input
#Define azimuths and inclination of slope, path to configuration file, path to
#       temporary DEM, path to output xyz.
#Create lava flow configuration file for MOLASSES.
#
#For each azimuth:
#  Create a DEM of azimuth, slope
#  Run Lava Flow Code
#  Read xyz output
#  Locate farthest point from vent (runout)
#  Define line from runout point to vent (runout path)
#  Locate farthest points on either side from runout path, find distance (width)
#  Report runout distance, runout width to output file
#
#Remove DEM, xyz output

import numpy as np
import sys, os
import warnings


if len(sys.argv) != 2:
	print("Lava Flow Code path required!")
	print("    Syntax: ./benchmarking_lvl1_test1.py <molasses-model-path>")
	sys.exit();


#DEFINE LAVA FLOW MODEL FROM USER INPUT##########
model = sys.argv[1]
print "          Benchmarking Trail Level 1, Test 1"
print "  Test for self-similarity of flow on different slopes\n"
print "  Lava Flow Model: "+model+"\n"

#DEFINE OTHER FILES##############################
tstdir     = "bmark_lvl1_t1"
configfile = tstdir+"/"+model+".bm_l1t1.cfg"
demfile    = tstdir+"/"+model+".bm_l1t1_DEM.asc"
outputfile = tstdir+"/"+model+".bm_l1t1_out.xyz"
resultsfile= tstdir+"/"+model+".bm_l1t1_results.dat"
ventfile   = tstdir+"/"+model+".bm_l1t1_vent.dat"

#reset test directory
#os.system("rm -rf "+tstdir)
#os.system("mkdir "+tstdir)


#DEFINE CONSTANT DEM PARAMETERS##################
ncols         = 750
nrows         = 750
xllcorner     = 0.0
yllcorner     = 0.0
cellsize      = 1.0
centroid_elev = 500.0
centroid_x    = (ncols*cellsize/2)+xllcorner     #DEM centroid location
centroid_y    = (nrows*cellsize/2)+yllcorner

#DEFINE MODEL PARAMETERS#########################
dem_slope     = 18.0  #degrees
degree_sweep  = 90.0  #total arc width to test
degree_step   =  5.0  #rotation resolution

#Create ASCII DEM header
dem_header  = "ncols   "+str(ncols)+"\n"
dem_header += "nrows   "+str(nrows)+"\n"
dem_header += "xllcorner       "+str(xllcorner)+"\n"
dem_header += "yllcorner       "+str(yllcorner)+"\n"
dem_header += "cellsize        "+str(cellsize)+"\n"

#DEFINE FLOW PARAMETERS##########################
residual_thickness = 1.0
vent_x             = centroid_x
vent_y             = centroid_y
pulse_volume       = (cellsize**2)*residual_thickness #"residual volume"
total_volume       = 3000

#write out vent coordinates
vf = open(ventfile,"w")
vf.write("%0.2f\t%0.2f" % (vent_x,vent_y))
vf.close()

#CREATE FLOW CONFIGURATION FILE##################
cfgf = open(configfile,"w")
cfgf.write("DEM_FILE = "+demfile+"\n")
cfgf.write("OUTFILE_A_THICKNESS = "+outputfile+"\n")
cfgf.write("MODAL_THICKNESS = "+str(residual_thickness)+"\n")
cfgf.write("NEW_VENT\n")
cfgf.write("VENT_EASTING = "+str(vent_x)+"\n")
cfgf.write("VENT_NORTHING = "+str(vent_y)+"\n")
cfgf.write("VENT_PULSE_VOLUME = "+str(pulse_volume)+"\n")
cfgf.write("VENT_TOTAL_VOLUME = "+str(total_volume))
cfgf.close()


#RUN FLOW CODE OVER ROTATING DEM#################

#define slope directions
azimuth_step = np.radians(degree_step) #radians
azimuths = np.arange(0,np.radians(degree_sweep),azimuth_step)

#define empty results array
flowdimensions = np.zeros((len(azimuths),4))

#define DEM coordinates and associated geometry
xs = np.tile(np.arange(xllcorner,(ncols*cellsize)+xllcorner,cellsize), nrows)
ys = np.repeat(np.arange((nrows*cellsize)+yllcorner,yllcorner,-1*cellsize), ncols)
dxs       = xs - centroid_x #negative indicates left of centroid
dys       = ys - centroid_y #negative indicates below centroid
distances = ((dxs)**2+(dys)**2)**0.5
warnings.filterwarnings("ignore") #if distance=0, the next line will cause a warning
dirs      = np.arccos(dys/distances) #deviation from north, blind to cw or ccw.
warnings.filterwarnings("default") #restart warning detection
dirs[np.where(distances==0)] = 0 #correct for location at centroid
dirs[np.where(dxs<0)] -= 2*np.pi #correct for ccw angles
dirs[np.where(dirs<0)] *= -1
ondip_dzs  = np.tan(np.radians(dem_slope))*distances 

act = 0
for azimuth in azimuths: ########################
	print("Dip Direction: %d degrees (from N)"  % np.degrees(azimuth))
	
	#Write DEM File################################
	demf = open(demfile,"w")
	#Write DEM header
	demf.write(dem_header)
	
	#calculate elevation given azimuth and coordinates
	dzs = ondip_dzs * np.cos(azimuth-dirs)
	elevs = centroid_elev - dzs
	
	#print elevations to DEM
	sys.stdout.write("\r Writing DEM...")
	ect = 0
	for e in elevs:
		demf.write("%0.2f\t" % e)
		ect+=1
		if ect%ncols==0:
			demf.write("\n") #new line for each row
	
	#close DEM file
	demf.close()
	sys.stdout.write("\r DEM written.                \n")
	
	#RUN LAVA FLOW CODE############################
	sys.stdout.write(" Running Lava Flow Code...\n ")
	os.system("./"+model+" "+configfile+"| grep 'ERROR'")
	
	
	#MAP FLOW OUTPUT###############################
	sys.stdout.write("Mapping Flow...\n")
	mapfile = tstdir+"/"+model+(".%03d.eps" % np.degrees(azimuth))
	os.system("grdcontour "+demfile+" -C50 -A100 `grdinfo -I1 "+demfile+"` -JX3i -Ba200WSen -K > "+mapfile)
	os.system("psxy " + outputfile + " -R -J -Sc0.05c -Gblack -O -K >> " + mapfile)
	os.system("psxy " + ventfile + " -R -J -O -St0.2c -G66 -Wwhite >> " + mapfile)
	os.system("ps2raster "+mapfile+" -A -Tg -P ")
	os.system("rm "+mapfile)
	
	#ANALYZE FLOW OUTPUT###########################
	#  Read xyz output
	flowxs,flowys,flowzs = np.loadtxt(outputfile,unpack=True)
	
	#  Locate farthest point from vent (runout)
	flow_dists = ((flowxs-vent_x)**2+(flowys-vent_y)**2)**0.5
	flowdimensions[act][0] = np.max(flow_dists)
	flowy_far = np.average(flowys[np.where(flow_dists==flowdimensions[act][0])])
	flowx_far = np.average(flowxs[np.where(flow_dists==flowdimensions[act][0])])
	
	#  Locate farthest points on either side from runout path, find distance (width)
	flow_dx = vent_x - flowx_far
	flow_dy = vent_y - flowy_far
	#Calculate distance off axis for each point
	flow_widths = flow_dy*flowxs - flow_dx*flowys + vent_x*flowy_far - vent_y*flowx_far
	flow_widths /= flowdimensions[act][0]
	flowdimensions[act][1] = np.max(flow_widths) - np.min(flow_widths)
	flowdimensions[act][2] = flow_dx
	flowdimensions[act][3] = flow_dy
	
	#increment azimuth count
	act+=1

#CALCULATE SELF-SIMILARITY#######################
#flow length variance
flow_ave_L = np.average(flowdimensions[:,0]) #mean
flow_std_L = np.std(flowdimensions[:,0])     #standard deviation
flow_CV_L = 100*flow_std_L/flow_ave_L        #coefficient of variance

#flow width variance
flow_ave_W = np.average(flowdimensions[:,1])
flow_std_W = np.std(flowdimensions[:,1])
flow_CV_W = 100*flow_std_W/flow_ave_W

#flow aspect ratio variance
aspectratios = flowdimensions[:,0]/flowdimensions[:,1]
flow_ave_AR = np.average(aspectratios)
flow_std_AR = np.std(aspectratios)
flow_CV_AR = 100*flow_std_AR/flow_ave_AR


#flow-DEM azimuth variance
flowdirection = np.arctan(flowdimensions[:,2]/flowdimensions[:,3])
flowvariance = np.degrees(flowdirection - azimuths)
flowvariance[np.where(flowvariance>=180.0)] -= 180.0
flow_ave_diroff = np.average(abs(flowvariance))

#WRITE OUTPUT####################################
resf = open(resultsfile,"w")
resf.write("#Benchmark LVL1 TEST1 Results for Lava Flow Model: "+model+"\n")
resf.write("#Azimuth(deg) FlowLength(m) FlowWidth(m) AspectRatio DirectionOffset\n")
resf.write("#Average \t%0.3f\t%0.3f\t%0.3f\t%0.3f\n" % (flow_ave_L,flow_ave_W,flow_ave_AR,flow_ave_diroff))
resf.write("#Std Dev \t%0.3f\t%0.3f\t%0.3f\n" % (flow_std_L,flow_std_W,flow_std_AR))
resf.write("#Variance\t%0.3f%%\t%0.3f%%\t%0.3f%%\n" % (flow_CV_L,flow_CV_W,flow_CV_AR))
for i in range(len(azimuths)):
	#write azimuth, flow length, flow width
	resf.write("%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n" % (np.degrees(azimuths[i]),flowdimensions[i][0],flowdimensions[i][1],aspectratios[i],flowvariance[i]))
	
resf.close()

print("\n  Flow Dimensions with respect to Slope Direction output to:")
print "    "+resultsfile

#Example code to plot the DEM and flow on top.
#os.system("grdimage bmark_lvl1_t1/bm_l1t1_DEM.asc -R0/750/0/750 -Ba100 -JX8i -Ccpt_bm.cpt -K -V -P > abm.eps")
#os.system("psxy -R -J widths.xyz -Sc0.1c -Cdist.cpt -O -V >> abm.eps")

#gmtselect bmark_lvl1_t1/bm_l1t1_results.dat -o0,2 > bmark_lvl1_t1/width.dat
#psxy bmark_lvl1_t1/bm_l1t1_results.dat -R0/90/0/150 -Ba22.5/10WS -JX6i -P -K > abm.eps
#psxy bmark_lvl1_t1/width.dat -R -J -O >> abm.eps 


#Create an animated gif and remove individual maps
print "  Creating rotating gif at:\n    "+tstdir+"/"+model
os.system("convert -dispose background -delay "+("%d" % (8.0*degree_step))+" -loop 0 "+tstdir+"/"+model+".*.png -resize 40% " + tstdir+"/"+model+".gif")


#clean benchmark directory
os.system("rm "+configfile+" "+demfile+" "+outputfile+" "+ventfile+" "+tstdir+"/"+model+".*.png")

print "\n***Done."
