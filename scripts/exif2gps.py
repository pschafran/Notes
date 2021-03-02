#! /usr/bin/python
# (C) Peter Schafran 8 March 2017


#This script extracts GPS coordinates and elevation from EXIF data in iPhone photos and writes CSV and GPX files containing these data. 
#Warnings: Input files must have unique file names. Output filenames should be changed to avoid being overwritten. Only tested for iPhone 5S photos
#Dependency: Requires ExifRead (https://pypi.python.org/pypi/ExifRead) and requires EXIF.py is in $PATH
## 30 Nov 2017 -- modified to prevent crashing when some photos don't have GPS or altitude info. Warning message is generated  in STDOUT giving the trouble file when this happens.

import sys, os, re, subprocess

filedict={}
filelist=sys.argv[1:]
exiflist=[]

exifsave = input("Save EXIF files? 1=Yes 2=No ")


#Create CSV and KML output files

CSV_outfile = open("iPhone_GPS.csv", "w")
GPX_outfile = open("iPhone_GPS.gpx", "w")

CSV_outfile.write("File_Name,Latitude,Longitude,Altitude(m)\n")
GPX_outfile.write('''<?xml version="1.0" encoding="UTF-8" standalone="no" ?><gpx>''')

#create dictionary keys based on unique ID in file name
for file in filelist:
	#print file
	splitfile = re.split('[.]', file)
	filedict[splitfile[0]]=[]
	file_out = open('%s.EXIF' %(splitfile[0]), 'w')
	subprocess.call(['EXIF.py', '%s' %(file)], stdout=file_out)
	file_out.close()
	exiflist.append('%s.EXIF' %(splitfile[0]))
	
exif_dict ={}
for exiffile in exiflist:
	open_file = open(exiffile, 'r')
	exif_dict[exiffile]={}
	for line in open_file:
		splitline = line.strip("\n").split(":")
		if len(splitline) >= 2:
			exif_dict[exiffile][splitline[0]]=splitline[1:]
	try:
		latdms = exif_dict[exiffile]["GPS GPSLatitude (Ratio)"]
		longdms = exif_dict[exiffile]["GPS GPSLongitude (Ratio)"]
		latref = exif_dict[exiffile]["GPS GPSLatitudeRef (ASCII)"]
		longref = exif_dict[exiffile]["GPS GPSLongitudeRef (ASCII)"]
		altitude = exif_dict[exiffile]["GPS GPSAltitude (Ratio)"]
	
		latdms_split=str(latdms).split(",")
		longdms_split=str(longdms).split(",")
	except:
		print "WARNING: Could not process %s" %(exiffile)
		latdms_split = "[' [0", ' 0', " 0/1]']"
		longdms_split = "[' [0", ' 0', " 0/1]']"
		latref = ["null"]
		longref = ["null"]
		altitude = [0]
	
	
####Latitude Manipulations and Calculations
	latdms_d = latdms_split[0]
	latdms_m = latdms_split[1]
	latdms_s = latdms_split[2]
	
	latdms_d2 = latdms_d[4:]
	latdms_s2 = latdms_s[:-3]
	
	latdms_s_split = latdms_s2.split("/")
	if len(latdms_s_split) == 2:
		latdms_s_decimal = float(latdms_s_split[0])/float(latdms_s_split[1])
	elif len(latdms_s_split) == 1:
		latdms_s_decimal = float(latdms_s_split[0])
	
	#Convert degrees, minutes, seconds into decimal degrees
	lat_decdeg = (((latdms_s_decimal/60)+float(latdms_m))/60)+float(latdms_d2)
	
	#If latitude reference is North, do not change; if South, make negative
	latref2 = str(latref[0])
	if "N" in latref2:
		lat_decdeg2 = lat_decdeg*1
	elif "S" in latref2:
		lat_decdeg2 = lat_decdeg*-1
	elif "null" in latref2:
		lat_decdeg2 = "null"



####Longitude Manipulations and Calculations
	longdms_d = longdms_split[0]
	longdms_m = longdms_split[1]
	longdms_s = longdms_split[2]
	
	longdms_d2 = longdms_d[4:]
	longdms_s2 = longdms_s[:-3]
	
	longdms_s_split = longdms_s2.split("/")
	if len(longdms_s_split) == 2:
		longdms_s_decimal = float(longdms_s_split[0])/float(longdms_s_split[1])
	elif len(latdms_s_split) == 1:
		longdms_s_decimal = float(longdms_s_split[0])
	
	#Convert degrees, minutes, seconds into decimal degrees
	long_decdeg = (((longdms_s_decimal/60)+float(longdms_m))/60)+float(longdms_d2)
	
	#If latitude reference is North, do not change; if South, make negative
	longref2 = str(longref[0])
	if "E" in longref2:
		long_decdeg2 = long_decdeg*1
	elif "W" in longref2:
		long_decdeg2 = long_decdeg*-1
	elif "null" in longref2:
		long_decdeg2 = "null"
		
	
####Altitude Manipulation
	altitude2 = str(altitude[0])
	
	altitude_split = altitude2.split("/")
	if len(altitude_split) == 2:
		altitude_decimal = float(altitude_split[0])/float(altitude_split[1])
	elif len(altitude_split) == 1:
		altitude_decimal = float(altitude_split[0])
	
	
	
	#print ("%s, %s, %s meters") %(lat_decdeg2, long_decdeg2, altitude2)
	exiffilename = exiffile.split(".")
	CSV_outfile.write("%s,%s,%s,%s\n" %(str(exiffilename[0]),lat_decdeg2,long_decdeg2,altitude_decimal))
	
	GPX_outfile.write('''<wpt lat="%s" lon="%s"><ele>%s</ele><name>%s</name><sym>Flag, Blue</sym></wpt>''' %(lat_decdeg2,long_decdeg2,altitude_decimal,str(exiffilename[0])))
	
	
	
	if exifsave == 2:
		subprocess.call(['rm','%s' %(exiffile)])
	
CSV_outfile.close()
GPX_outfile.write("</gpx>")
GPX_outfile.close()







