#!/usr/bin/env python
 
import time, os, requests, re
#################DOWNLOAD TIDE PREDICTION DATA#################
 
current_date= time.strftime("%Y%m%d %H:%M")	#store current date and time as a variable
print "Tide prediction data downloading from NOAA..."
#request the data we want from the following url using specified parameters: 'http://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=20141108&end_date=20141109&datum=MLLW&station=8638610&time_zone=GMT&units=english&interval=&format=csv'
payload = {'product': 'predictions', 'application': 'NOS.COOPS.TAC.WL', 'begin_date': current_date, 'range': '48', 'datum': 'NAVD', 'station': '8638610', 'time_zone': 'GMT', 'units': 'english', 'format': 'csv'}
r = requests.get("http://tidesandcurrents.noaa.gov/api/datagetter", params=payload)
#print (r.url)
#print (r.text)
print "Writing data to temp file..."
outfile = open('predict.csv', 'w')	#open file to write to. Need this temp file to make following for loop function correctly.
outfile.write(r.text)	#write the url curled values to a .csv called predict
outfile.close()	#close the file
infile = open('predict.csv', 'r')	#use predict.csv to read into the for loop line by line
outfile2 = open('SewellsPt_Predictions_48hr.csv', 'w')
#strippredict = predict.strip('\n')
#print strippredict
#####################DOWNLOAD SURGE DATA#####################
print "Surge data downloading from NWS..."
#request the data we want from the following url using specified parameters: 'http://www.nws.noaa.gov/mdl/etsurge/index.php?page=stn&region=me&datum=msl&list=&map=0-48&type=both&stn=vahamp'
payload = {'page': 'stn', 'region': 'me', 'datum': 'msl', 'list': '', 'map': '0-48', 'type': 'both', 'stn': 'vahamp'}
r2 = requests.get("http://www.nws.noaa.gov/mdl/etsurge/index.php", params=payload)
#print (r.url)
#print (r.text)
print "Writing data to temp file..."
outfile3 = open('~surge.csv', 'w')	#open file to write to. Need this temp file to make following for loop function correctly.
outfile3.write(r2.text)	#write the url curled values to a .csv called surge	
outfile3.close()	#close the file
infile2 = open('~surge.csv', 'r')	#use predict.csv to read into the for loop line by line
outfile4 = open('surges.csv', 'w')
for line in infile2:
	if re.match('^\d+/\d+', line):
		#print line
		outfile4.write(line)
	else:
		continue
outfile4.close()

#################MATCH UP AND COMBINE TIDE PREDICTIONS WITH SURGE DATA#################
linenum = 0	#start the counter
SewellsPtLat = "36.9467 N"	#set the lat
SewellsPtLong = "76.3300 W" #set the lon
outfile2.write('SewellsPtLat, SewellsPtLong, DateTime(GMT), NAVD(ft), ShortDate(Zulu), Surge(ft), FloodLevel(NAVD+Surge)\n')	#write header line
print "Reading in surge data..."
infile3 = open('surges.csv', 'r')
surgedict = {}
for surgeline in infile3:
	#print surgeline
	surgesplit = surgeline.strip(',\n').split(', ') #strip end of line characters and extra comma. split at comma+space
	#print surgesplit
	surgedate = surgesplit[0]
	surge = float(surgesplit[1]) #convert surge into a float
	surgedict[surgedate] = surge
	#print surgedate
	#print surge
	#print surgedict
	
print "Matching surge and tide prediction data..."	
for line in infile:	#loop through each line in the prediction file
	if linenum > 0:	#skip the header
		#print line
		month = line[5:7]
		day = line[8:10]
		hour = line[11:13]		
		#print month, day, hour
		datetime = '%s/%s %sZ' %(month, day, hour)
		splitline = line.strip('\n').split(',')	#strip and split the line
		#print splitline[0]
		#print splitline[1]
		NAVD = float(splitline[1])
		
		for key in surgedict.keys():
			if datetime == key:
				floodlevel = NAVD + surgedict[key]
				#print datetime, surgedict[key], floodlevel
		
				outfile2.write('%s,%s,%s,%s,%s,%.2f, %.2f\n' %(SewellsPtLat, SewellsPtLong, splitline[0], splitline[1], datetime, surgedict[key], floodlevel)) 	#write the lat, lon, date time, and NAVD
			else:
				continue
	else:
		linenum += 1
		

outfile2.close()
print "Finished."



