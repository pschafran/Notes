#! /usr/bin/python

#Description: Updates all pages on website with same section (e.g. header, sidebar, footer). Run in directory with site's html files

#Usage: websiteUpdater.py *.html

import sys, os

savePath = "./UpdatedFiles/"
switch = 0
websitefiles = sys.argv[1:]
for file in websitefiles:
	openfile = open(file, 'r')
	fileName = file.strip('\n').split('.')
	newfile = open(os.path.join(savePath, "%s.html" %(fileName[0])), 'w')
	newSideBar = open("sidebar.txt", 'r')
	for line in openfile:
		if switch == 0:
			if "<!-- Sidebar -->" in line:
				newfile.write(line)
				switch = 1
			else:
				newfile.write(line)
				#print line
		if switch == 1:
			if "</footer>" in line:
				newfile.write(line)
				switch = 0
			else:
				for newline in newSideBar:
					#print newline
					newfile.write(newline)
	openfile.close()
	newfile.close()


