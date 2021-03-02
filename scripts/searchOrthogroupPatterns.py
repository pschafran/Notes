#! /usr/bin/env python

import sys

stomataPresence = ["Anthoceros_agrestis","Anthoceros_angustus","Anthoceros_fusiformis","Anthoceros_punctatus","Leiosporoceros_dussii","Paraphymatoceros_pearsonii","Phaeoceros_carolinianus","Phaeoceros_sp","Phymatoceros_bulbiculosus"]
symbiosisPresence = ["Anthoceros_agrestis","Anthoceros_angustus","Anthoceros_fusiformis","Anthoceros_punctatus","Notothylas_orbicularis","Phaeoceros_carolinianus","Phaeoceros_sp","Phymatoceros_bulbiculosus"]
ccmPresence = ["Anthoceros_agrestis","Anthoceros_angustus","Anthoceros_fusiformis","Anthoceros_punctatus","Notothylus_orbicularis","Phaeoceros_carolinianus","Phaeoceros_sp","Phymatoceros_bulbiculosus"]


stomataOutfile = open("OGs_matching_stomata_pattern.txt", "a")
symbiosisOutfile = open("OGs_matching_symbiosis_pattern.txt", "a")
ccmOutfile = open("OGs_matching_CCM_pattern.txt", "a")

with open(sys.argv[1], "r") as infile:
	speciesPresent = []
	for line in infile:
		if line.startswith(">"):
			splitline = line.strip(">\n").split("_")
			species = "%s_%s" %(splitline[0], splitline[1])
			speciesPresent.append(species)
	if set(speciesPresent) == set(stomataPresence):
		stomataOutfile.write("%s\n" % sys.argv[1])
	if set(speciesPresent) == set(symbiosisPresence):
		symbiosisOutfile.write("%s\n" % sys.argv[1])
	if set(speciesPresent) == set(ccmPresence):
		ccmOutfile.write("%s\n" % sys.argv[1])
stomataOutfile.close()
symbiosisOutfile.close()
ccmOutfile.close()
