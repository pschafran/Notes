#! /bin/bash

ifconfig en0 ether 00:e2:e3:e4:e5:e6
DATE=`date +%Y-%m-%d`
dropbox="Dropbox_"
herbarium="ODUHerbarium_"
drafts="Drafts_"
onedrive="OneDriveODU_"
genbank="GenBank_"
isoetesdocs="WorkingIsoetesDocs_"
lilabbti="Li-Lab-BTI_"
lilabbtiwiki="Li-Lab-BTI.wiki_"
#print "Beginning backup, do not shut down until done!"
#print "**********************************************"
#print "Backing up Dropbox..."
tar -czf /Users/Peter/Desktop/Backups/${dropbox}${DATE}.tar.gz /Users/Peter/Dropbox\ \(Personal\)/
#print "Backing up Herbarium..."
tar -czf /Users/Peter/Desktop/Backups/${herbarium}${DATE}.tar.gz /Users/Peter/Dropbox\ \(Personal\)/ODU\ Herbarium/
#print "Backing up Drafts..."
tar -czf /Users/Peter/Desktop/Backups/${drafts}${DATE}.tar.gz /Users/Peter/Desktop/Drafts/
#print "Backing up OneDrive..."
tar -czf /Users/Peter/Desktop/Backups/${onedrive}${DATE}.tar.gz /Users/Peter/OneDrive\ -\ Old\ Dominion\ University/
#print "Done!"
tar -czf /Users/Peter/Desktop/Backups/${genbank}${DATE}.tar.gz /Users/Peter/Desktop/GenBank/
tar -czf /Users/Peter/Desktop/Backups/${isoetesdocs}${DATE}.tar.gz /Users/Peter/Documents/Working_Isoetes_Docs/
tar -czf /Users/Peter/Desktop/Backups/${lilabbti}${DATE}.tar.gz /Users/Peter/Li-Lab-BTI/
tar -czf /Users/Peter/Desktop/Backups/${lilabbtiwiki}${DATE}.tar.gz /Users/Peter/Li-Lab-BTI.wiki/