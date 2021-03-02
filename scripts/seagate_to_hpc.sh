#! /bin/bash

rsync -vr --ignore-existing --progress /Volumes/Seagate_Blue/School/Backups/ pscha005@turing.hpc.odu.edu:/RC/home/pscha005/Backups/
rsync -vr --ignore-existing --progress /Volumes/Seagate_Blue/Photos/ pscha005@turing.hpc.odu.edu:/RC/home/pscha005/Data/Photos
