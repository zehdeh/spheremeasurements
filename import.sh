#! /bin/sh

for i in `seq -f "%02g" 1 22`; do
	./cleanup_mesh.py ~/test/Scans/FullVolium/single_cam_reconstructions/$i/ res/final/onecam_fullvolume/$i/cleanedup/
	./eliminate_badscans.py res/final/onecam_fullvolume/$i/cleanedup/ --verbose
done
