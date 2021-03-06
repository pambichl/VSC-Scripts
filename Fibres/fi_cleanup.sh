if [ ! -d peakwavesfuncs ]; then
 mkdir peakwavefuncs
fi
mv pic.*.coeff.*.jpg peakwavefuncs/

if [ ! -d coefffiles ]; then
 mkdir coefffiles
fi
mv coeff.*.dat coefffiles/

if [ ! -d inputfiles ]; then
 mkdir inputfiles
fi
mv qpeak.*.xml inputfiles/
mv ypeak.xml inputfiles/
mv tdtpeak.*.*.xml inputfiles/
mv tdtpeak.xml inputfiles/

if [ ! -d submitfiles ]; then
 mkdir submitfiles
fi
mv submit.sh.* submitfiles/
mv qpeaksubmit.*.sh submitfiles/
mv qpeaksubmit.*.sh.* submitfiles/
mv ypeaksubmit.sh submitfiles/
mv ypeaksubmit.sh.* submitfiles/
mv tdtpeaksubmit.*.*.sh submitfiles/
mv tdtpeaksubmit.*.*.sh.* submitfiles/

if [ ! -d meigenstates ]; then
 mkdir meigenstates
fi
mv pic.*.streu.*.jpg meigenstates/

rm TR.*.dat
rm TR.*.txt
rm eigen.*.log
rm N_states
rm evals.txt
rm evecs.txt
rm fort.24
rm all_evals.txt
rm *.p
