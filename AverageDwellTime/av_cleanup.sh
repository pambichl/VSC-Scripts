if [ ! -d wavefuncs ]; then
 mkdir wavefuncs
fi
mv pic.*.streu.*.jpg wavefuncs/

if [ ! -d potentials ]; then
 mkdir potentials
fi
mv pic.potential.*.jpg potentials/

if [ ! -d submitfiles ]; then
 mkdir submitfiles
fi
mv submit.*.sh submitfiles/
mv submit.*.sh.* submitfiles/

if [ ! -d inputfiles ]; then
 mkdir inputfiles
fi
mv input.*.xml inputfiles/

if [ ! -d scatterdata ]; then
 mkdir scatterdata
fi
mv Smat.*.dat scatterdata

if [ ! -d refindices ]; then
 mkdir refindices
fi
mv obstacles.*.dat refindices

rm TR.*.dat
rm TR.*.txt
rm eigen.*.log
