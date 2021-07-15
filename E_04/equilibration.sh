#make
cp config.fcc config.0
./MolDyn_NVE.exe 0 1 0
for i in {1..7}; do cp config.final config.0; cp old.final old.0; ./MolDyn_NVE.exe 1 1 0; done

# <./MolDyn_NVE.exe> <readold> <rescale> <average>
