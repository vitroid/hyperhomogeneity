########## Preparations
BASE=/Volumes/workarea/venvs/hyperhomogeneity/bin#/net/jukebox4.local/u2/matto/venvs/hyperhomogeneity/bin
# BASE=/Users/matto/miniforge3/bin
BASE=/Users/matto/venvs/genice2/bin
PIP=$(BASE)/pip3
PYTHON=$(BASE)/python3
INKSCAPE=/Applications/Inkscape.app/Contents/MacOS/inkscape

prepare:
	$(PIP) install pairlist vapory yaplotlib matplotlib numpy networkx scipy cycless
	pip install --upgrade pip

fmodules: $(patsubst %.f95, %.cpython-37m-darwin.so, $(wildcard *.f95))
%.cpython-37m-darwin.so: %.f95
	f2py3 -c $< -m $*

########## General Rules

# Determine the weights by the quasi-inverse matrix method.
%.cycles5.pickle: q/%.q.nx3a
	$(PYTHON) cycles5.py $< $@

# Record all molecule-to-cycle interactions.
%.cyclesintr.pickle: q/%.q.nx3a %.cycles5.pickle
	$(PYTHON) cyclesintr.py $^ $@

# Classify the molecular arrangements in the cell.
# Used in Figure S5 only
%-1000.repr.pickle: r/%-1000.nx3a
	$(PYTHON) unique_orientations3.py $< $@
# Share the same table to use the unique labelings to the arrangements
	i=1001; while [ $$i -lt 1030 ]; do ln $*-1000.repr.pickle $*-$$i.repr.pickle; i=`expr $$i + 1`; done

# Molecule-cycle interaction classified by the molecular arrangements.
# Used in Figure S5 only.
%.cycles5stat.pickle: q/%.q.nx3a %.repr.pickle %.cycles5.pickle
	$(PYTHON) cycles5stat.py $^ $@

%.png: %.pdf
	$(INKSCAPE) $< -o $@
pngs:
	ls *.pdf | sed -e s/pdf/png/ | xargs make -j 4 -k

#ices.%:
#	for ice in 1h 3 5 6 7; do ls q/$$ice-1000.q.nx3a | sed -e "s/q.nx3a/$*/g" -e "s:q/::g"; done | xargs make -j 32 -k
ices.%:
	for ice in 1h 3 5 6 7; do ls *-*.cycles5.pickle | sed -e "s/cycles5.pickle/$*/g" ; done | xargs make -j 32 -k
#extend.%: $(wildcard *-1000.q.nx3a)
#	echo $^ | sed -e "s/q.nx3a/$*/g" -e "s:q/::g" | xargs make -k -j 32
#extendr.%: $(wildcard r/*-1000.nx3a)
#	echo $^ | sed -e "s/nx3a/$*/g" -e "s:r/::g" | xargs make -k -j 32
extend.%: $(wildcard *-*.cycles5.pickle)
	echo $^ | sed -e "s/cycles5.pickle/$*/g" | xargs make -k -j 32
extendr.%:
	make extend.$*

%.couhi.pickle: q/%.q.nx3a CoulombHist.py
	$(PYTHON) CoulombHist.py $< $@

########## Everything

everything: ices.hist.pickle extend.cycles5.pickle ices.cyclesintr.pickle \
	extend.cycles5stat.pickle extendr.repr.pickle ices.couhi.pickle \
	Figure1.pdf Figure2a.svg Figure2bc.svg Figure3.pdf Figure4.pdf FigureS2.yap \
	FigureS3.pdf FigureS4.pdf FigureS5.pdf FigureS1.pdf
	echo Done.


########## Figures

# Histogram of cumulative interactions on a molecule basis
%.hist.pickle: q/%.q.nx3a
	$(PYTHON) histogram.py $< $@

# Interaction distribution with nearby molecules
FigureA.pdf: FigureA.py # q/11.q.nx3a
	-make ices.hist.pickle
	$(PYTHON) FigureA.py

FigureE.pdf: $(wildcard q/1cs*.q.nx3a) FigureE.py
	$(PYTHON) FigureE.py


# Figure 1: Overview of the phenomenon.
Figure1.pdf: Figure1.py q/11.q.nx3a
	-make ices.hist.pickle
	$(PYTHON) Figure1.py

Figure2a.svg: Figure2a.py
	$(PYTHON) Figure2a.py

Figure2bc.svg: Figure2bc.py
	$(PYTHON) Figure2bc.py

# Figure 3: Comparison between molecule and cycle basis
Figure3.pdf: Figure3.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	$(PYTHON) Figure3.py

# Figure 4: Divergence of the interaction at the ice surface.
Figure4.pdf: q/1cs.q.nx3a Figure4.py
	$(PYTHON) Figure4.py


# Cumulative interaction of CO-like molecule
FigureS1.pdf: FigureS1.py
	$(PYTHON) FigureS1.py

# Illustration of cycles that passes through a water molecule in ice Ih.
# Render with yaplot (https://github.com/vitroid/Yaplot)
FigureS2.yap: FigureS2.py
	-make 1h-1000.cycles5.pickle q/1h-1000.q.nx3a
	$(PYTHON) FigureS2.py

# Interaction of a dipole with cycles that pass through the dipole,
# classified by the size of the cycle.
FigureS3.pdf: FigureS3.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	$(PYTHON) FigureS3.py

# Interactions with cycles.
FigureS4.pdf: FigureS4.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	$(PYTHON) FigureS4.py

# Cumulative Molecule-cycle interaction classified by the molecular arrangements.
FigureS5.pdf: FigureS5.py
	-make extend.cycles5.pickle extendr.repr.pickle extend.cycles5stat.pickle # 9.cycles5stat.pickle
	$(PYTHON) FigureS5.py

# Figure S7: Artifact on SD in a very small cell at the ice surface.
FigureS7.pdf: FigureS7.py
	$(PYTHON) FigureS7.py

########## Sync
sync:
	rsync -av q r *.repr.pickle bluebird1.local:/r7/matto/hyperhomogeneity/
syncback:
	rsync -av --include="*/" --include="*.pdf" --include="*.svg" --include="*.yap" --exclude="*" bluebird1.local:/r7/matto/hyperhomogeneity/* /Volumes/workarea/work/hyperhomogeneity
