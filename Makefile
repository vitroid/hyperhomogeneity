########## Preparations

# On my mac
BASE=/Users/matto/venvs/genice2/bin/
INKSCAPE=/Applications/Inkscape.app/Contents/MacOS/inkscape
CORES=8
# On linux
BASE=
CORES=32

PIP=$(BASE)pip3
PYTHON=$(BASE)python3

prepare:
	$(PIP) install pairlist vapory yaplotlib matplotlib numpy networkx scipy svgwrite cycless
	$(PIP) install --upgrade pip

fmodules:
	for f in *.f95; do f2py3 -c $$f -m `basename -s .f95 $$f`; done

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

prep: $(wildcard q/*[^L]-1???.q.nx3a)
	ls $^ | sed -e s@^q/@@ -e 's/.q.nx3a/.cycles5.pickle/' | xargs make -k #-j$(CORES) -k
ices.%:
	for ice in 1h 3 5 6 7; do ls q/$$ice-10??.q.nx3a | sed -e "s/q.nx3a/$*/g" -e "s:q/::g"; done | xargs make -j $(CORES) -k
icex.%:
	for ice in 1h 3 5 6 7 16 CS1 2D2 1c ; do ls q/$$ice-10??.q.nx3a | sed -e "s/q.nx3a/$*/g" -e "s:q/::g"; done | xargs make -j $(CORES) -k
icex1.%:
	for ice in 1h 3 5 6 7 16 CS1 2D2 1c ; do ls q/$$ice-10??.q.nx3a | sed -e "s/q.nx3a/$*/g" -e "s:q/::g"; done | xargs make -k

%.couhi.pickle: q/%.q.nx3a CoulombHist.py
	$(PYTHON) CoulombHist.py $< $@

########## Everything

everything:
	make -k fmodules ices.hist.pickle icex1.cycles5.pickle ices.cyclesintr.pickle \
	icex.cycles5stat.pickle icex.repr.pickle ices.couhi.pickle \
	Figure1.pdf Figure3a.svg Figure3bc.svg Figure4.pdf FigureS2.yap \
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

Figure3a.svg: Figure3a.py
	$(PYTHON) Figure3a.py

Figure3bc.svg: Figure3bc.py
	$(PYTHON) Figure3bc.py

# Figure 3: Comparison between molecule and cycle basis
Figure4.pdf: Figure4.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	$(PYTHON) Figure4.py

# Figure 4: Divergence of the interaction at the ice surface.
Figure5.pdf: $(wildcard q/1cs*.q.nx3a) Figure5.py
	$(PYTHON) Figure5.py


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


########## Remote job
REMOTE=172.23.78.15
sync:
	rsync -av *.f95 sd_ice r q *.py Makefile $(REMOTE):/r7/matto/hyperhomogeneity/
	#rsync -av --include="*/" --exclude="*.nx3a" *.f95 sd_ice r q *.py Makefile $(REMOTE):/r7/matto/hyperhomogeneity/
syncback:
	rsync -av --include="*/" --include="*.pdf" --include="*.svg" --include="*.yap" --exclude="*" $(REMOTE):/r7/matto/hyperhomogeneity/* /Volumes/workarea/work/hyperhomogeneity
