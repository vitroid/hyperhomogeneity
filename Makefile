########## Preparations

prepare:
	pip install pairlist vapory yaplotlib f2py matplotlib numpy networkx

fmodules: $(patsubst %.f95, %.cpython-37m-darwin.so, $(wildcard *.f95))
%.cpython-37m-darwin.so: %.f95
	f2py3 -c $< -m $*

########## General Rules

# Determine the weights by the quasi-inverse matrix method.
%.cycles5.pickle: q/%.q.nx3a
	python cycles5.py $< $@

# Record all molecule-to-cycle interactions.
%.cyclesintr.pickle: q/%.q.nx3a %.cycles5.pickle
	python cyclesintr.py $^ $@

# Classify the molecular arrangements in the cell.
# Used in Figure S5 only
%.repr.pickle: r/%.nx3a
	python unique_orientations3.py $< $@

# Molecule-cycle interaction classified by the molecular arrangements.
# Used in Figure S5 only.
%.cycles5stat.pickle: q/%.q.nx3a %.repr.pickle %.cycles5.pickle
	python3 cycles5stat.py $^ $@


ices.%:
	for ice in 1h 3 5 6 7; do ls q/$$ice-*.q.nx3a | sed -e "s/q.nx3a/$*/g" -e "s:q/::g"; done | xargs make -k -j 6
extend.%: $(wildcard q/*.q.nx3a)
	echo $^ | sed -e "s/q.nx3a/$*/g" -e "s:q/::g" | xargs make -k -j 6


%.couhi.pickle: q/%.q.nx3a
	python CoulombHist.py $< $@

########## Everything

everything: ices.hist.pickle ices.cycles5.pickle ices.cyclesintr.pickle ices.cycles5stat.pickle ices.repr.pickle H2.pickle ices.couhi.pickle Figure1.pdf Figure2a.svg Figure2bc.svg Figure3.pdf Figure4.pdf FigureS2.yap FigureS3.pdf FigureS4.pdf FigureS5.pdf FigureS6.pdf
	echo Done.


########## Figures

# Histogram of cumulative interactions on a molecule basis
%.hist.pickle: q/%.q.nx3a
	python histogram.py $< $@

# Figure 1: Overview of the phenomenon.
Figure1.pdf: Figure1.py
	-make ices.hist.pickle
	python Figure1.py

Figure2a.svg: Figure2a.py
	python Figure2a.py

Figure2bc.svg: Figure2bc.py
	python Figure2bc.py

# Figure 3: Comparison between molecule and cycle basis
Figure3.pdf: Figure3.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	python Figure3.py

# Figure 4: Interactions with cycles.
Figure4.pdf: Figure4.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle ices.couhi.pickle
	python Figure4.py

# Illustration of cycles that passes through a water molecule in ice Ih.
# Render with yaplot (https://github.com/vitroid/Yaplot)
FigureS2.yap: FigureS2.py
	-make 1h-1000.cycles5.pickle q/1h-1000.q.nx3a
	python FigureS2.py

# Interaction of a dipole with cycles that pass through the dipole,
# classified by the size of the cycle.
FigureS3.pdf: FigureS3.py
	-make ices.cycles5.pickle ices.cyclesintr.pickle
	python FigureS3.py

# Divergence of the interaction at the ice surface.
FigureS4.pdf: 1cs.nx3a FigureS4.py
	python FigureS4.py

# Cumulative Molecule-cycle interaction classified by the molecular arrangements.
FigureS5.pdf: FigureS5.py
	-make extend.cycles5.pickle extend.repr.pickle extend.cycles5stat.pickle 9.cycles5stat.pickle
	python FigureS5.py

# Cumulative interaction of CO-like molecule
FigureS6.pdf: FigureS6.py
	python FigureS6.py
