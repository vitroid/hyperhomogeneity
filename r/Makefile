all: 11.nx3a
	i=1000; while [ $$i -lt 1030 ]; do for ice in 1h 3 5 6 7 4 1c 16 CS1 2D2 1cs ; do echo $$ice-$$i.nx3a; done; i=`expr $$i + 1`; done | xargs make -j 8 -k

GENICE=genice2

1h-%.nx3a:
	$(GENICE) 1h -s $* -r 4 4 4 -f e > $@
11.nx3a:
	$(GENICE) 11 -r 4 4 4 -f e > $@
3-%.nx3a:
	$(GENICE) 3 -s $* -r 4 4 4 -f e > $@
5-%.nx3a:
	$(GENICE) 5R -s $* -r 3 4 1 -f e > $@
6-%.nx3a:
	$(GENICE) 6 -s $* -r 4 4 4 -f e --dens 1.389 > $@
7-%.nx3a:
	$(GENICE) 7 -s $* -r 4 4 4 -f e > $@
4-%.nx3a:
	$(GENICE) 4R -s $* -r 2 3 2 -f e > $@
1c-%.nx3a:
	$(GENICE) 1c -s $* -r 5 5 5 -f e > $@
16-%.nx3a:
	$(GENICE) 16 -s $* -r 2 2 2 -f e > $@
CS1-%.nx3a:
	$(GENICE) CS1 -s $* -r 3 3 3 -f e > $@
9.nx3a:
	$(GENICE) 9 -r 4 4 4 -f e > $@
2D2-%.nx3a:
	$(GENICE) 2D2 -s $* -r 8 8 1 -f e > $@
1cs-%.nx3a:
	$(GENICE) 1c -s $* -r 4 4 8 -f e | awk '/BOX3/{j=0}(j==1){$$3=100}{print;j++}'> $@
1cxs-%.nx3a:
	$(GENICE) 1c -s $* -r 2 2 2 -f e | ./double.py > $@
