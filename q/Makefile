all:
	i=1000; while [ $$i -lt 1030 ]; do for ice in 1h 3 5 6 7 4 1c 16 CS1 2D2 1cs ; do echo $$ice-$$i.q.nx3a; done; i=`expr $$i + 1`; done | xargs make -j 8 -k

ANALICE=analice2

%.q.nx3a: ../r/%.nx3a
	$(ANALICE) $< -w ice -f q | cat DEFR RCOA13 ID08TIP4Pice - | ~/GoogleDrive/gitwork/simplemd/simplequench 1000 | python3 nx4a2nx3a.py > $@
