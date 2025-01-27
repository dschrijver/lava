CHOSEN_TESTCASES = poiseuille-body poiseuille-pressure stefan channel poiseuille-dual


all: poiseuille-body poiseuille-pressure stefan channel poiseuille-dual

selection: $(CHOSEN_TESTCASES)


poiseuille-body:
	make -C Poiseuille_Body_Force test

poiseuille-pressure:
	make -C Poiseuille_Pressure_Boundary_Conditions test

stefan:
	make -C Stefan test

channel:
	make -C Channel test

 poiseuille-dual:
	make -C Poiseuille_Dual_Component test
