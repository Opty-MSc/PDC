all:
	cd Serial && make
	cd OMP && make
	cd MPI && make

clean:
	@echo Cleaning...
	cd Serial && make $@
	cd OMP && make $@
	cd MPI && make $@
