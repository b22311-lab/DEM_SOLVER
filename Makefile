FC ?= gfortran
FFLAGS ?= -O3 -Wall -Wextra
OMPFLAGS ?= -fopenmp

EXE := demo.exe

OBJS := \
	input_reader.o \
	domain_mkr.o \
	particle_init.o \
	write_output.o \
	main_loop.o

.PHONY: all clean clean-build clean-artifacts run

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $(OMPFLAGS) -o $@ $(OBJS)

input_reader.o: input_reader.f90
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $<

domain_mkr.o: domain_mkr.f90 input_reader.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $<

particle_init.o: particle_init.f90 input_reader.o domain_mkr.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $<

write_output.o: write_output.f90 input_reader.o particle_init.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $<

main_loop.o: main_loop.f90 input_reader.o domain_mkr.o particle_init.o write_output.o
	$(FC) $(FFLAGS) $(OMPFLAGS) -c $<

run: $(EXE)
	./$(EXE) inputs/input_config.in

clean-build:
	rm -f *.o *.mod $(EXE)

clean-artifacts:
	rm -f output_*.vtk
	rm -f *.png *.csv *.log
	rm -f *.aux *.out *.toc *.lof *.lot *.fls *.fdb_latexmk missfont.log
	rm -f *.bbl *.blg *.synctex.gz ieee_dem_whitepaper.pdf
	rm -f part18_findings.tex part19_2_findings.tex
	rm -f strong_scaling_findings.tex weak_scaling_findings.tex
	rm -rf perf_logs part10_outputs part12_profile
	rm -rf part18_logs part18_outputs part19_configs part19_outputs
	rm -rf parallel_validation_outputs strong_scaling_logs weak_scaling_logs weak_scaling_configs
	rm -rf latex_build

clean:
	$(MAKE) clean-build
	$(MAKE) clean-artifacts
