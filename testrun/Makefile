.PHONY: all
all: run-stamp pwgNLO.top plots-stamp

run-stamp: powheg.input-save
	./run.sh 0
	touch run-stamp

pwgNLO.top: run-stamp
	(echo 1 ; ls -c1 pwg-*-NLO.top ; echo "") | ../shell/mergedata
	mv fort.12 pwgNLO.top

plots-stamp: pwgNLO.top
	cd plots; ../../shell/genplots.sh ../pwgNLO.top ../pwgNLO.top NLO; cd -;
	cd plots; ../../shell/gnuplotsplit.sh genplots.gp; cd -;
	touch plots-stamp

.PHONY: clean veryclean

# Remove all files related to plotting
clean:
	rm -f plots-stamp
	rm -f plots/*.*
	rm -f pwgNLO.top

# Remove all POWHEG/plotting generated files
veryclean: clean
	rm -f run-stamp
	rm -fr pwg*stat.dat
	rm -fr pwg*grid.top
	rm -fr BadPoints
	rm -f pwgcounters-st*
	rm -f run-*-*.log
	rm -f run-*-*-*.log
	rm -f pwggrid*
	rm -f pwggridinfo-btl-xg*
	rm -f powheg.input-2
	rm -f fort.7
	rm -f FlavRegList
	rm -f pwhg_checklimits
	rm -f pwgboundviolations-*.dat
	rm -f pwgubound-*.dat
	rm -f pwgxgrid.dat
	rm -f pwghistnorms.top
	rm -f pwgborngrid.top
	rm -f pwgcounters.dat
	rm -f pwg-*-NLO.top
	rm -f *equiv*
	ls powheg.input* | grep -v '\-save' |  xargs rm -f
	rm -f log-*.log
	rm -f run-*.log
	rm -f pwgcounters0*.dat
	rm -f pwgLHEF_analysis-*-W*.top
	rm -f pwgLHEF_analysis-max.top
	rm -f pwgLHEF_analysis-min.top
	rm -f pwgLHEF_analysis.top
	rm -f pwgLHEF_analysis-W*.top
	rm -f pwgbtildeupb-*.dat
	rm -f pwgremnupb-*.dat
	rm -f pwgevents*.lhe
	rm -f pwgfullgrid*.dat
	rm -f testbnd*.top
	rm -f mint_*.top
	rm -f pwgpwhgalone-*.top
	rm -f lock *.lock
