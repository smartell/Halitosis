#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=vonBH
TPL=$(EXEC).tpl
DAT=$(EXEC).dat
ARG=


all: $(EXEC) $(EXEC).rep dust

$(EXEC): $(TPL)
	admb $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT) $(ARG)

mcmc: runmc mceval dust

runmc: $(EXEC) $(EXEC).rep
	./$(EXEC) -ind $(DAT) -mcmc 100000 -mcsave 250 -nosdmcmc -mcscale 50000

mceval: $(EXEC) $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval



dust:
	rm -f *.log *.rpt *.htp admodel.* variance *.bar *.mcm *.[prb][0123]*

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.p* *.b*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance