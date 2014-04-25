#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=HGM
TPL=$(EXEC).tpl
DAT=LenAgeKey2011.dat


all: $(EXEC) $(EXEC).rep dust

$(EXEC): $(TPL)
	admb $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT)

mcmc: runmc evalmc dust

runmc: $(EXEC) $(EXEC).rep
	./$(EXEC) -ind $(DAT) -mcmc 500000 -mcsave 200 -nosdmcmc

evalmc: $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval



dust:
	rm -f *.log *.rpt *.htp admodel.* variance *.bar *.mcm *.[prb][0123456789]*

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.p* *.b*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance