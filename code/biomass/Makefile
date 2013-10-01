.PHONY : tables
.PHONY : extract

#Default variable values (in case none were specified on the command line):
outdir ?= ~/git/paleon/output/sept-2013
root ?= ~/misc/paleon
cluster ?= 18



#extract the tarred-and-gzipped simulation results, import them to R, and save the resulting Rdata.
extract: 
	mkdir -p $(root)/output
	sh extract.sh $(root) $(cluster)

