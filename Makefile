#################
# Interactome3d
#################

#Organism
ORGANISM = human

#Main directories
SRCDIR = $(CURDIR)/src
DATADIR = $(CURDIR)/data
ORGFILE = $(CURDIR)/organisms.tab
RESULTSDIR = $(CURDIR)/results

#Tools
NACCESS ?= $(shell which naccess)
WGET ?= $(shell which wget)
SED ?= $(shell which gsed || which sed)
CURL ?= $(shell which curl)
EGREP ?= $(shell which egrep)
SORT ?= $(shell which sort)
UNIQ ?= $(shell which uniq)
TAR ?= $(shell which tar)
PERL ?= $(shell which perl)
RSCRIPT ?= $(shell which Rscript) 

#Uniprot
UNIPROTFTP ?= ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/
UNIPROTFASTA = $(DATADIR)/uniprot_$(ORGANISM).fasta

#Interactome3d
INT3URL ?= 'http://interactome3d.irbbarcelona.org/user_data/$(ORGANISM)/download/representative'
FILELISTPAGE ?= 'http://interactome3d.irbbarcelona.org/downloadset.php?queryid=$(ORGANISM)&release=current&path=representative'

# Necessary for variables parsing
NACCESSCUT = $(shell echo $1 | sed -r 's/([^\.]+).*/\1/')
PARENTPDB = $(shell echo $1 | sed -r 's/(.+)_[AB].pdb/\1/')
CSVCUT = $(shell sed "1d" $(ORGFILE) | grep -E '^$(1)\s' | cut -f$(2))

#List of tgz to download
FILELIST ?= $(DATADIR)/filelist.txt
#Interactions
INTERACTOME3DINTSDIR ?= $(DATADIR)/interactions
INTERACTOME3DINTSTARDIR ?= $(INTERACTOME3DINTSDIR)/tars
INTERACTOME3DINTSPDBDIR ?= $(INTERACTOME3DINTSDIR)/pdbs
INTERACTOME3DINTSRSADIR ?= $(INTERACTOME3DINTSDIR)/rsas
INTERACTOME3DINTSLOGDIR ?= $(INTERACTOME3DINTSDIR)/logs
INTERACTOME3DINTSASADIR ?= $(INTERACTOME3DINTSDIR)/asas
INTERACTOME3DINTFILE ?= $(DATADIR)/interactions.dat
#Proteins
INTERACTOME3DPROTFILE ?= $(DATADIR)/proteins.dat
INTERACTOME3DPROTSDIR ?= $(DATADIR)/proteins
INTERACTOME3DPROTSPDBDIR ?= $(INTERACTOME3DPROTSDIR)/pdbs
INTERACTOME3DPROTSRSADIR ?= $(INTERACTOME3DPROTSDIR)/rsas
INTERACTOME3DPROTSLOGDIR ?= $(INTERACTOME3DPROTSDIR)/logs
INTERACTOME3DPROTSASADIR ?= $(INTERACTOME3DPROTSDIR)/asas

#Output
INTERFACESFILE ?= $(RESULTSDIR)/interfaces.tab
INTERFACESFILE_NOTREMAPPED ?= $(RESULTSDIR)/interfaces_notmapped.tab
INTERACTPROTACC ?= $(RESULTSDIR)/accessibilities.tab

#pdb mappings
PDBPOSMAPPINGS ?= $(DATADIR)/pdb_mappings.rds

#Interactome3d
INTERACTOME3DINTSTARFILES = $(foreach FILE,$(shell grep interactions $(FILELIST) | grep tgz),$(INTERACTOME3DINTSTARDIR)/$(basename $(FILE)).tgz)
INTERACTOME3DINTPDBS = $(shell grep pdb $(INTERACTOME3DINTFILE) | cut -f 22)
INTERACTOME3DINTPDBSFILENAMES = $(addprefix $(INTERACTOME3DINTSPDBDIR)/,$(INTERACTOME3DINTPDBS))
INTERACTOME3DINTRSAS = $(INTERACTOME3DINTPDBS:%.pdb=%.rsa)
INTERACTOME3DINTRSASFILENAMES = $(addprefix $(INTERACTOME3DINTSRSADIR)/,$(INTERACTOME3DINTRSAS))

INTERACTOME3DPROTPDBS = $(addsuffix _A.pdb, $(basename $(shell grep pdb $(INTERACTOME3DINTFILE) | cut -f 22))) $(addsuffix _B.pdb,$(basename $(shell grep pdb $(INTERACTOME3DINTFILE) | cut -f 22)))
INTERACTOME3DPROTRSAS = $(INTERACTOME3DPROTPDBS:%.pdb=%.rsa)
INTERACTOME3DPROTRSASFILENAMES = $(addprefix $(INTERACTOME3DPROTSRSADIR)/,$(INTERACTOME3DPROTRSAS))

.PRECIOUS: $(FILELIST) $(INTERACTOME3DINTFILE) $(INTERACTOME3DINTRSAS)
.PHONY : int3Ddirectories clean-download int3Drsas

prepare: $(UNIPROTFASTA) download-indexes

all: $(INTERFACESFILE_NOTREMAPPED) $(INTERFACESFILE) $(INTERACTPROTACC)

download-tars: $(INTERACTOME3DINTSTARFILES)

download-indexes: int3Ddirectories $(FILELIST) $(INTERACTOME3DINTFILE) 

test: $(UNIPROTFASTA)

int3Ddirectories:
	mkdir -p $(DATADIR)
	mkdir -p $(RESULTSDIR)
	mkdir -p $(INTERACTOME3DPROTSDIR)
	mkdir -p $(INTERACTOME3DPROTSPDBDIR)
	mkdir -p $(INTERACTOME3DPROTSRSADIR)
	mkdir -p $(INTERACTOME3DPROTSLOGDIR)
	mkdir -p $(INTERACTOME3DPROTSASADIR)
	mkdir -p $(INTERACTOME3DINTSDIR)
	mkdir -p $(INTERACTOME3DINTSPDBDIR)
	mkdir -p $(INTERACTOME3DINTSRSADIR)
	mkdir -p $(INTERACTOME3DINTSTARDIR)
	mkdir -p $(INTERACTOME3DINTSLOGDIR)
	mkdir -p $(INTERACTOME3DINTSASADIR)

int3Drsas: download-tars $(INTERACTOME3DINTRSASFILENAMES) $(INTERACTOME3DPROTRSASFILENAMES)

clean:
	-rm $(INTERACTOME3DINTFILE)
	-rm -r $(FILELIST)
	-rm -rvf $(DATADIR)/*
	-rm $(INTERFACESFILE)
	-rm $(INTERACTPROTACC)

#Download fasta from uniprot
$(UNIPROTFASTA):
	$(WGET) -P $(DATADIR) $(UNIPROTFTP)$(call CSVCUT,$(ORGANISM),2)/$(call CSVCUT,$(ORGANISM),3).fasta.gz
	gunzip $(DATADIR)/$(call CSVCUT,$(ORGANISM),3).fasta.gz
	mv $(DATADIR)/$(call CSVCUT,$(ORGANISM),3).fasta $@

#Download list of files on the 
$(FILELIST):
	$(WGET) -O - $(FILELISTPAGE) | $(EGREP) -o 'proteins_[0-9]+\.tgz|interactions_[0-9]+\.tgz' | $(SORT) | $(UNIQ)  > $@

#Download interactions.dat
$(INTERACTOME3DINTFILE):
	$(WGET) $(INT3URL)/interactions.dat -O $@

# #Downloading tgz
$(INTERACTOME3DINTSTARDIR)/%.tgz:
	$(WGET) $(INT3URL)/$*.tgz -O $@
	$(TAR) -xzf $@ -C $(INTERACTOME3DINTSPDBDIR)

#Creating interaction rsa using naccess
$(INTERACTOME3DINTSRSADIR)/%.rsa: $(INTERACTOME3DINTSPDBDIR)/%.pdb
	-$(NACCESS) $(INTERACTOME3DINTSPDBDIR)/$*.pdb
	if test -f $(call NACCESSCUT,$*).rsa; then \
		mv $(call NACCESSCUT,$*).rsa $(INTERACTOME3DINTSRSADIR)/$*.rsa ;\
		mv $(call NACCESSCUT,$*).log $(INTERACTOME3DINTSLOGDIR)/$*.rsa ;\
		mv $(call NACCESSCUT,$*).asa $(INTERACTOME3DINTSASADIR)/$*.rsa ;\
	fi

#Individual proteins extracted from interactome3d complexes
$(INTERACTOME3DPROTSPDBDIR)/%.pdb: 
	$(PERL) $(SRCDIR)/splitInteraction.pl $(INTERACTOME3DINTSPDBDIR)/$(call PARENTPDB, $*.pdb).pdb $(INTERACTOME3DPROTSPDBDIR)

#Creating individual protein rsa using naccess
$(INTERACTOME3DPROTSRSADIR)/%.rsa: $(INTERACTOME3DPROTSPDBDIR)/%.pdb
	-$(NACCESS) $(INTERACTOME3DPROTSPDBDIR)/$*.pdb
	if test -f $(call NACCESSCUT,$*).rsa; then \
		mv $(call NACCESSCUT,$*).rsa $(INTERACTOME3DPROTSRSADIR)/$*.rsa ;\
		mv $(call NACCESSCUT,$*).log $(INTERACTOME3DPROTSLOGDIR)/$*.rsa ;\
		mv $(call NACCESSCUT,$*).asa $(INTERACTOME3DPROTSASADIR)/$*.rsa ;\
	fi

$(INTERFACESFILE_NOTREMAPPED): int3Drsas
	$(PERL) $(SRCDIR)/interactionsParser.pl \
	$(INTERACTOME3DINTFILE) \
	$(DATADIR) > $@


$(PDBPOSMAPPINGS): $(INTERACTOME3DINTFILE)
	$(RSCRIPT) $(SRCDIR)/getPdbPositions.R \
	$(INTERACTOME3DINTFILE) \
	$(INTERACTOME3DINTSPDBDIR) \
	$(UNIPROTFASTA) \
	$@

$(INTERFACESFILE): $(PDBPOSMAPPINGS)
	$(RSCRIPT) $(SRCDIR)/rempaPdbPositions.R \
	$(INTERFACESFILE_NOTREMAPPED) \
	$(PDBPOSMAPPINGS) \
	$@

#this file is not remapped
$(INTERACTPROTACC): int3Drsas
	$(PERL) $(SRCDIR)/intProteinsAccParser.pl \
	$(INTERACTOME3DINTFILE) \
	$(DATADIR) > $@
