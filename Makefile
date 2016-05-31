all: none

none:


.PHONY: export


VERSIONSTR = `git rev-parse HEAD | cut -c 1-8`
VERSIONSTRF = `git rev-parse HEAD`
TARFILE = solarv-$(VERSIONSTR)

export: 
	git archive --prefix=solarv/ HEAD > $(TARFILE).tar
	bzip2 $(TARFILE).tar
