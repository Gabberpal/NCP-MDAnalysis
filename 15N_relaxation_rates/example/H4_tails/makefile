SHELL := /bin/bash

.PHONY: all clean

default: all

SUBDIRS := $(shell for i in ./*/makefile ; do dirname $${i} ; done )

$(SUBDIRS)::
	$(MAKE) -C $@ $(MAKECMDGOALS)

all clean : $(SUBDIRS)