# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = ./source
BUILDDIR      = ./




# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
clean:
	rm -rf docs*
	rm -r doctrees*
	
%: Makefile
	rm -rf docs*
	rm -rf doctrees*

	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) 
	mv html docs
	touch ./docs/.nojekyll

	./encode.sh docs/secret

	rm -rf doctrees/secret
	rm -rf docs/_sources/secret
