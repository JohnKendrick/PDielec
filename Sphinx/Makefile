# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = PDielec
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


.PHONY: html
# Do a clean build by deleting autoapi and _build
html:
	rm -rf autoapi
	rm -rf _build
	sphinx-build -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)
	rsync -r _build/html/* ../docs

.PHONY: help Makefile
# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

livehtml:
	sphinx-autobuild -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) 

