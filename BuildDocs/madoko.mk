MADOKO ?= madoko
BUILD_DIR ?= out

MARKDOWN := $(TARGETS:%=%.md)
PDF := $(TARGETS:%=%.pdf)
HTML := $(TARGETS:%=%.html)

# Shortcuts.
.PHONY: default
default: $(PDF) $(HTML)

# Build PDF via LaTeX.
%.pdf: %.md $(DEPS)
	$(MADOKO) --odir=$(BUILD_DIR) --pdf $<
	cp $(BUILD_DIR)/$@ .

# Build Web page.
%.html: %.md $(DEPS)
	$(MADOKO) --odir=$(BUILD_DIR) $<
	cp $(BUILD_DIR)/$@ .

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

.PHONY: deploy
deploy: $(PDF)
	scp $< $(DEST)
	@echo http://$(DEST_PATH)/$(notdir $(PDF))
deploy-html: $(HTML)
	scp $< $(DEST)
	@echo http://$(DEST_PATH)/$(notdir $(HTML))

.PHONY: install
install:	$(PDF) $(HTML)
		cp $(PDF)  ../Documentation
		cp $(HTML) ../Documentation
		cp $(HTML) ~/Data/GitHub/JohnKendrick.github.io


# View products.

OS=$(shell uname -s)
ifeq ($(OS),Darwin)
OPEN ?= open
else
OPEN ?= xdg-open
endif

.PHONY: view view-html
view: $(PDF)
	$(OPEN) $(PDF)
view-html: $(HTML)
	$(OPEN) $(HTML)


# Auto-build based on `livereload`.

.PHONY: watch watch-pdf watch-html

LIVESERVE_ARGS := $(MARKDOWN:%=-w %) $(DEPS:%=-w %)

watch-pdf:
	liveserve $(LIVESERVE_ARGS) -x 'make pdf' -S

watch-html:
	liveserve $(LIVESERVE_ARGS) -x 'make html' $(BUILD_DIR)

watch: watch-pdf
