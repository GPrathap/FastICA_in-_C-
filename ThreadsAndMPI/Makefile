# Include system specific values
ifndef SYSTEM

$(error SYSTEM must be specified)
endif

# tools to build
TOOLS = hFastICA

API = API-1.1.alpha

BIN=$(API)/bin
DOC=$(API)/doc

API_TOOLS = $(addprefix $(BIN)/,$(TOOLS))

all: APIClean objects tool $(API) $(BIN) $(DOC) $(API_TOOLS) docs eeglabs

objects:
	cd src; make all

tool:
	cd tools; make all

APIClean:
	rm -rf $(API)

$(API):
	mkdir $(API)

$(BIN):
	mkdir $(BIN)

$(BIN)/%: bin/%
	cp $< $@

clean:
	cd src; make clean; cd ../tools; make clean; cd ..; rm -rf $(API); make clean
