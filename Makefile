SOURCES_NO_MAIN= partitions.cpp tree.cpp labeled_tree.cpp weightbasis.cpp niceliegroup.cpp liegroupsfromdiagram.cpp gauss.cpp log.cpp niceeinsteinliegroup.cpp ricci.cpp filter.cpp permutations.cpp weightmatrix.cpp implicitmetric.cpp antidiagonal.cpp adinvariantobstruction.cpp parsetree.cpp automorphisms.cpp

SOURCES=nice.cpp $(SOURCES_NO_MAIN)

INCLUDES=arrow.h labeled_tree.h partitions.h liegroupsfromdiagram.h  permutations.h diagramprocessor.h linearinequalities.h ricci.h double_arrows_tree.h linearsolve.h taskrunner.h filter.h log.h tree.h gauss.h niceeinsteinliegroup.h weightbasis.h horizontal.h niceliegroup.h weightmatrix.h  xginac.h tree.hpp matrixbuilder.h options.h implicitmetric.h antidiagonal.h nicediagramsinpartition.h adinvariantobstruction.h includes.h diagramanalyzer.h parsetree.h automorphisms.h

DIST=Makefile COPYING README

CXXFLAGS=-g -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes -std=c++17 

LIBS=-lginac -lwedge -lcln -lgmp -lpthread -lboost_filesystem -lboost_system -lboost_program_options

INCLUDEDIR=-I/usr/local/include/wedge-0.3

.PHONY: debug
debug: $(SOURCES) $(INCLUDES) clean_log
	g++ $(INCLUDEDIR) $(SOURCES) $(LIBS) -O0 $(CXXFLAGS) -o debug
 
.PHONY: release
release: $(SOURCES) $(INCLUDES) clean_log
	g++ $(INCLUDEDIR) $(SOURCES) $(LIBS) -O3 -DNDEBUG  $(CXXFLAGS) -o release

.PHONY: dot2ps
dot2ps:
	for dir in $$(ls output* -d); do \
		cd $$dir ; \
		rm -f all.dot ; \
		rm -f *.ps ; \
		for i in $$( ls part*.dot); do \
			cat $$i >> all.dot ; \
		done ; \
		for i in $$( ls graph*.dot); do \
			dot -Tps $$i >$$i.ps ; \
		done ; \
		if [ -e all.dot ] ;\
		then \
			dot -Tps2 all.dot >all.ps ; \
		fi ; \
		cd .. ; \
	done

.PHONY: plugin
plugin: $(SOURCES_NO_MAIN) $(INCLUDES) $(source) clean_log
	g++ $(INCLUDEDIR) $(SOURCES_NO_MAIN) $(source) $(LIBS) -O3 $(CXXFLAGS) -o `basename $(source) .cpp`

.PHONY: debugplugin
debugplugin: $(SOURCES_NO_MAIN) $(INCLUDES) $(source) clean_log
	g++ $(INCLUDEDIR) $(SOURCES_NO_MAIN) $(source) $(LIBS) -O0 $(CXXFLAGS) -o `basename $(source) .cpp`

.PHONY: dist
dist: $(SOURCES) $(INCLUDES) $(DIST)
	tar -cz $(SOURCES) $(INCLUDES) $(DIST) -f DEMONbLAST.tar.gz
	
.PHONY: clean
clean: clean_log
	for dir in $$(ls output* -d); do \
		cd $$dir ; \
		rm -f *.dot *.ps ; \
		cd .. ; \
	done

.PHONY: clean_log
clean_log:
	rm log/* -f

.PHONY: test
test: $(SOURCES) $(INCLUDES)
	./test.sh
