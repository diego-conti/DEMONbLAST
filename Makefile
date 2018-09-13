SOURCES=nice.cpp partitions.cpp tree.cpp labeled_tree.cpp weightbasis.cpp niceliegroup.cpp liegroupsfromdiagram.cpp gauss.cpp log.cpp niceeinsteinliegroup.cpp ricci.cpp filter.cpp permutations.cpp weightmatrix.cpp

INCLUDES=arrow.h labeled_tree.h partitions.h liegroupsfromdiagram.h  permutations.h diagramprocessor.h linearinequalities.h ricci.h double_arrows_tree.h linearsolve.h taskrunner.h filter.h log.h tree.h gauss.h niceeinsteinliegroup.h weightbasis.h horizontal.h niceliegroup.h weightmatrix.h  xginac.h tree.hpp matrixbuilder.h

DIST=Makefile COPYING README

CXXFLAGS=-g -ffor-scope -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes -std=c++17

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
		for i in $$( ls *.dot); do \
			dot -Tps2 $$i >$$i.ps ; \
			cat $$i >> all.dot ; \
		done ; \
		if [ -e all.dot ] ;\
		then \
			dot -Tps2 all.dot >all.ps ; \
		fi ; \
		cd .. ; \
	done

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
