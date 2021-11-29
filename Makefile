CPP1 = $(CURDIR)/cpp/extract_discordant
CPP2 = $(CURDIR)/cpp/extract_unmapped
CPP3 = $(CURDIR)/cpp/convert_rep_to_2bit_k11
CPP4 = $(CURDIR)/cpp/remove_multimapping_reads_from_fa
CPP5 = $(CURDIR)/cpp/save_redundant_kmers
LHTS = $(CURDIR)/external/htslib

all: $(CPP1).so $(CPP2).so $(CPP3).so $(CPP4).so $(CPP5).so

$(CPP1).so:
	g++ -shared -fPIC -I $(LHTS) -L $(LHTS) -pthread -O2 -o $@ $(CPP1).cpp -lhts -Wl,-rpath=$(LHTS)

$(CPP2).so:
	g++ -shared -fPIC -I $(LHTS) -L $(LHTS) -pthread -O2 -o $@ $(CPP2).cpp -lhts -Wl,-rpath=$(LHTS)

$(CPP3).so:
	g++ -shared -fPIC -O2 -o $@ $(CPP3).cpp

$(CPP4).so:
	g++ -shared -fPIC -O2 -o $@ $(CPP4).cpp

$(CPP5).so:
	g++ -shared -fPIC -O2 -o $@ $(CPP5).cpp

clean:
	rm -f $(CPP1).so $(CPP2).so $(CPP3).so $(CPP4).so $(CPP5).so

re: clean all
