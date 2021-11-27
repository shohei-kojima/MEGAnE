CPP1 = ./cpp/extract_discordant
CPP2 = ./cpp/extract_unmapped
CPP3 = ./cpp/convert_rep_to_2bit_k11
CPP4 = ./cpp/remove_multimapping_reads_from_fa
CPP5 = ./cpp/save_redundant_kmers
LHTS = ./external/htslib


all: $(CPP1).so $(CPP2).so $(CPP3).so $(CPP4).so $(CPP5).so

$(CPP1).so:
	gcc -xc++ -lstdc++ -shared-libgcc -shared -fPIC -I $(LHTS) -L $(LHTS) -pthread -O2 -o $@ $(CPP1).cpp

$(CPP2).so:
	gcc -xc++ -lstdc++ -shared-libgcc -shared -fPIC -I $(LHTS) -L $(LHTS) -pthread -O2 -o $@ $(CPP2).cpp

$(CPP3).so:
	gcc -xc++ -lstdc++ -shared-libgcc -shared -fPIC -O2 -o $@ $(CPP3).cpp

$(CPP4).so:
	gcc -xc++ -lstdc++ -shared-libgcc -shared -fPIC -O2 -o $@ $(CPP4).cpp

$(CPP5).so:
	gcc -xc++ -lstdc++ -shared-libgcc -shared -fPIC -O2 -o $@ $(CPP5).cpp

clean:
	rm -f $(CPP1).so $(CPP2).so $(CPP3).so $(CPP4).so $(CPP5).so

re: clean all
