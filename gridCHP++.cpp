/* Copyright 2023 Daniel Grier and Luke Schaeffer

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. */

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <chrono>
#include <unistd.h>
#include "clifford.hpp"
#include "clustergen.hpp"

using namespace std;


void ClusterProg::readOutcomes(const string & fn) {
   ifstream file(fn);
   string line;
   if (file.is_open()) {
      while (getline(file,line)) {
         for(char c : line) {
            if(c == '0') random_bits.push_back(ZERO);
            if(c == '1') random_bits.push_back(ONE);
         }
      }
      file.close();
   } else {
      cerr << "Unable to open outcome file \"" << fn << "\"" << endl;
      exit(EXIT_FAILURE);
   }

   // reverse so that the first bit to read is at the end;
   reverse(random_bits.begin(), random_bits.end());
}

void ClusterProg::readBases(const string & fn) {
   ifstream file(fn);
   string line;
   if (file.is_open()) {
      while (getline(file,line)) {
         for(char c : line) {
            if (c >= '0' && c <= '2') {
                bases.push_back((Basis)(c-'0'));
            }
         }
      }
      file.close();
   } else {
      cerr << "Unable to open basis file \"" << fn << "\"" << endl;
      exit(EXIT_FAILURE);
   }
}

int main(int argc, char **argv) {

   bool printCHP = false;
   srand(time(0));

   if (argc==1) {
      cout << "usage:" << endl;
      cout << "    " << argv[0] << " [options] <grid length>" << endl;
      cout << "where the options are: " << endl;
      cout << "    -b <filename>   Set the measurement bases" << endl;
      cout << "    -r <filename>   Set the measurement outcomes" << endl;
      cout << "    -c              Print CHP circuit" << endl;
      cout << "    -l              Print *only* the running time" << endl;
      cout << "    -t              Display the running time" << endl;
      cout << "    -s              Hide (silence) measurement outcomes" << endl;
      cout << "    -h              Use row-by-row holographic algorithm " << endl;
      cout << "    -n              Use naive brute force algorithm " << endl;
      cout << "    -o              Use naive brute force algorithm with random ordering of measurements" << endl;
      cout << "    -v              Brute force circuit for verification " << endl;
      cout << "    -d              Simulation without destabilizers " << endl;
      exit(0);
   }

   ClusterProg cp;
   cp.setFlag(SILENT);

   int opt;
   opterr = 0; // don't complain about invalid options, we'll handle them ourselves
   string outcomeFile;
   string basesFile;
   while ((opt = getopt(argc, argv, "tslr:b:chnovd")) != -1) {
      switch (opt) {
         case 't':
            cp.setFlag(DISPTIME);
            break;
         case 's':
            cp.clearFlag(PRINT_OUTCOMES);
            cp.setFlag(SILENT);
            break;
         case 'l':
            cp.clearFlag(PRINT_OUTCOMES);
            cp.setFlag(LACONIC);
            cp.setFlag(SILENT);
            break;
         case 'r':
            cp.setFlag(FIXEDMEASUREMENTS);
            outcomeFile = string(optarg);
            break;
         case 'b':
            cp.setFlag(FIXEDBASES);
            basesFile = string(optarg);
            break;
         case 'c':
            printCHP = true;
            cp.clearFlag(SPACE_OPT);
            break;
         case 'h':
            cp.setFlag(ROW_BY_ROW);
            break;
         case 'n':
            cp.setFlag(NAIVE);
            break;
         case 'o':
            cp.setFlag(RANDOM_ORDER);
            break;
         case 'v':
            cp.setFlag(VERIFY);
            break;
         case 'd':
            cp.setFlag(NO_DESTABILIZERS);
            break;
         default:
            cout << "Unexpected option: " << opt << endl;
            exit(EXIT_FAILURE);
      }
   }
   if (optind >= argc) { // we still need grid size
      cout << "Missing grid size" << endl;
      exit(EXIT_FAILURE);
   }

   try {
      int n = stoi(string(argv[optind]));
      if (n <= 1) {
         cout << "Need grid length to be positive; got " << n << endl;
         exit(EXIT_FAILURE);
      }
      cp.len = n;
   } catch(std::invalid_argument& e){
      cout << "Could not parse grid length: " << string(argv[optind]) << endl;
      exit(EXIT_FAILURE);
   } catch(std::out_of_range& e){
      cout << "Grid length too large: " << string(argv[optind]) << endl;
      exit(EXIT_FAILURE);
   }

   // read fixed basis file for -b option
   if (cp.isFlagSet(FIXEDBASES)) {
      cp.readBases(basesFile);
   }
   // read fixed measurement outcome file for -r option
   if (cp.isFlagSet(FIXEDMEASUREMENTS)) {
      cp.readOutcomes(outcomeFile);
   }

   auto before = chrono::high_resolution_clock::now();

   // construct the cluster
   cp.generateCluster();


   if (printCHP) { // if -c option then we're just printing the CHP representation
      cp.printCHP();
      exit(0);
   }

   cp.run();

   if (cp.isFlagSet(PRINT_OUTCOMES)) {
      cp.printOutcomes();
   }

   auto after = chrono::high_resolution_clock::now();
   chrono::duration<double> elapsed = after - before;
   if (cp.isFlagSet(DISPTIME)) {
      cout << "Total time: " << elapsed.count() << "seconds" << endl;
   }
   if (cp.isFlagSet(LACONIC)) {
      cout << elapsed.count();
   }

   return 0;
}
