/*
 * This file is part of the BSGS distribution (https://github.com/JeanLucPons/BSGS).
 * Copyright (c) 2020 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "Timer.h"
#include "BSGS.h"
#include "SECPK1/SECP256k1.h"
#include <fstream>
#include <string>
#include <string.h>
#include <stdexcept>

#define RELEASE "1.2 Block Bloom Filter"

using namespace std;

#define CHECKARG(opt,n) if(a>=argc-1) {::printf(opt " missing argument #%d\n",n);exit(0);} else {a++;}

// ------------------------------------------------------------------------------------------

void printUsage() {

  printf("BSGS [-v] [-t nbThread] inFile\n");
  printf(" -v: Print version\n");
  printf(" -t nbThread: Secify number of thread\n");
  printf(" -rand: random startkey\n");
  printf(" -m maxStep: rand mode, number of operations before next search,default is 32 (2^32).\n");
  printf(" inFile: intput configuration file\n");
  exit(0);

}

// ------------------------------------------------------------------------------------------

int getInt(string name,char *v) {

  int r;

  try {

    r = std::stoi(string(v));

  } catch(std::invalid_argument&) {

    printf("Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

  return r;

}

// ------------------------------------------------------------------------------------------

double getDouble(string name,char *v) {

  double r;

  try {

    r = std::stod(string(v));

  } catch(std::invalid_argument&) {

    printf("Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

  return r;

}

// ------------------------------------------------------------------------------------------

// Default params
static bool randomFlag = false;
static double maxStep = 32.0;

// ------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  // Global Init
  Timer::Init();
  rseed(Timer::getSeed32());

  // Init SecpK1
  Secp256K1 *secp = new Secp256K1();
  secp->Init();

  int a = 1;
  int nbCPUThread = Timer::getCoreNumber();
  string configFile = "";

  while (a < argc) {

    if(strcmp(argv[a], "-t") == 0) {
      CHECKARG("-t",1);
      nbCPUThread = getInt("nbCPUThread",argv[a]);
      a++;
    } else if (strcmp(argv[a], "-rand") == 0) {
      randomFlag = true;
      a++;
    } else if (strcmp(argv[a], "-m") == 0) {
      CHECKARG("-m",1);
      maxStep = getDouble("maxStep",argv[a]);
      a++;
    } else if (strcmp(argv[a], "-h") == 0) {
      printUsage();
    } else if (a == argc - 1) {
      configFile = string(argv[a]);
      a++;
    } else {
      printf("Unexpected %s argument\n",argv[a]);
      exit(-1);
    }

  }

  printf("BSGS v" RELEASE "\n");

  BSGS *v = new BSGS(secp,randomFlag,maxStep);
  if( !v->ParseConfigFile(configFile) )
    exit(-1);
  v->Run(nbCPUThread);

  return 0;

}
