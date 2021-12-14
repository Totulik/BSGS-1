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

#include "BSGS.h"
#include <fstream>
#include "SECPK1/IntGroup.h"
#include "Timer.h"
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#ifndef WIN64
#include <pthread.h>
#endif

extern "C" {
#include "filter/block.h"
}

using namespace std;

// Baby steps
Point Gn[CPU_GRP_SIZE / 2];
Point _2Gn;

// Giant steps
Point GSn[CPU_GRP_SIZE / 2];
Point _2GSn;

libfilter_block bloomFilter;
libfilter_block keysToSearchBloomFilter;

// ----------------------------------------------------------------------------

BSGS::BSGS(Secp256K1 *secp,bool randomFlag,double maxStep) {

  this->secp = secp;
  this->randomFlag = randomFlag;
  this->maxStep = maxStep;

  // Compute Generator table G[n] = (n+1)*G (Baby steps group adding table)
  Point g = secp->G;
  Gn[0] = g;
  g = secp->DoubleDirect(g);
  Gn[1] = g;
  for(int i = 2; i < CPU_GRP_SIZE / 2; i++) {
    g = secp->AddDirect(g,secp->G);
    Gn[i] = g;
  }
  // _2Gn = CPU_GRP_SIZE*G
  _2Gn = secp->DoubleDirect(Gn[CPU_GRP_SIZE / 2 - 1]);

  // Init mutex
#ifdef WIN64
  ghMutex = CreateMutex(NULL,FALSE,NULL);
#else
  ghMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

}

// ----------------------------------------------------------------------------

bool BSGS::ParseConfigFile(std::string fileName) {

  // Check file
  FILE *fp = fopen(fileName.c_str(),"rb");
  if(fp == NULL) {
    ::printf("Error: Cannot open %s %s\n",fileName.c_str(),strerror(errno));
    return false;
  }
  fclose(fp);

  // Get lines
  vector<string> lines;
  int nbLine = 0;
  string line;
  ifstream inFile(fileName);
  while(getline(inFile,line)) {

    // Remove ending \r\n
    int l = (int)line.length() - 1;
    while(l >= 0 && isspace(line.at(l))) {
      line.pop_back();
      l--;
    }

    if(line.length() > 0) {
      lines.push_back(line);
      nbLine++;
    }

  }

  if(lines.size()<4) {
    ::printf("Error: %s not enough arguments\n",fileName.c_str());
    return false;
  }

  bsSize = ::stoll(lines[0],NULL,16);
  rangeStart.SetBase16((char *)lines[1].c_str());
  rangeEnd.SetBase16((char *)lines[2].c_str());
  for(int i=3;i<(int)lines.size();i++) {
    
    Point p;
    bool isCompressed;
    if( !secp->ParsePublicKeyHex(lines[i],p,isCompressed) ) {
      ::printf("%s, error line %d: %s\n",fileName.c_str(),i,lines[i].c_str());
    }
    keysToSearch.push_back(p);

  }

  nbKeysToSearch = keysToSearch.size();

#ifdef WIN64
  ::printf("BabyStep:0x%016I64X (2^%.2f)\n",bsSize,log2((double)bsSize));
#else
  ::printf("BabyStep:0x%" PRIx64 " (2^%.2f)\n",bsSize,log2((double)bsSize));
#endif

  ::printf("Start:%s\n",rangeStart.GetBase16().c_str());
  ::printf("Stop :%s\n",rangeEnd.GetBase16().c_str());
  ::printf("Keys :%d\n",nbKeysToSearch);

  return true;

}


// ----------------------------------------------------------------------------

void BSGS::FillBabySteps(TH_PARAM *ph) {

  // Global init
  int thId = ph->threadId;
  counters[thId] = 0;
  uint64_t nbStep = kPerThread / CPU_GRP_SIZE;

  // CPU Thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE / 2 + 1);

  // Group Init
  Point startP;

  Int dx[CPU_GRP_SIZE / 2 + 1];
  Point pts[CPU_GRP_SIZE];

  Int dy;
  Int dyn;
  Int _s;
  Int _p;
  Point pp;
  Point pn;
  grp->Set(dx);

  Int km(&ph->startKey);
  km.Add((uint64_t)CPU_GRP_SIZE / 2 - 1);
  startP = secp->ComputePublicKey(&km);
  startP = secp->AddDirect(startP, keyToSearch);

  ph->hasStarted = true;

#ifdef WIN64
    ::printf("BabyStep Thread %d: 0x%016I64X -> 0x%016I64X\n",ph->threadId,ph->startKey.bits64[0],ph->startKey.bits64[0]+kPerThread-1);
#else
    ::printf("BabyStep Thread %d: 0x%" PRIx64 " -> 0x%" PRIx64 "\n",ph->threadId,ph->startKey.bits64[0],ph->startKey.bits64[0]+kPerThread-1);
#endif

  // Baby Step Hashtable contains G,2G,3G,.....,(bsSize).G

  for(uint64_t s=0;s<nbStep;s++) {

    // Fill group
    int i;
    int hLength = (CPU_GRP_SIZE / 2 - 1);

    for(i = 0; i < hLength; i++) {
      dx[i].ModSub(&Gn[i].x,&startP.x);
    }
    dx[i].ModSub(&Gn[i].x,&startP.x);  // For the first point
    dx[i + 1].ModSub(&_2Gn.x,&startP.x); // For the next center point

    // Grouped ModInv
    grp->ModInv();

    // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
    // We compute key in the positive and negative way from the center of the group

    // center point
    pts[CPU_GRP_SIZE / 2] = startP;

    for(i = 0; i<hLength; i++) {

      pp = startP;
      pn = startP;

      // P = startP + i*G
      dy.ModSub(&Gn[i].y,&pp.y);

      _s.ModMulK1(&dy,&dx[i]);        // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pp.x.ModNeg();
      pp.x.ModAdd(&_p);
      pp.x.ModSub(&Gn[i].x);           // rx = pow2(s) - p1.x - p2.x;

#if 0
      pp.y.ModSub(&Gn[i].x,&pp.x);
      pp.y.ModMulK1(&_s);
      pp.y.ModSub(&Gn[i].y);           // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
      dyn.Set(&Gn[i].y);
      dyn.ModNeg();
      dyn.ModSub(&pn.y);

      _s.ModMulK1(&dyn,&dx[i]);      // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pn.x.ModNeg();
      pn.x.ModAdd(&_p);
      pn.x.ModSub(&Gn[i].x);          // rx = pow2(s) - p1.x - p2.x;

#if 0
      pn.y.ModSub(&Gn[i].x,&pn.x);
      pn.y.ModMulK1(&_s);
      pn.y.ModAdd(&Gn[i].y);          // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
      pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;

    }

    // First point (startP - (GRP_SZIE/2)*G)
    pn = startP;
    dyn.Set(&Gn[i].y);
    dyn.ModNeg();
    dyn.ModSub(&pn.y);

    _s.ModMulK1(&dyn,&dx[i]);
    _p.ModSquareK1(&_s);

    pn.x.ModNeg();
    pn.x.ModAdd(&_p);
    pn.x.ModSub(&Gn[i].x);

#if 0
    pn.y.ModSub(&Gn[i].x,&pn.x);
    pn.y.ModMulK1(&_s);
    pn.y.ModAdd(&Gn[i].y);
#endif

    pts[0] = pn;

    // Add to table
    LOCK(ghMutex);
    for(uint64_t i=0;i<CPU_GRP_SIZE;i++) {
      libfilter_block_add_hash(pts[i].x.bits64[0], &bloomFilter);
      // hashTable.Add(&pts[i].x, (uint64_t)( ph->startKey.bits64[0] + s*CPU_GRP_SIZE + i) );
      // ::printf("x[%" PRIu64 "]:%s\n", i, pts[i].x.GetBase16().c_str());
    }
    UNLOCK(ghMutex);
    counters[thId] += CPU_GRP_SIZE;

    // Next start point (startP + GRP_SIZE*G)
    pp = startP;
    dy.ModSub(&_2Gn.y,&pp.y);

    _s.ModMulK1(&dy,&dx[i + 1]);
    _p.ModSquareK1(&_s);

    pp.x.ModNeg();
    pp.x.ModAdd(&_p);
    pp.x.ModSub(&_2Gn.x);

    pp.y.ModSub(&_2Gn.x,&pp.x);
    pp.y.ModMulK1(&_s);
    pp.y.ModSub(&_2Gn.y);
    startP = pp;

  }

  delete grp;
  ph->isRunning = false;

}

// ----------------------------------------------------------------------------

int64_t BSGS::BinarySearch(Point &p, Int pk) {
  int64_t max = bsSize - 1;
  int64_t min = 0;
  int64_t middle;
  int64_t ret = 0;

  Point pointMiddle;
  Point pointCurrent;
  Point pkP;
  Point tagetP;

  pkP = secp->ComputePublicKey(&pk);

  while (min <= max) {
    middle = min + (max - min) / 2;
    // ::printf("min:%ld,max:%ld,mid:%ld\n", min, max, middle);
    Int current(middle);
    pointMiddle = secp->ComputePublicKey(&current);
    pointMiddle.y.ModNeg();
    tagetP = secp->AddDirect(pkP, pointMiddle);
    // ::printf("target:%s\n", tagetP.x.GetBase16().c_str());
    // ::printf("current:%" PRIx64 " middle:%s\n", middle, pointMiddle.x.GetBase16().c_str());
    if ( libfilter_block_find_hash(tagetP.x.bits64[0], &keysToSearchBloomFilter) ) {
      ret = middle;
      break;
    } else {
      pointCurrent = secp->AddDirect(p, pointMiddle);

      if (libfilter_block_find_hash(pointCurrent.x.bits64[0], &bloomFilter)) {
        min = middle + 1;
      } else {
        max = middle - 1;
      }
    }
  }

  // ::printf("ret2:%ld\n", ret);
  return ret;
}

// ----------------------------------------------------------------------------

void BSGS::SolveKey(TH_PARAM *ph) {

  // Global init
  int thId = ph->threadId;
  counters[thId] = 0;
  vector<uint64_t> off;

  // CPU Thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE / 2 + 1);

  // Group Init
  Point startP;

  Int dx[CPU_GRP_SIZE / 2 + 1];
  Point pts[CPU_GRP_SIZE];

  Int dy;
  Int dyn;
  Int _s;
  Int _p;
  Point pp;
  Point pn;
  grp->Set(dx);

  // Substart startRange to the point to solve
  Int km(&ph->startKey);
  km.Neg();
  km.Add(&secp->order);
  km.Sub((uint64_t)(CPU_GRP_SIZE/2)*bsSize);
  startP = secp->ComputePublicKey(&km);
  // startP = secp->AddDirect(keyToSearch,startP);

  ph->hasStarted = true;

  // if(keyIdx==0)
  ::printf("GiantStep Thread %d: %s\n",ph->threadId,ph->startKey.GetBase16().c_str());

  // Substart ((s*CPU_GRP_SIZE+i)*bsSize).G to the point to solve and look for a match into the hashtable

  Int s((uint64_t)0);
  while( s.IsLower(&ph->nbStep) && !endOfSearch ) {

    // Fill group
    int i;
    int hLength = (CPU_GRP_SIZE / 2 - 1);

    for(i = 0; i < hLength; i++) {
      dx[i].ModSub(&GSn[i].x,&startP.x);
    }
    dx[i].ModSub(&GSn[i].x,&startP.x);  // For the first point
    dx[i + 1].ModSub(&_2GSn.x,&startP.x); // For the next center point

    // Grouped ModInv
    grp->ModInv();

    // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
    // We compute key in the positive and negative way from the center of the group

    // center point
    pts[CPU_GRP_SIZE / 2] = startP;

    for(i = 0; i<hLength; i++) {

      pp = startP;
      pn = startP;

      // P = startP + i*G
      dy.ModSub(&GSn[i].y,&pp.y);

      _s.ModMulK1(&dy,&dx[i]);        // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pp.x.ModNeg();
      pp.x.ModAdd(&_p);
      pp.x.ModSub(&GSn[i].x);           // rx = pow2(s) - p1.x - p2.x;

#if 0
      pp.y.ModSub(&GSn[i].x,&pp.x);
      pp.y.ModMulK1(&_s);
      pp.y.ModSub(&GSn[i].y);           // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
      dyn.Set(&GSn[i].y);
      dyn.ModNeg();
      dyn.ModSub(&pn.y);

      _s.ModMulK1(&dyn,&dx[i]);       // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pn.x.ModNeg();
      pn.x.ModAdd(&_p);
      pn.x.ModSub(&GSn[i].x);          // rx = pow2(s) - p1.x - p2.x;

#if 0
      pn.y.ModSub(&GSn[i].x,&pn.x);
      pn.y.ModMulK1(&_s);
      pn.y.ModAdd(&GSn[i].y);          // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
      pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;

    }

    // First point (startP - (GRP_SZIE/2)*G)
    pn = startP;
    dyn.Set(&GSn[i].y);
    dyn.ModNeg();
    dyn.ModSub(&pn.y);

    _s.ModMulK1(&dyn,&dx[i]);
    _p.ModSquareK1(&_s);

    pn.x.ModNeg();
    pn.x.ModAdd(&_p);
    pn.x.ModSub(&GSn[i].x);

#if 0
    pn.y.ModSub(&GSn[i].x,&pn.x);
    pn.y.ModMulK1(&_s);
    pn.y.ModAdd(&GSn[i].y);
#endif

    pts[0] = pn;

    // Check key
    for (uint64_t i = 0; i < CPU_GRP_SIZE; i++) {

      if (libfilter_block_find_hash(pts[i].x.bits64[0], &bloomFilter)) {
        // ::printf("x[S:%d][I:%" PRIu64 "]:%s\n", s.GetInt32(), i, pts[i].x.GetBase16().c_str());
        // ::printf("y[S:%d][I:%" PRIu64 "]:%s\n", s.GetInt32(), i, pts[i].y.GetBase16().c_str());

        int64_t offset;
        Int pk(bsSize);
        Int bigS(&s);
        bigS.Mult((uint64_t)(CPU_GRP_SIZE));
        bigS.Add(i);
        pk.Mult(&bigS);
        pk.Add(&ph->startKey);
        pts[i].y = secp->GetY(pts[i].x, false);
        offset = BinarySearch(pts[i], pk);
        if (offset <= 0) {
          pts[i].y.ModNeg();
          offset = BinarySearch(pts[i], pk);
        }
        if (offset > 0) {
          // ::printf("pk:0x%s\n", pk.GetBase16().c_str());
          Int bigO(offset);
          pk.Sub(&bigO);
          // pk.Add(&ph->startKey);
          // Check
          Point p = secp->ComputePublicKey(&pk);
          // ::printf("pk2:0x%s\n", pk.GetBase16().c_str());
          for (keyIdx = 0; keyIdx < keysToSearch.size(); keyIdx++) {
            keyToSearch = keysToSearch[keyIdx];
            if (p.x.IsEqual(&keyToSearch.x)) {
              // Key solved
              FILE *f = fopen("KEYFOUNDKEYFOUND.txt", "a+");
              if (f != NULL) {
                ::fprintf(f, "\nKey#%2d Pub: 0x%s \n", keyIdx, secp->GetPublicKeyHex(true, p).c_str());
                ::fprintf(f, "     Priv: 0x%s \n", pk.GetBase16().c_str());
                fclose(f);
              }
              ::printf("\nKey#%2d Pub: 0x%s \n", keyIdx, secp->GetPublicKeyHex(true, p).c_str());
              ::printf("      Priv: 0x%s \n", pk.GetBase16().c_str());
              nbFindedKey++;
              keysToSearch.erase(keysToSearch.begin() + keyIdx);
              endOfSearch = true;
            }
          }
        }
      }
    }
    counters[thId] += CPU_GRP_SIZE;

    // Next start point (startP += (bsSize*GRP_SIZE).G)
    pp = startP;
    dy.ModSub(&_2GSn.y,&pp.y);

    _s.ModMulK1(&dy,&dx[i + 1]);
    _p.ModSquareK1(&_s);

    pp.x.ModNeg();
    pp.x.ModAdd(&_p);
    pp.x.ModSub(&_2GSn.x);

    pp.y.ModSub(&_2GSn.x,&pp.x);
    pp.y.ModMulK1(&_s);
    pp.y.ModSub(&_2GSn.y);
    startP = pp;

    s.AddOne();
  }

  delete grp;
  ph->isRunning = false;

}

// ----------------------------------------------------------------------------

#ifdef WIN64
DWORD WINAPI _FillBS(LPVOID lpParam) {
#else
void *_FillBS(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->FillBabySteps(p);
  return 0;
}

// ----------------------------------------------------------------------------

#ifdef WIN64
DWORD WINAPI _SolveKey(LPVOID lpParam) {
#else
void *_SolveKey(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->SolveKey(p);
  return 0;
}

// ----------------------------------------------------------------------------

void BSGS::Run(int nbThread) {

  double t0 = Timer::get_tick();

  nbCPUThread = nbThread;
  endOfSearch = false;

  TH_PARAM *params = (TH_PARAM *)malloc(nbCPUThread * sizeof(TH_PARAM));
  THREAD_HANDLE *thHandles = (THREAD_HANDLE *)malloc(nbCPUThread * sizeof(THREAD_HANDLE));

  memset(params, 0, nbCPUThread * sizeof(TH_PARAM));
  memset(counters, 0, sizeof(counters));
  ::printf("Number of CPU thread: %d\n", nbCPUThread);

  // Check input parameters
  if( bsSize > (1ULL<<36) ) {
    ::printf("Error, BS size cannot exceed 2^36 !\n");
    exit(0);
  }
  uint64_t k = 1;
  kPerThread = bsSize/nbThread;
  if( bsSize%nbThread != 0 ) {
    int gSize = CPU_GRP_SIZE;
    ::printf("Warning, BS size is not a multiple of nbThread\n");
  }
  if(kPerThread % CPU_GRP_SIZE != 0) {
    int gSize = CPU_GRP_SIZE;
    ::printf("Warning, BSSize/nbThread is not a multiple of %d\n",gSize);
  }

  /* init bloom filter */
  unsigned need_bytes = libfilter_block_bytes_needed(bsSize*keysToSearch.size(), 0.0001);
  libfilter_block_init(need_bytes, &bloomFilter);
  double bloomFilterSize = (libfilter_block_size_in_bytes(&bloomFilter) / (1024.0 * 1024.0));
  ::printf("Bloom filter size:%.1f MB\n", bloomFilterSize);
  unsigned keysToSearchNeedBytes = libfilter_block_bytes_needed(keysToSearch.size(), 0.0001);
  libfilter_block_init(keysToSearchNeedBytes, &keysToSearchBloomFilter);


  /* Fill baby step */
  for (keyIdx = 0; keyIdx < keysToSearch.size(); keyIdx++) {
    keyToSearch = keysToSearch[keyIdx];
    libfilter_block_add_hash(keyToSearch.x.bits64[0], &keysToSearchBloomFilter);
    // Launch Baby Step threads
    for (int i = 0; i < nbCPUThread; i++) {
      params[i].threadId = i;
      params[i].isRunning = true;
      params[i].startKey.bits64[0] = k;
      thHandles[i] = LaunchThread(_FillBS, params + i);
      k += kPerThread;
    }

    // Wait for end of baby step calculation
    Process(params, "MKey/s");
    JoinThreads(thHandles, nbCPUThread);
    FreeHandles(thHandles, nbCPUThread);
    k = 1;
  }

  // Compute Giant steps adding table GSn[n] = -(n+1)*BS
  Int bs(bsSize);
  bs.Neg();
  bs.Add(&secp->order);

  Point bsP = secp->ComputePublicKey(&bs);
  Point g = bsP;
  GSn[0] = g;
  g = secp->DoubleDirect(g);
  GSn[1] = g;
  for(int i = 2; i < CPU_GRP_SIZE / 2; i++) {
    g = secp->AddDirect(g,bsP);
    GSn[i] = g;
  }
  // _2GSn = -CPU_GRP_SIZE*BS
  _2GSn = secp->DoubleDirect(GSn[CPU_GRP_SIZE / 2 - 1]);

  nbFindedKey == 0;
  while (nbKeysToSearch != nbFindedKey) {
    // Compute range per thread
    Int bs(bsSize);
    Int nbTh;
    Int r;
    Int sk;
    nbTh.SetInt32(nbThread);
    Int rgPerTh(&rangeEnd);
    if (randomFlag) {
      Int randStart;
      randStart.Rand(&rangeStart, &rangeEnd);
      rgPerTh.Sub(&randStart);
      sk.Set(&randStart);
    }
    else {
      rgPerTh.Sub(&rangeStart);
      sk.Set(&rangeStart);
    }
    rgPerTh.AddOne();
    rgPerTh.Div(&nbTh, &r);
    if (!r.IsZero()) {
      ::printf("Warning, range is not a multiple of nbThread\n");
    }
    Int stepPerThred(&rgPerTh);
    Int grpSize;
    grpSize.SetInt32(CPU_GRP_SIZE);
    stepPerThred.Div(&bs, &r);
    if (!r.IsZero())
      stepPerThred.AddOne();
    stepPerThred.Div(&grpSize, &r);
    if (!r.IsZero()) {
      ::printf("Warning, range is not a multiple of nbThread*%d\n", CPU_GRP_SIZE);
    }

    // Lanch Giant Step threads
    endOfSearch = false;
    // Int sk(&rangeStart);
    for(int i = 0; i < nbCPUThread; i++) {
      params[i].threadId = i;
      params[i].isRunning = true;
      params[i].startKey.Set(&sk);
      params[i].nbStep.Set(&stepPerThred);
      thHandles[i] = LaunchThread(_SolveKey,params + i);
      sk.Add(&rgPerTh);
    }

    // Wait for end
    Process(params,"MKey/s");
    JoinThreads(thHandles,nbCPUThread);
    FreeHandles(thHandles,nbCPUThread);

  }

  double t1 = Timer::get_tick();

  ::printf("\nDone: Total time %s \n" , GetTimeStr(t1-t0).c_str());

}


