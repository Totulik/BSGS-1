# Baby Step Giant Step for SECPK1

A simple Baby Step Giant Step program for SecpK1.
This BSGS version is block bloom filter base https://github.com/jbapple/libfilter.

in wsl2 test passed.
55.txt
```
4000000
60000000000000
80000000000000
0385a30d8413af4f8f9e6312400f2d194fe14f02e719b24c3f83bf1fd233a8f963
```
block bloom filter mem use 210.7 MB
```
BSGS v1.2 Block Bloom Filter
BabyStep:0x4000000 (2^26.00)
Start:60000000000000
Stop :80000000000000
Keys :1
Number of CPU thread: 4
Bloom filter size:210.7 MB
BabyStep Thread 0: 0x1 -> 0x1000000
BabyStep Thread 1: 0x1000001 -> 0x2000000
BabyStep Thread 3: 0x3000001 -> 0x4000000
BabyStep Thread 2: 0x2000001 -> 0x3000000
[9.47 MKey/s][Count 2^26.00][06s]
Warning, range is not a multiple of nbThread
GiantStep Thread 0: 60000000000000
GiantStep Thread 1: 68000000000000
GiantStep Thread 3: 78000000000000
GiantStep Thread 2: 70000000000000
[5.02 MKey/s][Count 2^25.26][08s]
Key# 0 Pub: 0x0385A30D8413AF4F8F9E6312400F2D194FE14F02E719B24C3F83BF1FD233A8F963
      Priv: 0x6ABE1F9B67E114

Done: Total time 15s
```

original BSGS mem use 1323.6MB:
```
BSGS v1.1
BabyStep:0x4000000 (2^26.00)
Start:60000000000000
Stop :80000000000000
Keys :1
Number of CPU thread: 4
BabyStep Thread 0: 0x1 -> 0x1000000
BabyStep Thread 1: 0x1000001 -> 0x2000000
BabyStep Thread 2: 0x2000001 -> 0x3000000
BabyStep Thread 3: 0x3000001 -> 0x4000000
[1.89 MKey/s][Cnt 2^26.00][28s][1323.6MB]
Sort Thread 0: 00000000 -> 00800000
Sort Thread 2: 01000000 -> 01800000
Sort Thread 3: 01800000 -> 02000000
Sort Thread 1: 00800000 -> 01000000
[21.61 MSort/s][Cnt 2^25.00][01s][1323.6MB]
Warning, range is not a multiple of nbThread
GiantStep Thread 0: 60000000000000
GiantStep Thread 2: 70000000000000
GiantStep Thread 1: 68000000000000
GiantStep Thread 3: 78000000000000
[4.31 MKey/s][Cnt 2^25.42][10s][1323.6MB]
Key# 0 Pub:  0x0385A30D8413AF4F8F9E6312400F2D194FE14F02E719B24C3F83BF1FD233A8F963
       Priv: 0x6ABE1F9B67E114

Done: Total time 41s
```
