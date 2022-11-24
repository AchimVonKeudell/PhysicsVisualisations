unit piperdat;

interface

var
  a0,
  amu,
  el,
  m1,m2,z1,z2,
  lambda,
  m,
  energy,ecutoff,
  q                        : single;

  dth,kel,
  n,stotal,atf             : single;
  NB                       : integer;

  edisplacement            : integer;

  ptNB                     : integer;
  c                        : array[1..100] of integer;
  x,y,z                    : array[1..100,0..1000] of single;
  e,vx,vy,vz               : array[1..100] of single;

implementation


begin
  a0:=0.0529*1e-9;
  el:=1.6e-19;
  amu:=6e-27;
  lambda:=1.309;
  m:=0.333;
  q:=0.667;

  m1:=12;
  m2:=28;

  z1:=6;
  z2:=14;
  energy:=10000;

  ecutoff:=1;
  NB:=1000;

  edisplacement:=3;
end.
