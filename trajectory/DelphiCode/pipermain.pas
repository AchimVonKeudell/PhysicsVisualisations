unit pipermain;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls,math, Plotwin, ExtCtrls, ComCtrls, ToolWin, Buttons, Menus;

type
  TFormMain = class(TForm)
    StatusBar1: TStatusBar;
    CoolBar1: TCoolBar;
    ToolBar1: TToolBar;
    Button1: TButton;
    Panel1: TPanel;
    Plot1: TPlot;
    SpeedButton1: TSpeedButton;
    Splitter1: TSplitter;
    Plot2: TPlot;
    ToolButton1: TToolButton;
    SpeedButton2: TSpeedButton;
    ToolButton2: TToolButton;
    SpeedButton3: TSpeedButton;
    MainMenu1: TMainMenu;
    File1: TMenuItem;
    N1: TMenuItem;
    procedure Button1Click(Sender: TObject);
    procedure SpeedButton1Click(Sender: TObject);
    procedure Plot1ScaleChange(Sender: TObject; Ymin, Ymax: Single;
      SpectrumID: Integer);
    procedure SpeedButton2Click(Sender: TObject);
    procedure SpeedButton3Click(Sender: TObject);
    procedure N1Click(Sender: TObject);
  private
    AbortTRUE: Boolean;
    { Private-Deklarationen }
  public
    procedure GeneratePovrayFile;
    procedure DisplayCascade;

    procedure GenerateTrajectory(pt:integer;x0,y0,z0,vx0,vy0,vz0,en:single);
  end;

var
  FormMain: TFormMain;

implementation

{$R *.DFM}

uses piperdat, piperparameter, piperamk;

function f_func(t:single):single;
begin
  result:=lambda*Power(t,0.5-m)*
          Power(
                1+
                  Power(2*lambda*Power(t,1-m),q),

                -1/q);
end;

function t_func(th,eps:single):single;
begin
  result:=sqr(eps*sin(th/2));
end;

function ds_func(atf,t:single):single;
var f : double;
begin
  f:=f_func(t);
  result:=pi*sqr(atf)*0.5*power(t,-1.5)*f;
end;

function eps_func(atf,z1,z2,m1,m2,e:single):single;
begin
  result:=atf*m2/(z1*z2*sqr(el)*(m1+m2))*e*(4*pi*8.85e-12);
end;

function atf_func(z1,z2:single):single;
begin
  result:=a0*0.88534*power(sqrt(z1)+sqrt(z2),-2/3);
end;

function thlab_func(m1,m2,th:single):single;
begin
  result:=arccos((m1+m2*cos(th))/sqrt(sqr(m1)+2*m1*m2*cos(th)+sqr(m2)));
end;

function kelec_func(z1,z2,m1,m2:single):single;
begin
  result:=3.83*power(z1,7/6)*z2/
          (Power(power(z1,2/3)+power(z2,2/3),3/2)*power(m1,1/2));
end;

procedure TFormMain.Button1Click(Sender: TObject);
var i,j : integer;
label ausgang;
begin
  { initialisierung }
  FormParameter.GetParameter;
  AbortTRUE:=FALSE;
  for i:=1 to 100 do
    begin
      Plot2.Spectrum[1,i]:=i;
      Plot2.Spectrum[2,i]:=0;
    end;
  Plot2.SpectrumNB[1]:=100;
  Plot2.SpectrumNB[2]:=100;

  kel:=kelec_func(z1,z2,m1,m2);
  Plot1.Spectrum[3,1]:=0;
  Plot1.Spectrum[4,1]:=-100;
  Plot1.Spectrum[3,2]:=0;
  Plot1.Spectrum[4,2]:=100;

  Plot1.SpectrumNB[3]:=2;
  Plot1.SpectrumNB[4]:=2;
  Plot1.PaintDataPoint(4,2);

  
  atf:=atf_func(z1,z2);
  n:=1e23*1e6;
  stotal:=power(n,-2/3);
  dth:=0.0001;

  for i:=1 to NB do
    begin
      StatusBar1.Panels[0].Text:=' Particle : '+IntToStr(i);
      Application.ProcessMessages;
      PtNB:=1;
      GenerateTrajectory(PtNB,0,0,0,0,0,1,energy);
      DisplayCascade;
      GeneratePovRayFile;
      if AbortTRUE then goto ausgang;
    end;
  ausgang:

end;

    procedure direction(var vx1,vy1,vz1: single; vx0,vy0,vz0,t,ph: single);
    var srat,sine,cx2,cy2,cz2,un  : single;
    begin
      sine:=sqrt(vy0*vy0+vz0*vz0);
      srat:=sin(t)/sine;
      cx2:=cos(t)*vx0+sin(t)*sine*cos(ph);
      cy2:=cos(t)*vy0-srat*(vy0*vx0*cos(ph)-vz0*sin(ph));
      cz2:=cos(t)*vz0-srat*(vy0*vx0*cos(ph)-vz0*sin(ph));
      un:=1/sqrt(cx2*cx2+cy2*cy2+cz2*cz2);
      vx1:=cx2*un;
      vy1:=cy2*un;
      vz1:=cz2*un+1e-12;
    end;



procedure TFormMain.GenerateTrajectory(pt:integer;x0,y0,z0,vx0,vy0,vz0,en:single);
var t,r1,eps,th,s,
    etransfer,
    sgesamt,thc,deltae,
    thmin              : double;

    datei              : textfile;

    vxn,vyn,vzn        : single;
    vxt,vyt,vzt        : single;
    thl,phi            : single;

    i                  : integer;


label ausgang;

begin
  e[pt]:=en*el;


  { Startbedingungen }
  s:=0;

  c[pt]:=0;
  x[pt,0]:=x0;
  y[pt,0]:=y0;
  z[pt,0]:=z0;
  vx[pt]:=vx0;
  vy[pt]:=vy0;
  vz[pt]:=vz0;

  repeat

    inc(c[pt]);

    {----------------------------------------------------

      nuclear stopping

     ---------------------------------------------------}

    r1:=random*0.97;

    th:=pi/2;
    s:=0;
    repeat
      th:=th-dth;
      eps:=eps_func(atf,z1,z2,m1,m2,e[pt]);
      t:=t_func(th,eps);
      if eps*eps-t>=0 then
         s:=s+ds_func(atf,t)*sqrt(t)*sqrt(eps*eps-t)*dth;
    until (s>=stotal) or (th=0);
    thmin:=th;
    sgesamt:=s;

    th:=thmin;
    s:=0;
    repeat
      th:=th+dth;
      eps:=eps_func(atf,z1,z2,m1,m2,e[pt]);
      t:=t_func(th,eps);
      if eps*eps-t>=0 then
        s:=s+ds_func(atf,t)*sqrt(t)*sqrt(eps*eps-t)*dth;
    until s/sgesamt>=r1;

    if th>pi/2 then
      begin
        th:=0;{pi/2}
      end;

    { --------------------------------
      nuclear stopping
      -------------------------------}
    deltae:=4*(m1*m2)/sqr(m1+m2)*e[pt]*sqr(sin(th/2));
    etransfer:=deltae/el;

    { --------------------------------
      electronic stopping
      -------------------------------}
    deltae:=deltae+kel*sqrt(e[pt]/el/1e3)*power(n,-1/3)*1e15*el*1e-4;

    e[pt]:=e[pt]-deltae;

    thl:=thlab_func(m1,m2,th);
    phi:=2*pi*random;

    vxn:=vx[pt];
    vyn:=vy[pt];
    vzn:=vz[pt];

    direction(vx[pt],vy[pt],vz[pt],vxn,vyn,vzn,thl,phi);

    x[pt,c[pt]]:=x[pt,c[pt]-1]+vx[pt];
    y[pt,c[pt]]:=y[pt,c[pt]-1]+vy[pt];
    z[pt,c[pt]]:=z[pt,c[pt]-1]+vz[pt];

    { --------------------------------
      cascade
      -------------------------------}
    if (etransfer>edisplacement) and (ptnb<1) then
      begin
        inc(PtNb);
        StatusBar1.Panels[1].Text:=IntToStr(PtNB);
        StatusBar1.Refresh;
        direction(vxt,vyt,vzt,vxn,vyn,vzn,thl+pi/2,phi-pi);
        GenerateTrajectory(ptnb,x[pt,c[pt]-1],y[pt,c[pt]-1],z[pt,c[pt]-1],
                                vxt,vyt,vzt,etransfer);
      end;


     if z[pt,c[pt]]<0 then goto ausgang;
    {----------------------------------------------------

      electronic stopping

     ---------------------------------------------------}
  until e[pt]<ecutoff*el;

  ausgang:

end;


procedure TFormMain.DisplayCascade;
var k,ipt,i : integer;
begin
  Plot1.Spectrum[1,1]:=0;
  Plot1.Spectrum[2,1]:=0;
  ipt:=1;
  for k:=1 to PtNB do
    begin
       for i:=0 to c[k] do
         begin
           inc(ipt);
           Plot1.Spectrum[1,ipt]:=z[k,i];
           Plot1.Spectrum[2,ipt]:=y[k,i];
           Plot1.PaintDataPoint(2,ipt);
         end;
    end;

  Plot1.SpectrumNB[1]:=ipt;
  Plot1.SpectrumNB[2]:=ipt;
  if (z[1,c[1]]>=1) and (z[1,c[1]]<100) then
    begin
      Plot2.Spectrum[2,trunc(z[1,c[1]])]:=Plot2.Spectrum[2,trunc(z[1,c[1]])]+1;
      Plot2.UpdatePlot;
    end;
end;

procedure TFormMain.Plot1ScaleChange(Sender: TObject; Ymin, Ymax: Single;
  SpectrumID: Integer);
begin
  if SpectrumID=1 then
    begin
      Plot1.SetSpectrumRange(1,Ymin,Ymax);
      Plot1.SetSpectrumRange(3,Ymin,Ymax);
    end;
  if SpectrumID=2 then
    begin
      Plot1.SetSpectrumRange(2,Ymin,Ymax);
      Plot1.SetSpectrumRange(4,Ymin,Ymax);
    end;
end;


procedure TFormMain.SpeedButton1Click(Sender: TObject);
begin
  AbortTRUE:=TRUE;
end;


procedure TFormMain.GeneratePovrayFile;
var datei     : textfile;
    k,j         : longint;
begin
  assignfile(datei,ExtractFilePAth(Application.ExeName)+'piper.pov');
  rewrite(datei);

  writeln(datei,'  #include "colors.inc"');
  writeln(datei,'  #include "textures.inc"');
  writeln(datei,'#include "shapes.inc"');

  writeln(datei,'light_source { < 10, 10, -10 >');
  writeln(datei,'color White');
  writeln(datei,'}             ');

  writeln(datei,'light_source { < 10, 10, 10 >');
  writeln(datei,'color White                     ');
  writeln(datei,'}                                  ');

  writeln(datei,'background { color White }');

  writeln(datei,'camera {');
  writeln(datei,' right      < -1.33, 0, 0 > ');
  writeln(datei,' up         < 0, 1, 0 >');
  writeln(datei,' direction  < 0, 0, 1 >');
  writeln(datei,' location   < 15, 15, 45 >');
  writeln(datei,' look_at    < 0,0, 20 >');
  writeln(datei,'}');

  write(datei,'cylinder { <0,0,0>,<0,0,100>,');
  writeln(datei,'           0.1 pigment { color Blue } }');

  for k:=1 to PtNB do
    begin
       for j:=1 to c[k] do
         begin
           write(datei,'cylinder { <');
           writeln(datei,'           ',x[k,j-1]:5:3,',',y[k,j-1]:5:3,',',z[k,j-1]:5:3,'>,');
           writeln(datei,'           <',x[k,j]:5:3,',',y[k,j]:5:3,',',z[k,j]:5:3,'>,');
           writeln(datei,'           0.1 pigment { color Red } }');
         end;
    end;

  closefile(datei);
end;

procedure TFormMain.SpeedButton2Click(Sender: TObject);
begin
  close;
end;

procedure TFormMain.SpeedButton3Click(Sender: TObject);
begin
  FormPArameter.Visible:=SpeedButton3.Down;
end;

procedure TFormMain.N1Click(Sender: TObject);
begin
  FormAbout.ShowModal;
end;

end.
