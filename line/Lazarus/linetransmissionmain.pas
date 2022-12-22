unit LineTransmissionMain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, TAGraph, TASeries, Forms, Controls, Graphics,
  Dialogs, ExtCtrls, Buttons, StdCtrls, LCLType, Menus;

type

  { TForm1 }

  TForm1 = class(TForm)
    Button1: TButton;
    Chart1: TChart;
    Chart1AreaSeries1: TAreaSeries;
    Chart1AreaSeries2: TAreaSeries;
    Chart1LineSeries1: TLineSeries;
    Chart1LineSeries2: TLineSeries;
    Chart1LineSeries3: TLineSeries;
    Chart1LineSeries4: TLineSeries;
    Chart2: TChart;
    Chart2LineSeries1: TLineSeries;
    Chart2LineSeries2: TLineSeries;
    CheckBox1: TCheckBox;
    ComboBoxBoundaryRight: TComboBox;
    ComboBoxBoundaryLeft: TComboBox;
    EditXsize: TEdit;
    Editdx: TEdit;
    EditBoundary2: TEdit;
    EditBoundary1: TEdit;
    EditStoreEachBMPDt: TEdit;
    EditDrawEachDt: TEdit;
    EditRefractiveIndexPlug: TEdit;
    Edittend: TEdit;
    EditTimeFaktor: TEdit;
    EditDamping: TEdit;
    EditRefractiveIndexMedium: TEdit;
    EditRefractiveIndexAmbient: TEdit;
    GroupBox1: TGroupBox;
    GroupBox2: TGroupBox;
    GroupBox3: TGroupBox;
    GroupBox4: TGroupBox;
    GroupBox5: TGroupBox;
    Label1: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    Label12: TLabel;
    Label13: TLabel;
    Label14: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    MainMenu1: TMainMenu;
    MenuItem1: TMenuItem;
    MenuItemSaveSignal: TMenuItem;
    MenuItemReadData: TMenuItem;
    OpenDialog1: TOpenDialog;
    Panel1: TPanel;
    Panel2: TPanel;
    Panel3: TPanel;
    Panel4: TPanel;
    SaveDialog1: TSaveDialog;
    Splitter1: TSplitter;
    procedure Button1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure MenuItem1Click(Sender: TObject);
    procedure MenuItemSaveSignalClick(Sender: TObject);
    procedure MenuItemReadDataClick(Sender: TObject);
    procedure Panel3Click(Sender: TObject);
  private
    { private declarations }
  public
    tend,
    dt             : single;
    DrawEachStep,
    DrawEachBMP    : integer;
    dx             : single;
    boundary1      : single;
    boundary2      : single;
    XSize          : integer;
    line  : array [1..10000,0..2] of double;
    dline : array [1..10000] of double;
    alphavector : array [1..10000] of double;
    signal         : array [1..10000] of single;
    signalreflected: array [1..10000] of single;

    TimeTraceNB      : integer;
    timetraceBCS,timetraceElectrode : array [1..1000,1..2] of single;

    DataTRUE       : boolean;
    SimulationTRUE : boolean;

    RefractiveIndexPlug,
    RefractiveIndexMedium,
    RefractiveIndexAmbient,
    gamma          : single;

    BoundaryConditionLeft,
    BoundaryConditionRight   : integer;

    procedure ReadData(Filename:string;tforwardbegin,tforwardend,treflectedbegin,treflectedend,startposition: single);
    procedure GetParameter;
    procedure RunSimulation;
    { public declarations }
  end;

var
  Form1: TForm1;

implementation

uses linetransmissionreaddata;

{$R *.lfm}

procedure TForm1.Button1Click(Sender: TObject);
begin
  if SimulationTRUE=FALSE then
    begin
      SimulationTRUE:=TRUE;
      Button1.Caption:='Stop';
      RunSimulation;
      SimulationTRUE:=FALSE;
      Button1.Caption:='Start';
    end
  else
    begin
      SimulationTRUE:=FALSE;
      Button1.Caption:='Start';
    end;
end;

procedure TForm1.FormCreate(Sender: TObject);
begin
  ComboBoxBoundaryLeft.Itemindex:=0;
  XSize:=10000;
  dx:=0.002;

  dt:=1000*1e-15;
  DrawEachBMP:=100;
  DrawEachStep:=1000;

  Editdx.Text:=FloatToStrF(dx,ffgeneral,5,3);
  EditXsize.Text:=IntToStr(Xsize);

  EditTimeFaktor.Text:=IntTostr(round(dt/1e-15));
  EditStoreEachBMPDt.Text:=IntToStr(round(DrawEachBMP*dt/1e-15));
  EditDrawEachDt.Text:=IntToStr(round(DrawEachStep*dt/1e-15));

  boundary1:=12;
  boundary2:=12.2;

  EditBoundary1.Text:=FloatToStrF(boundary1,fffixed,5,2);
  EditBoundary2.Text:=FloatToStrF(boundary2,fffixed,5,2);

  Chart1LineSeries2.Clear;
  Chart1AreaSeries1.Clear;

  Chart1.LeftAxis.Range.Min:=-2;
  Chart1.LeftAxis.Range.Max:=2;

  Chart1LineSeries2.AddXY(0,0);
  Chart1LineSeries2.AddXY(XSize*dx,0);
  Chart1AreaSeries1.AddXY(boundary1,2);
  Chart1AreaSeries1.AddXY(boundary2,2);
  Chart1AreaSeries2.AddXY(boundary2,2);
  Chart1AreaSeries2.AddXY(XSize*dx,2);

  DataTRUE:=FALSE;
  SimulationTRUE:=FALSE;

  BoundaryConditionLeft:=1;
  BoundaryConditionRight:=2;

  ComboBoxBoundaryLeft.Itemindex:=BoundaryConditionLeft;
  ComboBoxBoundaryRight.Itemindex:=BoundaryConditionRight;

end;

procedure TForm1.MenuItem1Click(Sender: TObject);
begin

end;

procedure TForm1.MenuItemSaveSignalClick(Sender: TObject);
var i     : integer;
    datei : textfile;
begin
  {-----------------------------------------------------------------------
   Save Reflected Signal
   ----------------------------------------------------------------------}
{  SaveDialog1.Title:='Save reflected Signal';
  if SaveDialog1.Execute then
    begin
      assignfile(datei,SaveDialog1.FileName);
      rewrite(datei);
      writeln(datei,'t reflected');
      for i:=1 to Xsize do writeln(datei,(10-i*dx)*RefractiveIndexAmbient/3e8/1e-9,' ',line[i,1]);
      closefile(datei);
    end;  }
  {-----------------------------------------------------------------------
   Save Transmitted Signal
   ----------------------------------------------------------------------}
{  SaveDialog1.Title:='Save transmitted Signal';
  if SaveDialog1.Execute then
    begin
      assignfile(datei,SaveDialog1.FileName);
      rewrite(datei);
      writeln(datei,'t transmitted');
      for i:=1 to Xsize do writeln(datei,(12.5-i*dx)*RefractiveIndexMedium/3e8/1e-9,' ',line[i,1]);
      closefile(datei);
    end;  }

  {-----------------------------------------------------------------------
   Save Probe Pickup Details
   ----------------------------------------------------------------------}
  SaveDialog1.Title:='Save Signals';
  if SaveDialog1.Execute then
    begin
      assignfile(datei,SaveDialog1.FileName);
      rewrite(datei);
      writeln(datei,'t BCS Electrode');
      for i:=1 to TimeTraceNB do writeln(datei,TimeTraceBCS[i,1],' ',TimeTraceBCS[i,2],' ',TimeTraceElectrode[i,2]);
      closefile(datei);
    end;


end;

procedure TForm1.MenuItemReadDataClick(Sender: TObject);
begin
  Form2.ShowModal;
end;

procedure TForm1.ReadData(Filename:string;tforwardbegin,tforwardend,treflectedbegin,treflectedend,startposition: single);
var t,t0,
    dtsignal               : single;
    max,x,y                : single;
    i,j,DtPt               : integer;
    rawsignal              : array [1..10000,1..2] of single;
    rawsignalNB            : integer;
    datei                  : textfile;

begin
  {tforwardbegin:=-2;
  tforwardend:=30;

  startposition:=9;

  treflectedbegin:=57;
  treflectedend:=90; }


  dtsignal:=1/3e8*dx*1/1e-9;
  for i:=1 to XSize do signal[i]:=0;
  DtPt:=0;

       assignfile(datei,FileName);
       reset(datei);
       readln(datei,x,y);
       t0:=x;
       DtPt:=0;
       repeat
         readln(datei,x,y);
         inc(DtPt);
         rawsignal[DtPt,1]:=x;
         rawsignal[DtPt,2]:=y;
       until eof(datei);
       rawsignalNB:=DtPt;
       closefile(datei);

       { forward signal }
       for i:=1 to XSize do
         begin
           x:=startposition-i*dx;
           t:=1/3e8*x/1e-9;
           for j:=1 to RawSignalNB do
             begin
               if (t>(rawsignal[j,1]-tforwardbegin)) and
                  (t<(rawsignal[j+1,1]-tforwardbegin)) and
                  (t< tforwardend-tforwardbegin) then
                 begin
                   signal[i]:=rawsignal[j,2];
                 end;
             end;
         end;
       max:=0;
       for i:=1 to XSize do
         begin
           if signal[i]>max then max:=signal[i];
         end;
       for i:=1 to Xsize do
           begin
             signal[i]:=signal[i]/max-0.15;
           end;
       for i:=trunc(startposition/dx) to Xsize do
           begin
             signal[i]:=0;
           end;


       { reflected signal }
       for i:=1 to XSize do
         begin
           x:=startposition-i*dx;
           t:=1/3e8*x/1e-9;
           for j:=1 to RawSignalNB do
             begin
               if (t>(rawsignal[j,1]-treflectedbegin)) and
                  (t<(rawsignal[j+1,1]-treflectedbegin)) and
                  (t< treflectedend-treflectedbegin) then
                 begin
                   signalreflected[i]:=rawsignal[j,2]*(-1);
                 end;
             end;
         end;
       for i:=1 to Xsize do signalreflected[i]:=signalreflected[i]/max-0.15;
       for i:=trunc(startposition/dx) to Xsize do
           begin
             signalreflected[i]:=0;
           end;

       Chart1LineSeries1.Clear;
       Chart1LineSeries3.Clear;
       Chart1LineSeries4.Clear;
       for i:=1 to Xsize do
         begin
           Chart1LineSeries1.AddXY(i*dx,signal[i]);
           Chart1LineSeries3.AddXY(i*dx,signalreflected[i]);
           Chart1LineSeries4.AddXY(i*dx,signal[i]);
         end;
       Chart1.LeftAxis.Range.Min:=-2;
       Chart1.LeftAxis.Range.Max:=2;

       DataTRUE:=TRUE;

end;

procedure TForm1.Panel3Click(Sender: TObject);
begin

end;

procedure TForm1.GetParameter;
begin
  dx:=StrToFloat(Editdx.Text);
  Xsize:=StrToInt(EditXsize.Text);

  dt:=StrToFloat(EditTimeFaktor.Text)*1e-15;

  DrawEachBMP:=trunc(StrToFloat(EditStoreEachBMPDt.Text)*1e-15/dt);
  DrawEachStep:=trunc(StrToFloat(EditDrawEachDt.Text)*1e-15/dt);

  BoundaryConditionLeft:=ComboBoxBoundaryLeft.Itemindex;
  BoundaryConditionRight:=ComboBoxBoundaryRight.Itemindex;

  boundary1:=StrToFloat(EditBoundary1.Text);
  boundary2:=StrToFloat(EditBoundary2.Text);

  tend:=StrToFloat(EditTend.Text)*1e-9;

  RefractiveIndexAmbient:=StrToFloat(EditRefractiveIndexAmbient.Text);
  RefractiveIndexMedium:=StrToFloat(EditRefractiveIndexMedium.Text);
  RefractiveIndexPlug:=StrToFloat(EditRefractiveIndexPlug.Text);
  gamma:=StrToFloat(EditDamping.Text);

end;


{---------------------------------------------------------------------

 Solution via Euler Scheme FTDT

 ---------------------------------------------------------------------}
procedure TForm1.RunSimulation;
var t   : single;
    c2,
    alpha,
    alpha1,
    alpha2,
    vphase     : double;
    i : integer;
    k,l,
    DtPt : integer;
    DrawDt : double;

    TimeFaktor        : single;

    gamma0            : single;

    Bitmap1           : TBitMap;

    TimeTracePt,
    BMPCounter,
    BmpPt             : integer;

    DestRect,SourceRect : TRect;
    datei               : textfile;

    BCSpt,ElectrodePt   : integer;

begin
  GetParameter;

  BmpPt:=0;
  BMPCounter:=0;
  t:=0;
  vphase:=3e8;
  gamma0:=0;
  alpha:=vphase*dt/dx;

  assignfile(datei,'timetraceboundary2.dat');
  rewrite(datei);

  Bitmap1 := TBitmap.Create;
  Bitmap1.Width:=Chart1.Width;
  Bitmap1.Height:=Chart1.Height;

  Chart1LineSeries2.Clear;
  Chart1AreaSeries1.Clear;
  Chart1AreaSeries2.Clear;

  Chart1LineSeries2.AddXY(0,0);
  Chart1LineSeries2.AddXY(XSize*dx,0);
  Chart1AreaSeries1.AddXY(boundary1,2);
  Chart1AreaSeries1.AddXY(boundary2,2);
  Chart1AreaSeries2.AddXY(boundary2,2);
  Chart1AreaSeries2.AddXY(Xsize*dx,2);


  Chart2LineSeries1.Clear;
  Chart2LineSeries2.Clear;

  {Check von Neumann Diffusion condition }
  if alpha < 1 then
    SimulationTRUE:=TRUE
  else
    begin
      Application.MessageBox('von Neumann condition violated - reduce time steps','Simulation Stopped',MB_OK);
      SimulationTRUE:=FALSE;
      Button1.Caption:='Start';
      closefile(datei);
      exit;
    end;



  for i:=1 to Xsize do
    begin
      line[i,0]:=0;
      line[i,1]:=0;
      line[i,2]:=0;
      dline[i]:=0;
    end;

  for i:=1 to Xsize do
    begin
      if i<=boundary1/dx then
        alphavector[i]:=alpha
      else
        alphavector[i]:=alpha*1/(RefractiveIndexMedium);
    end;

  { Starting Condition Amplitude }
  if DataTRUE then
    begin
      for i:=1 to Xsize do line[i,0]:=signal[i];
    end
  else
    begin
      for i:=1 to Xsize do line[i,0]:=exp(-sqr(i-Xsize/4)/sqr(0.3/dx));
    end;

  { Starting Condition Velocity }
  for i:=2 to Xsize-1 do
    begin
      dline[i]:=-vphase/RefractiveIndexAmbient*(line[i+1,0]-line[i-1,0])/(2*dx);
    end;

 { Starting Condition Ghost Point in Time }
  for i:=1 to Xsize do
    begin
      line[i,1]:=line[i,0]+0.5*sqr(alpha)*(line[i+1,0]-2*line[i,0]+line[i-1,0]) + dline[i]*dt;
    end;

  DtPt:=0;
  TimeTracePt:=0;
  repeat
    t:=t+dt;
    inc(DtPt);

    { inner nodes }
    for i:=2 to XSize-1 do
      begin
{       line[i,2]:= -line[i,0] + 2*line[i,1] + sqr(alpha)*(line[i+1,1]-2*line[i,1]+line[i-1,1]); }

        if i < trunc(boundary1/dx) then
          begin
            alpha1:=alpha*1/(RefractiveIndexAmbient);
            alpha2:=alpha*1/(RefractiveIndexAmbient);
            gamma0:=0;
          end
        else if i = trunc(boundary1/dx) then
          begin
            alpha1:=alpha*1/(RefractiveIndexAmbient);
            alpha2:=alpha*1/(RefractiveIndexPlug);
            gamma0:=0;
          end
        else if (i > trunc(boundary1/dx)) and (i < trunc(boundary2/dx)) then
          begin
            alpha1:=alpha*1/(RefractiveIndexPlug);
            alpha2:=alpha*1/(RefractiveIndexPlug);
            gamma0:=0;
          end
        else if i = trunc(boundary2/dx) then
          begin
            alpha1:=alpha*1/(RefractiveIndexplug);
            alpha2:=alpha*1/(RefractiveIndexMedium);
            gamma0:=0;
          end
        else if i > trunc(boundary2/dx) then
          begin
            alpha1:=alpha*1/(RefractiveIndexMedium);
            alpha2:=alpha*1/(RefractiveIndexMedium);
            gamma0:=gamma;
          end;

{       line[i,2]:= 1/(1+0.5*gamma*dt)*(
                             2*line[i,1] -line[i,0] +
                             gamma*dt/2*line[i,0]
                             + 0.5*c2*(alphavector[i+1]+alphavector[i])*(line[i+1,1]-line[i,1])
                             - 0.5*c2*(alphavector[i-1]+alphavector[i])*(line[i,1]-line[i-1,1])
                             ); }

        line[i,2]:= 1/(1+0.5*gamma0*dt)*(
                          2*line[i,1] -line[i,0] +
                          gamma0*dt/2*line[i,0]
                                     + sqr(alpha2)*(line[i+1,1]-line[i,1])
                                     - sqr(alpha1)*(line[i,1]-line[i-1,1])
                                       );


      end;
    { outer nodes }

    case BoundaryConditionLeft of
     0 : begin  { fixed ends }
          line[1,2]:=0;
        end;
     1 : begin  { open ends }
           alpha1:=alpha*1/(RefractiveIndexAmbient);
           line[1,2]:=-line[1,0] + 2*line[1,1] + sqr(alpha1)*(line[2,1]-line[1,1]);
         end;
     2 : begin  { absorbing ends }
           alpha1:=alpha*1/(RefractiveIndexAmbient);
           line[1,2]:=line[2,1] + (alpha1-1)/(alpha1+1)*(line[2,2]-line[1,1]);
         end;
      end;

    case BoundaryConditionRight of
      0 : begin  { fixed ends }
           line[Xsize,2]:=0;
         end;
      1 : begin  { open ends }
            alpha2:=alpha*1/(RefractiveIndexMedium);
            line[Xsize,2]:=-line[Xsize,0] + 2*line[Xsize,1] - sqr(alpha2)*(line[Xsize,1]-line[Xsize-1,1]);
          end;
      2 : begin  { absorbing ends }
            alpha2:=alpha*1/(RefractiveIndexMedium);
            line[Xsize,2]:=line[Xsize-1,1] + (alpha2-1)/(alpha2+1)*(line[Xsize-1,2]-line[XSize,1]);
          end;
      end;

    { write array back }
    for i:=1 to Xsize do
      begin
        line[i,0]:=line[i,1];
        line[i,1]:=line[i,2];
      end;
    if (DtPt mod DrawEachStep = 0) then
      begin
        DtPt:=0;

        Application.ProcessMessages;
        Chart1LineSeries1.Clear;
        for i:=1 to Xsize do
          begin
            Chart1LineSeries1.AddXY(i*dx,line[i,1]);
          end;
        Chart1.LeftAxis.Range.Min:=-2;
        Chart1.LeftAxis.Range.Max:=2;
        Chart1.Canvas.TextOut(100,20,FloatToStrF(t/1e-9,ffgeneral,5,0)+' ns');

        if SimulationTRUE=FALSE then
          begin
            Button1.Caption:='Start';
            closefile(datei);
            exit;
          end;
        Panel4.Caption:=' time: '+FloatToStrF(t/1e-9,ffgeneral,5,2)+ ' ns';

        writeln(datei,t/1e-9,' ',line[trunc(boundary2/dx),1]);

        if CheckBox1.Checked then
          begin
            inc(BmpPt);
            if (BmpPt mod DrawEachBmp = 0) then
              begin
                BmpPt:=0;
                inc(BMPCounter);
                DestRect.Top:=0;
                DestRect.Left:=0;
                DestRect.Bottom:=Chart1.Height;
                DestRect.Right:=Chart1.Width;
                SourceRect:=DestRect;
                Bitmap1.Canvas.CopyRect(DestRect,Chart1.Canvas,SourceRect);
                Bitmap1.SaveToFile('c:\temp\linetransmission\map'+IntToStr(BMPCounter)+'.bmp');
              end;
          end;

        { BCS signal }
         inc(TimeTracePt);
         BCSPT:=round(6/dx);
         Chart2LineSeries1.AddXY(t/1e-9,line[BCSPt,1]);
         TimeTraceBCS[TimeTracePt,1]:=t/1e-9;
         TimeTraceBCS[TimeTracePt,2]:=line[BCSPt,1];

        { Electrode signal }
         ElectrodePT:=round(boundary2/dx);
         Chart2LineSeries2.AddXY(t/1e-9,line[ElectrodePt,1]);
         TimeTraceElectrode[TimeTracePt,1]:=t/1e-9;
         TimeTraceElectrode[TimeTracePt,2]:=line[ElectrodePt,1];


      end;

  until (t>tend) or (SimulationTRUE=FALSE);
  TimeTraceNB:=TimeTracePt;
  closefile(datei);
  Bitmap1.Destroy;
end;


end.

