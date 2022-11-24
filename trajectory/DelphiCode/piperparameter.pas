unit piperparameter;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls;

type
  TFormParameter = class(TForm)
    GroupBox1: TGroupBox;
    GroupBox2: TGroupBox;
    Editm1: TEdit;
    Editz1: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Edite: TEdit;
    Editz2: TEdit;
    Editm2: TEdit;
    Label4: TLabel;
    Label5: TLabel;
    GroupBox3: TGroupBox;
    Editecutoff: TEdit;
    EditNB: TEdit;
    Label6: TLabel;
    Label7: TLabel;
    procedure FormCreate(Sender: TObject);
  private
    { Private-Deklarationen }
  public
    procedure SetParameter;
    procedure GetParameter;
    { Public-Deklarationen }
  end;

var
  FormParameter: TFormParameter;

implementation

uses piperdat;

procedure TFormParameter.SetParameter;
begin
  Editz1.Text:=FloatToStrF(z1,fffixed,5,3);
  Editm1.Text:=FloatToStrF(m1,fffixed,5,3);
  Editz2.Text:=FloatToStrF(z2,fffixed,5,3);
  Editm2.Text:=FloatToStrF(m2,fffixed,5,3);
  Edite.Text:=FloatToStrF(energy,fffixed,5,3);

  Editecutoff.Text:=FloatToStrF(ecutoff,fffixed,5,3);
  EditNB.Text:=IntToStr(NB);
end;

procedure TFormParameter.GetParameter;
begin
  z1:=StrToFloat(Editz1.Text);
  m1:=StrToFloat(Editm1.Text);
  z2:=StrToFloat(Editz2.Text);
  m2:=StrToFloat(Editm2.Text);
  energy:=StrToFloat(Edite.Text);

  ecutoff:=StrToFloat(Editecutoff.Text);
  NB:=StrToInt(EditNB.Text);
end;

{$R *.DFM}

procedure TFormParameter.FormCreate(Sender: TObject);
begin
  SetParameter;
end;

end.
