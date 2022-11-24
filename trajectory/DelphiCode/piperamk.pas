unit piperamk;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, jpeg, ExtCtrls;

type
  TFormAbout = class(TForm)
    Image1: TImage;
    Label1: TLabel;
    Button1: TButton;
    procedure Button1Click(Sender: TObject);
  private
    { Private-Deklarationen }
  public
    { Public-Deklarationen }
  end;

var
  FormAbout: TFormAbout;

implementation

{$R *.DFM}

procedure TFormAbout.Button1Click(Sender: TObject);
begin
  close;
end;

end.
