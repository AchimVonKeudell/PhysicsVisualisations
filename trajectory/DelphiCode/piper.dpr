program piper;

uses
  Forms,
  pipermain in 'pipermain.pas' {FormMain},
  piperdat in 'piperdat.pas',
  piperparameter in 'piperparameter.pas' {FormParameter},
  piperamk in 'piperamk.pas' {FormAbout};

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TFormMain, FormMain);
  Application.CreateForm(TFormParameter, FormParameter);
  Application.CreateForm(TFormAbout, FormAbout);
  Application.Run;
end.
