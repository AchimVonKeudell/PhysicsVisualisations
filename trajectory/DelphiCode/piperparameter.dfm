object FormParameter: TFormParameter
  Left = 625
  Top = 224
  Width = 279
  Height = 319
  Caption = 'Parameter...'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  FormStyle = fsStayOnTop
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 271
    Height = 105
    Align = alTop
    Caption = ' Projectile '
    TabOrder = 0
    object Label1: TLabel
      Left = 152
      Top = 20
      Width = 14
      Height = 13
      Caption = 'm1'
    end
    object Label2: TLabel
      Left = 156
      Top = 48
      Width = 11
      Height = 13
      Caption = 'z1'
    end
    object Label3: TLabel
      Left = 112
      Top = 76
      Width = 54
      Height = 13
      Caption = 'energy (eV)'
    end
    object Editm1: TEdit
      Left = 176
      Top = 16
      Width = 73
      Height = 21
      TabOrder = 0
      Text = 'Editm1'
    end
    object Editz1: TEdit
      Left = 176
      Top = 44
      Width = 73
      Height = 21
      TabOrder = 1
      Text = 'Editz1'
    end
    object Edite: TEdit
      Left = 176
      Top = 72
      Width = 73
      Height = 21
      TabOrder = 2
      Text = 'Edite'
    end
  end
  object GroupBox2: TGroupBox
    Left = 0
    Top = 105
    Width = 271
    Height = 105
    Align = alTop
    Caption = ' Target '
    TabOrder = 1
    object Label4: TLabel
      Left = 152
      Top = 20
      Width = 14
      Height = 13
      Caption = 'm2'
    end
    object Label5: TLabel
      Left = 156
      Top = 48
      Width = 11
      Height = 13
      Caption = 'z2'
    end
    object Editz2: TEdit
      Left = 176
      Top = 44
      Width = 73
      Height = 21
      TabOrder = 0
      Text = 'Edit2'
    end
    object Editm2: TEdit
      Left = 176
      Top = 16
      Width = 73
      Height = 21
      TabOrder = 1
      Text = 'Edit1'
    end
  end
  object GroupBox3: TGroupBox
    Left = 0
    Top = 210
    Width = 271
    Height = 79
    Align = alTop
    Caption = ' Simulation '
    TabOrder = 2
    object Label6: TLabel
      Left = 112
      Top = 48
      Width = 57
      Height = 13
      Caption = 'NB particles'
    end
    object Label7: TLabel
      Left = 84
      Top = 20
      Width = 84
      Height = 13
      Caption = 'cutoff energy (eV)'
    end
    object Editecutoff: TEdit
      Left = 176
      Top = 16
      Width = 73
      Height = 21
      TabOrder = 0
      Text = 'Editecutoff'
    end
    object EditNB: TEdit
      Left = 176
      Top = 44
      Width = 73
      Height = 21
      TabOrder = 1
      Text = 'EditNB'
    end
  end
end
