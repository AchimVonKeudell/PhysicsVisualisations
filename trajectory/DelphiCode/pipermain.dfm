object FormMain: TFormMain
  Left = 196
  Top = 93
  Width = 1023
  Height = 859
  Caption = 'Ion-Solid Interaction: Piper'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  Menu = MainMenu1
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object StatusBar1: TStatusBar
    Left = 0
    Top = 794
    Width = 1015
    Height = 19
    Panels = <
      item
        Width = 150
      end
      item
        Width = 50
      end>
    SimplePanel = False
  end
  object CoolBar1: TCoolBar
    Left = 0
    Top = 0
    Width = 1015
    Height = 37
    Bands = <
      item
        Control = ToolBar1
        ImageIndex = -1
        MinHeight = 27
        Width = 1011
      end>
    object ToolBar1: TToolBar
      Left = 9
      Top = 0
      Width = 998
      Height = 27
      AutoSize = True
      ButtonHeight = 25
      Caption = 'ToolBar1'
      EdgeBorders = []
      TabOrder = 0
      object SpeedButton2: TSpeedButton
        Left = 0
        Top = 2
        Width = 23
        Height = 25
        Glyph.Data = {
          76010000424D7601000000000000760000002800000020000000100000000100
          04000000000000010000120B0000120B00001000000000000000000000000000
          800000800000008080008000000080008000808000007F7F7F00BFBFBF000000
          FF0000FF000000FFFF00FF000000FF00FF00FFFF0000FFFFFF00333333333333
          33333333333333333333333333C3333333333333337F3333333333333C0C3333
          333333333777F33333333333C0F0C3333333333377377F333333333C0FFF0C33
          3333333777F377F3333333CCC0FFF0C333333373377F377F33333CCCCC0FFF0C
          333337333377F377F3334CCCCCC0FFF0C3337F3333377F377F33C4CCCCCC0FFF
          0C3377F333F377F377F33C4CC0CCC0FFF0C3377F3733F77F377333C4CCC0CC0F
          0C333377F337F3777733333C4C00CCC0333333377F773337F3333333C4CCCCCC
          3333333377F333F7333333333C4CCCC333333333377F37733333333333C4C333
          3333333333777333333333333333333333333333333333333333}
        NumGlyphs = 2
        OnClick = SpeedButton2Click
      end
      object ToolButton1: TToolButton
        Left = 23
        Top = 2
        Width = 8
        Caption = 'ToolButton1'
        Style = tbsSeparator
      end
      object Button1: TButton
        Left = 31
        Top = 2
        Width = 75
        Height = 25
        Caption = 'Start'
        TabOrder = 0
        OnClick = Button1Click
      end
      object SpeedButton1: TSpeedButton
        Left = 106
        Top = 2
        Width = 70
        Height = 25
        Caption = 'Stop'
        OnClick = SpeedButton1Click
      end
      object ToolButton2: TToolButton
        Left = 176
        Top = 2
        Width = 8
        Caption = 'ToolButton2'
        ImageIndex = 0
        Style = tbsSeparator
      end
      object SpeedButton3: TSpeedButton
        Left = 184
        Top = 2
        Width = 105
        Height = 25
        AllowAllUp = True
        GroupIndex = 1
        Caption = 'Parameter...'
        OnClick = SpeedButton3Click
      end
    end
  end
  object Panel1: TPanel
    Left = 0
    Top = 37
    Width = 1015
    Height = 757
    Align = alClient
    BevelInner = bvLowered
    BevelOuter = bvNone
    Caption = 'Panel1'
    TabOrder = 2
    object Splitter1: TSplitter
      Left = 1
      Top = 413
      Width = 1013
      Height = 8
      Cursor = crVSplit
      Align = alTop
      Beveled = True
    end
    object Plot1: TPlot
      Left = 1
      Top = 1
      Width = 1013
      Height = 412
      Color = clWhite
      Align = alTop
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Arial'
      Font.Style = []
      GrafikLeft = 10
      GrafikTop = 5
      GrafikWidth = 80
      GrafikHeight = 80
      GrafikColor = clBtnFace
      AxisColor = clBlack
      Spectra.NB = 4
      Spectra.Spectrum = (
        4
        1
        1
        False
        False
        True
        1
        0
        255
        'Times New Roman'
        14
        True
        True
        False
        False
        'atom distances'
        '-5'
        '100'
        '-5'
        '100'
        'Spectrum_1'
        'Legende Spectrum_1'
        False
        2
        1
        False
        True
        True
        3
        0
        255
        'Times New Roman'
        14
        True
        False
        False
        False
        'atom distances'
        '-20'
        '20'
        '-20'
        '20'
        'Spectrum_2'
        'Legende Spectrum_2'
        False
        3
        3
        False
        False
        False
        1
        0
        255
        'Times New Roman'
        14
        True
        True
        False
        False
        'atom distances'
        '-5'
        '100'
        '-5'
        '100'
        'Spectrum_3'
        'Legende Spectrum_3'
        False
        4
        3
        False
        True
        False
        3
        0
        16711680
        'Times New Roman'
        14
        True
        False
        False
        False
        'N(z)'
        '-20'
        '20'
        '-20'
        '20'
        'Spectrum_4'
        'Legende Spectrum_4'
        False)
      OnScaleChange = Plot1ScaleChange
    end
    object Plot2: TPlot
      Left = 1
      Top = 421
      Width = 1013
      Height = 335
      Color = clWhite
      Align = alClient
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Arial'
      Font.Style = []
      GrafikLeft = 10
      GrafikTop = 5
      GrafikWidth = 80
      GrafikHeight = 80
      GrafikColor = clBtnFace
      AxisColor = clBlack
      Spectra.NB = 2
      Spectra.Spectrum = (
        2
        1
        1
        False
        False
        True
        1
        0
        255
        'Times New Roman'
        14
        True
        True
        False
        False
        'atom distances'
        '-5'
        '100'
        '-5'
        '100'
        'Spectrum_1'
        'Legende Spectrum_1'
        False
        2
        1
        False
        True
        True
        3
        0
        16711680
        'Times New Roman'
        14
        False
        False
        True
        False
        'N(z)'
        '0'
        '10'
        '0'
        '10'
        'Spectrum_2'
        'Legende Spectrum_2'
        False)
    end
  end
  object MainMenu1: TMainMenu
    Left = 136
    Top = 81
    object File1: TMenuItem
      Caption = '&File'
    end
    object N1: TMenuItem
      Caption = '?'
      OnClick = N1Click
    end
  end
end
