function varargout=sc1_sc2_imana(what,varargin)

% UPDATE INFO IN SECTION 2 WITH EVERY NEW SUBJECT, SESSION

numDummys = 3;                                                              % per run
numTRs    = 601;                                                            % per run (includes dummies)

%========================================================================================================================
% (1) Directories
baseDir          = '/Users/maedbhking/Documents/Cerebellum_Cognition';
% baseDir          = '/Volumes/Seagate Backup Plus Drive';
% baseDir        = '/Volumes/MotorControl/data/super_cerebellum_new';
% baseDir        = '/Users/jdiedrichsen/Data/super_cerebellum_new';

studyDir{1}     =fullfile(baseDir,'sc1');
studyDir{2}     =fullfile(baseDir,'sc2');
restDir         =fullfile(baseDir,'restingState');
behavDir        ='/data';
imagingDir      ='/imaging_data';
imagingDirNA    ='/imaging_data_nonaggr';
imagingDirA     ='/imaging_data_aggr';
imagingDirRaw   ='/imaging_data_raw';
dicomDir        ='/imaging_data_dicom';
anatomicalDir   ='/anatomicals';
freesurferDir   ='/surfaceFreesurfer';
suitDir         ='/suit';
caretDir        ='/surfaceCaret';
regDir          ='/RegionOfInterest/';
eyeDirRaw       ='/eyeTracking_raw';
eyeDir          ='/eyeTracking';
physioDir       ='/physio';
encodeDir       ='/encoding';

%========================================================================================================================
% (2) PRE-PROCESSING: Study AND Subject AND Session specific info. EDIT W/ EVERY NEW SUBJECT AND W/ EVERY NEW
% SESSION
% Subjs 1-22 were preprocessed using this script. The naming convention for
% the dicom files changed after s22 so subjs 23-31 were preprocessed with
% sc1_sc2_imana_ch (exact same preprocessing steps)
% this script was used to process resting state data

% MANUAL INPUT REQUIRED FOR THE FOLLOWING VARIABLES:
% DicomName:     find this in imaging_data_dicom
% rs_DicomName:  find this in restingState
% NiiRawName:    find this in imaging_data_dicom (numbers before series number)
% rs_NiiRawName: find this in restingState (numbers before series number)
% saccRun:       runs to be included for saccade analysis
% anatNum:       one scan per subject (collected in session 1 of study 1 - usually series 2)

subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};

% returnSubjs=[2,3,4,6,7,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];
returnSubjs=[2,3,4,6,8,9,10,12,14,15,17,18,19,20,21,22,24,25,26,27,28,29,30,31];

DicomName{1}  = {'2016_03_15_HP10.MR.DIEDRICHSEN_FEB2016',...%s01
    '2016_03_23_HP10.MR.DIEDRICHSEN_FEB2016',... %s01
    '2016_04_13_I011.MR.DIEDRICHSEN_FEB2016',... %s02
    '2016_04_156_I011.MR.DIEDRICHSEN_FEB2016',...%s02
    '2016_04_18_UH11.MR.DIEDRICHSEN_FEB2016',... %s03
    '2016_04_19_UH11.MR.DIEDRICHSEN_FEB2016',... %s03
    '2016_04_26_CC23.MR.DIEDRICHSEN_FEB2016',... %s04
    '2016_04_29_CC23.MR.DIEDRICHSEN_FEB2016',... %s04
    '2016_05_18_EK05.MR.DIEDRICHSEN_FEB2016',... %s05
    '2016_05_19_EK05.MR.DIEDRICHSEN_FEB2016',... %s05
    '2016_05_24_M028.MR.DIEDRICHSEN_FEB2016',... %s06
    '2016_05_30_M028.MR.DIEDRICHSEN_FEB2016',... %s06
    '2016_05_25_HI15.MR.DIEDRICHSEN_FEB2016',... %s07
    '2016_05_26_HI15.MR.DIEDRICHSEN_FEB2016',... %s07
    '2016_05_31_NL20.MR.DIEDRICHSEN_FEB2016',... %s08
    '2016_06_01_NL20.MR.DIEDRICHSEN_FEB2016',... %s08
    '2016_06_02_OH14.MR.DIEDRICHSEN_FEB2016',... %s09
    '2016_06_03_OH14.MR.DIEDRICHSEN_FEB2016',... %s09
    '2016_07_07_RV14.MR.DIEDRICHSEN_FEB2016',... %s10
    '2016_07_08_RV14.MR.DIEDRICHSEN_FEB2016',... %s10
    '2016_07_15_AC15.MR.DIEDRICHSEN_FEB2016',... %s11
    '2016_07_22_AC15.MR.DIEDRICHSEN_FEB2016',... %s11
    '2016_07_20_RM20.MR.DIEDRICHSEN_FEB2016',... %s12
    '2016_07_21_RM20.MR.DIEDRICHSEN_FEB2016',... %s12
    '2016_07_28_AW12.MR.DIEDRICHSEN_FEB2016',... %s13
    '2016_07_29_AW12.MR.DIEDRICHSEN_FEB2016',... %s13
    '2016_07_28_SR12.MR.DIEDRICHSEN_FEB2016',... %s14
    '2016_07_29_SR12.MR.DIEDRICHSEN_FEB2016',... %s14
    '2016_08_04_RT25.MR.DIEDRICHSEN_FEB2016',... %s15
    '2016_08_05_RT25.MR.DIEDRICHSEN_FEB2016',... %s15
    '2016_08_04_FW03.MR.DIEDRICHSEN_FEB2016',... %s16
    '2016_08_05_FW03.MR.DIEDRICHSEN_FEB2016',... %s16
    '2016_08_22_AR14.MR.DIEDRICHSEN_FEB2016',... %s17
    '2016_08_23_AR14.MR.DIEDRICHSEN_FEB2016',... %s17
    '2017_02_14_AO21.MR.DIEDRICHSEN_MAEBDH',...  %s18
    '2017_02_22_AO21.MR.DIEDRICHSEN_MAEBDH',...  %s18
    '2017_03_09_HA31.MR.DIEDRICHSEN_MAEBDH',...  %s19
    '2017_03_13_HA31.MR.DIEDRICHSEN_MAEBDH',...  %s19
    '2017_03_15_ML24.MR.DIEDRICHSEN_MAEBDH',...  %s20
    '2017_03_20_ML24.MR.DIEDRICHSEN_MCC',...     %s20
    '2017_03_13_AZ18.MR.DIEDRICHSEN_MAEBDH',...  %s21
    '2017_03_20_AZ18.MR.DIEDRICHSEN_MCC',...     %s21
    '2017_03_15_VE14.MR.DIEDRICHSEN_MAEBDH',...  %s22
    '2017_03_16_VE14.MR.DIEDRICHSEN_MCC',...     %s22
    };

DicomName{2} = {'',...%s01
    '',...             %s01
    '2017_05_25_IO11.MR.DIEDRICHSEN_MCC',...%s02
    '2017_05_26_IO11.MR.DIEDRICHSEN_MCC',...%s02
    '2017_03_27_UH11.MR.DIEDRICHSEN_MCC',...%s03
    '2017_03_28_UH11.MR.DIEDRICHSEN_MCC',...%s03
    '2017_06_12_CC23.MR.DIEDRICHSEN_MCC',...%s04
    '2017_06_13_CC23.MR.DIEDRICHSEN_MCC',...%s04
    '',...%s05
    '',...%s05
    '2017_04_12_MO28.MR.DIEDRICHSEN_MCC',...%s06
    '2017_04_20_MO28.MR.DIEDRICHSEN_MCC',...%s06
    '2017_05_16_HI15.MR.DIEDRICHSEN_MCC',...%s07
    '2017_05_17_HI15.MR.DIEDRICHSEN_MCC',...%s07
    '2017_04_24_NL20.MR.DIEDRICHSEN_MCC',...%s08
    '2017_05_01_NL20.MR.DIEDRICHSEN_MCC',...%s08
    '2017_04_19_OH14.MR.DIEDRICHSEN_MCC',...%s09
    '2017_04_20_OH14.MR.DIEDRICHSEN_MCC',...%s09
    '2017_04_27_RV14.MR.DIEDRICHSEN_MCC',...%s10
    '2017_04_28_RV14.MR.DIEDRICHSEN_MCC',...%s10
    '',...%s11
    '',...%s11
    '2017_04_24_RM20.MR.DIEDRICHSEN_MCC',...%s12
    '2017_04_25_RM20.MR.DIEDRICHSEN_MCC',...%s12
    '',...%s13
    '',...%s13
    '2017_06_12_SR12.MR.DIEDRICHSEN_MCC',...%s14
    '2017_06_15_SR12.MR.DIEDRICHSEN_MCC',...%s14
    '2017_06_14_RT25.MR.DIEDRICHSEN_MCC',...%s15
    '2017_06_15_RT25.MR.DIEDRICHSEN_MCC',...%s15
    '',...%s16
    '',...%s16
    '2017_03_24_AR14.MR.DIEDRICHSEN_MCC',...%s17
    '2017_03_28_AR14.MR.DIEDRICHSEN_MCC',...%s17
    '2017_05_09_AO21.MR.DIEDRICHSEN_MCC',...%s18
    '2017_05_10_AO21.MR.DIEDRICHSEN_MCC',...%s18
    '2017_05_25_HA31.MR.DIEDRICHSEN_MCC',...%s19
    '2017_05_26_HA31.MR.DIEDRICHSEN_MCC',...%s19
    '2017_04_12_ML24.MR.DIEDRICHSEN_MCC',...%s20
    '2017_04_23_ML24.MR.DIEDRICHSEN_MCC',...%s20
    '2017_04_25_AZ18.MR.DIEDRICHSEN_MCC',...%s21
    '2017_04_28_AZ18.MR.DIEDRICHSEN_MCC',...%s21
    '2017_04_26_VE14.MR.DIEDRICHSEN_MCC',...%s22
    '2017_04_28_VE14.MR.DIEDRICHSEN_MCC',...%s22
    };

rs_DicomName = {'',... %s01
    '2017_09_07_IO11.MR.DIEDRICHSEN_MCC',...%s02
    '2017_09_13_UH11.MR.DIEDRICHSEN_MCC',...%s03
    '',...%s04
    '',...%s05
    '2017_08_24_MO28.MR.DIEDRICHSEN_MCC',...%s06
    '2017_09_05_HI15.MR.DIEDRICHSEN_MCC',...%s07
    '2017_09_07_NL20.MR.DIEDRICHSEN_MCC',...%s08
    '',...%s09
    '2017_08_15_RV14.MR.DIEDRICHSEN_MCC',...%s10
    '',...%s11
    '2017_08_11_RM20.MR.DIEDRICHSEN_MCC',...%s12
    '',...%s13
    '2017_09_07_SR12.MR.DIEDRICHSEN_MCC',...%s14
    '',...%s15
    '',...%s16
    '',...%s17
    '2017_09_07_AO21.MR.DIEDRICHSEN_MCC',...%s18
    '',...%s19
    '2017_08_24_ML24.MR.DIEDRICHSEN_MCC',...%s20
    '',...%s21
    '2017_09_07_VE14.MR.DIEDRICHSEN_MCC',...%s22
    '2017_09_12_HO21.MR.Diedrichsen_MCC',...%s23
    '2017_11_01_NA14.MR.Diedrichsen_MCC',...%s24
    '2017_09_12_DA18.MR.Diedrichsen_MCC',...%s25
    '2017_11_03_ME01.MR.Diedrichsen_MCC',...%s26
    '2017_10_25_VI14.MR.Diedrichsen_MCC',...%s27
    '2017_11_24_ME06.MR.Diedrichsen_MCC',...%s28
    '2017_11_02_MA19.MR.Diedrichsen_MCC',...%s29
    '2017_11_01_AV10.MR.Diedrichsen_MCC',...%s30
    '2017_11_03_MA18.MR.Diedrichsen_MCC',...%s31
    };


NiiRawName{1} = {'160315140823DST131221107524367007',...%s01
    '160323140933DST131221107524367007',...%s01
    '160413101535DST131221107524367007',...%s02
    '160415100919DST131221107524367007',...%s02
    '160418140549DST131221107524367007',...%s03
    '160419152813DST131221107524367007',...%s03
    '160426135807DST131221107524367007',...%s04
    '160429150526DST131221107524367007',...%s04
    '160518145625DST131221107524367007',...%s05
    '160519140503DST131221107524367007',...%s05
    '160524143451DST131221107524367007',...%s06
    '160530143056DST131221107524367007',...%s06
    '160525145854DST131221107524367007',...%s07
    '160526150947DST131221107524367007',...%s07
    '160531141052DST131221107524367007',...%s08
    '160601140606DST131221107524367007',...%s08
    '160602145638DST131221107524367007',...%s09
    '160603144457DST131221107524367007',...%s09
    '160707141119DST131221107524367007',...%s10
    '160708125723DST131221107524367007',...%s10
    '160715123900DST131221107524367007',...%s11
    '160722143518DST131221107524367007',...%s11
    '160720110707DST131221107524367007',...%s12
    '160721105741DST131221107524367007',...%s12
    '160728100746DST131221107524367007',...%s13
    '160729100321DST131221107524367007',...%s13
    '160728120942DST131221107524367007',...%s14
    '160729120849DST131221107524367007',...%s14
    '160804140003DST131221107524367007',...%s15
    '160805130528DST131221107524367007',...%s15
    '160804100807DST131221107524367007',...%s16
    '160805100423DST131221107524367007',...%s16
    '160822095244DST131221107524367007',...%s17
    '160823100118DST131221107524367007',...%s17
    '170214144242STD131221107524367007',...%s18
    '170222143758STD131221107524367007',...%s18
    '170309121041STD131221107524367007',...%s19
    '170313120111DST131221107524367007',...%s19
    '170315121648DST131221107524367007',...%s20
    '170320120712DST131221107524367007',...%s20
    '170313095950DST131221107524367007',...%s21
    '170320100001DST131221107524367007',...%s21
    '170315144047DST131221107524367007',...%s22
    '170316144945DST131221107524367007',...%s22
    };


NiiRawName{2} =  {'',...%s01
    '',...           %s01
    '170525100719DST131221107524367007',...%s02
    '170525100719DST131221107524367007',...%s02
    '170327151145DST131221107524367007',...%s03
    '170328124001DST131221107524367007',...%s03
    '170612140815DST131221107524367007',...%s04
    '170612140815DST131221107524367007',...%s04
    '',...%s05
    '',...%s05
    '170412101047DST131221107524367007',...%s06
    '170420150758DST131221107524367007',...%s06
    '2017_05_16_HI15',...%s07 (strange formatting for this subject)
    '2017_05_16_HI15',...%s07
    '170424100628DST131221107524367007',...%s08
    '170501083431DST131221107524367007',...%s08
    '170419115131DST131221107524367007',...%s09
    '170420120648DST131221107524367007',...%s09
    '170427131043DST131221107524367007',...%s10
    '170427131043DST131221107524367007',...%s10
    '',...%s11
    '',...%s11
    '170424140431DST131221107524367007',...%s12
    '170424140431DST131221107524367007',...%s12
    '',...%s13
    '',...%s13
    '170612120235DST131221107524367007',...%s14
    '170615145716DST131221107524367007',...%s14
    '170614140812DST131221107524367007',...%s15
    '170614140812DST131221107524367007',...%s15
    '',...%s16
    '',...%s16
    '170324151234DST131221107524367007',...%s17
    '170328100924DST131221107524367007',...%s17
    '170509140556DST131221107524367007',...%s18
    '170509140556DST131221107524367007',...%s18
    '170525151250DST131221107524367007',...%s19
    '170525151250DST131221107524367007',...%s19
    '170412121416DST131221107524367007',...%s20
    '170417125335DST131221107524367007',...%s20
    '170425155510DST131221107524367007',...%s21
    '170425155510DST131221107524367007',...%s21
    '170427102053DST131221107524367007',...%s22
    '170427102053DST131221107524367007',...%s22
    };

rs_NiiRawName = {'',...%s01
    'IO11',...%s02
    'UH11',...%s03
    '',...%s04
    '',...%s05
    'MO28',...%s06
    '2017_09_05_HI15',...%s07 (strange formatting for this subject)
    'NL20',...%s08
    '',...%s09
    'RV14',...%s10
    '',...%s11
    'RM20',...%s12
    '',...%s13
    'SR12',...%s14
    '',...%s15
    '',...%s16
    '',...%s17
    'AO21',...%s18
    '',...%s19
    'ML24',...%s20
    '',...%s21
    'VE14',...%s22
    'HO21',...%s23
    '171101121712DST131221107524367007',...%s24
    'DA18',...%s25
    '171103143826DST131221107524367007',...%s26
    'VI14',...%s27
    '171124115601STD131221107524367007',...%s28
    'MA19',...%s29
    '171101113420DST131221107524367007',...%s30
    '171103100256DST131221107524367007',...%s31
    };

fscanNum{1} = {[3,4,5,6,8,9,10,11],...%s01
    [2,3,4,5,6,7,8,9],...    %s01
    [3,4,5,6,8,10,11,12],... %s02
    [2,3,4,5,6,7,8,9],...    %s02
    [3,4,5,6,7,8,9,10],...   %s03
    [3,4,5,6,8,9,10,11],...  %s03
    [4,7,8,9,10,11,12,13],...%s04
    [2,3,4,5,6,7,8,9],...    %s04
    [3,4,5,6,7,8,9,10],...   %s05
    [2,3,4,5,6,7,8,9],...    %s05
    [3,4,5,6,7,8,9,10],...   %s06
    [2,3,4,5,6,7,8,9],...    %s06
    [3,4,5,6,7,8,9,10],...   %s07
    [2,3,4,5,6,7,8,9],...    %s07
    [3,4,5,6,7,8,9,10],...   %s08
    [2,3,4,5,6,7,8,9],...    %s08
    [3,4,5,6,7,8,9,10],...   %s09
    [2,3,4,5,6,7,8,9],...    %s09
    [3,4,5,6,7,8,9,10],...   %s10
    [2,3,4,5,6,7,8,9],...    %s10
    [4,5,6,7,8,9,10,11],...  %s11
    [3,4,5,6,7,8,9,10],...   %s11
    [3,4,5,7,8,9,10,11],...  %s12
    [2,3,4,5,6,7,8,9],...    %s12
    [4,5,6,7,8,9,10,11],...  %s13
    [2,3,4,5,6,7,9,10],...   %s13
    [4,5,6,7,8,9,10,11],...  %s14
    [2,3,4,5,6,7,8,9],...    %s14
    [3,4,5,6,7,8,9,10],...   %s15
    [2,3,4,6,7,8,9,10],...   %s15
    [3,4,5,6,7,8,9,10],...   %s16
    [2,3,4,5,6,7,8,9],...    %s16
    [3,4,5,6,7,8,9,10],...   %s17
    [2,3,4,5,6,7,8,9],...    %s17
    [3,4,5,6,7,8,9,10],...   %s18
    [2,4,5,6,7,8,10,11],...  %s18
    [3,4,5,6,7,8,9,10],...   %s19
    [2,3,4,5,6,7,8,9],...    %s19
    [3,4,5,6,7,8,9,10],...   %s20
    [2,3,4,5,6,7,8,9],...    %s20
    [3,4,5,6,7,9,10,11],...  %s21
    [3,4,5,7,8,9,10,11],...  %s21
    [3,4,5,6,7,8,9,10],...   %s22
    [2,3,4,5,6,7,8,9],...    %s22
    };

fscanNum{2}   = {[],...%s01 functional scans (series number)
    [],...%s01
    [3,5,7,9,11,13,15,17],...%s02
    [2,4,6,8,10,12,14,16],...%s02
    [3,5,7,9,11,13,15,17],...% s03_1
    [2,4,6,8,10,12,14,16],...% s03_2
    [3,5,7,9,11,13,15,17],...%s04
    [2,4,6,8,10,12,14,16],...%s04
    [],...%s05
    [],...%s05
    [4,6,8,10,12,14,16,18],...%s06_1
    [2,4,6,8,10,12,14,16],...%s06_2
    [4,6,8,10,12,14,16,18],...%s07_1
    [2,4,6,8,10,12,14,16],...% s07_2
    [3,5,7,9,11,13,15,17],...% s08_1
    [2,4,6,8,10,12,14,16],...% s08_2
    [3,5,7,9,11,13,15,17],...% s09_1
    [2,4,6,8,10,12,14,16],...% s09_2
    [3,5,7,9,11,13,15,17],...% s10_1
    [3,5,7,9,11,13,15,17],...% s10_2
    [],...%s11
    [],...%s11
    [3,5,7,9,11,13,15,17],...% s12_1
    [2,4,6,8,10,12,14,16],...% s12_2
    [],...%s13
    [],...%s13
    [3,5,7,9,11,13,15,17],...%s14
    [2,4,6,8,10,12,14,16],...%s14
    [3,5,7,9,11,13,15,17],...%s15
    [2,4,6,8,10,12,14,16],...%s15
    [],...%s16
    [],...%s16
    [3,5,7,9,11,13,15,17],...% s17_1
    [2,4,6,8,10,12,14,16],...% s17_2
    [3,5,7,9,11,13,15,17],...% s18_1
    [2,4,6,8,10,12,14,16],...% s18_2
    [3,5,7,9,11,13,15,17],...% s19_1
    [2,4,6,8,10,12,14,16],...% s19_2
    [3,5,7,9,11,13,15,17],...% s20_1
    [2,4,6,8,10,12,14,16],...% s20_2
    [3,5,7,9,11,13,15,17],...% s21_1
    [2,4,6,8,10,12,14,16],...% s21_2
    [3,5,7,9,11,13,15,17],...% s22_1
    [2,4,6,8,10,12,14,16],...% s22_2
    };

rs_fscanNum   = {[],...%s01 functional scans (series number)
    [2,4],...%s02
    [2,4],...%s03
    [],...%s04
    [],...%s05
    [2,4],...%s06
    [2,4],...%s07
    [2,4],...%s08
    [],...%s09
    [2,4],...%s10
    [],...%s11
    [2,4],...% s12
    [],...%s13
    [2,4],...%s14
    [],...%s15
    [],...%s16
    [],...%s17
    [2,4],...%s18
    [],...%s19
    [2,4],...% s20
    [],...% s21
    [2,4],...% s22
    [2,4],...%s23
    [2,4],...%s24
    [2,4],...%s25
    [23,25],...%s26
    [2,4],...%s27
    [20,22],...%s28
    [],...%s29 % need to fix this one
    [3,5],...%s30
    [2,4],...%s31
    };


saccRun{1} = {[],...%s01
    [1],...      %s02
    [],...       %s03
    [],...       %s04
    [1,2,3,4],...%s05
    [1,2,3,4],...%s06
    [3,4],...    %s07
    [1,2,3],...  %s08
    [1,2,3,4],...%s09
    [3,4],...    %s10
    [2,3,4],...  %s11
    [1,2],...    %s12
    [2],...      %s13
    [1,2,3],...  %s14
    [1,2],...    %s15
    [2],...      %s16
    [1,2],...    %s17
    [1,2],...    %s18
    [1,2],...    %s19
    [1],...      %s20
    [1,2],...    %s21
    [1],...      $s22
    };

saccRun{2} =   {[],... %s01
    [],...   %s02 [1,2]
    [],...   %s03
    [],...   %s04
    [],...   %s05
    [1,2],...%s06
    [1,2],...%s07
    [],...   %s08
    [2],...  %s09
    [1,2],...%s10
    [],...   %s11
    [],...   %s12
    [],...   %s13
    [1,2],...%s14
    [],...   %s15
    [],...   %s16
    [],...   %s17
    [1,2],...%s18
    [1,2],...%s19
    [1,2],...%s20
    [1,2],...%s21
    [1,2],...%s22
    };

anatNum = {7,...%s01
    2,...       %s02
    2,...       %s03
    2,...       %s04
    2,...       %s05
    2,...       %s06
    2,...       %s07
    2,...       %s08
    2,...       %s09
    2,...       %s10
    2,...       %s11
    2,...       %s12
    2,...       %s13
    3,...       %s14
    2,...       %s15
    2,...       %s16
    2,...       %s17
    2,...       %s18
    2,...       %s19
    2,...       %s20
    2,...       %s21
    2,...       %s22
    };

% The values of loc_AC should be acquired manually prior to the preprocessing
%   Step 1: open .nii file with MRIcron and manually find AC and read the xyz coordinate values
%   (note: these values are not [0 0 0] in the MNI coordinate)
%   Step 2: set those values into loc_AC (subtract from zero)

loc_AC = {[-108,-128,-149],...%s01
    [-84, -124,-131],...      %s02
    [-93,-132,-151],...       %s03
    [-91,-136,-147],...       %s04
    [-89, -133,-125],...      %s05
    [-93,-139,-136],...       %s06
    [-90,-130,-149],...       %s07
    [-90,-132,-145],...       %s08
    [-98,-141,-127],...       %s09
    [-93,-132,-122],...       %s10
    [-102,-132,-137],...      %s11
    [-96,-146,-147],...       %s12
    [-95,-141,-143],...       %s13
    [-90,-138,-146],...       %s14
    [-98,-147,-142],...       %s15
    [-93,-126,-157],...       %s16
    [-92,-131,-134],...       %s17
    [-94,-139,-137],...       %s18
    [-92,-137,-127],...       %s19
    [-95,-140,-152],...       %s20
    [-94,-139,-137],...       %s21
    [-94,-143,-145],...       %s22
    };
%========================================================================================================================
% (3) GLM.

% GLM Directories - change glmNum when appropriate

funcRunNum = [51,66];  % first and last behavioural run numbers (16 runs per subject)

run        = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16'};

runB = [51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66];  % Behavioural labelling of runs

sess = [1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2]; % session number

runs{1}{1}=1:8;   % study 1, session 1
runs{1}{2}=9:16;  % study 1, session 2
runs{2}{1}=17:24; % study 2, session 1
runs{2}{2}=25:32; % study 2, session 2

PW = {'betasNW','betasUW'}; % methods of noise normalisation (none and univariate noise normalisation)

% (4) Hemisphere and Region Names
hem        = {'lh','rh'};
regSide    = 1:2;
atlasA    = 'x';
atlasname = 'fsaverage_sym';
hemName   = {'LeftHem','RightHem'};

switch(what)
    
    case 'ANAT:dicom_import'                 % STEP 1.1: Import anat dicom
        % STUDY 1 ONLY
        % converts dicom to nifti files w/ spm_dicom_convert
        % anatomical always collected in study 1, sess 1
        % example: sc1_sc2_imana('ANAT:dicom_import',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            dircheck(fullfile(studyDir{1},dicomDir,[subj_name{sn(s)},'_1']));
            cd(fullfile(studyDir{1},dicomDir,[subj_name{sn(s)},'_1']));
            
            for i=1:length(anatNum{sn(s)})
                r=anatNum{sn(s)}(i);
                ss=(sn(s)*2)-1;
                DIR=dir(sprintf('%s.%4.4d.*.IMA',DicomName{ss},r)); % Get DICOM FILE NAMES
                Names=vertcat(DIR.name);
                if (~isempty(Names))
                    HDR=spm_dicom_headers(Names,1);  % Load dicom headers
                    dirname{r}=sprintf('series%2.2d',r);
                    if (~exist(dirname{r}))
                        mkdir(dirname{r}); %make dir for series{r} for .nii file output
                    end;
                    dircheck(fullfile(studyDir{1},dicomDir,[subj_name{sn(s)},'_1'],dirname{r}));
                    cd(fullfile(studyDir{1},dicomDir,[subj_name{sn(s)},'_1'],dirname{r}));
                    spm_dicom_convert(HDR,'all','flat','nii');%,dirname{r});  % Convert the data to nifti
                    cd ..
                end;
                fprintf('Series %d done \n',anatNum{sn(s)}(i))
            end;
            fprintf('Anatomical runs have been imported for %s Rename the .nii file as ''anatomical_raw.nii'' and move to the anatomical folder.\n',subj_name{sn(s)})
        end
        
        % MAKE SURE TO FIND ANATOMICAL AFTER THIS STEP AND MOVE TO FOLDER
        % Rename to 'anatomical_raw
    case 'ANAT:reslice_LPI'                  % STEP 1.2: Reslice anatomical image within LPI coordinate systems
        % STUDY 1 ONLY
        sn  = varargin{1}; % subjNum
        % example: sc1_sc2_imana('ANAT:reslice_LPI',1)
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            source  = fullfile(studyDir{1},anatomicalDir,subj_name{sn(s)},['anatomical_raw','.nii']);
            dest    = fullfile(studyDir{1},anatomicalDir,subj_name{sn(s)},['anatomical','.nii']);
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            display 'Manually retrieve the location of the anterior commissure (x,y,z) before continuing'
        end
    case 'ANAT:centre_AC'                    % STEP 1.3: Re-centre AC
        % STUDY 1 ONLY
        % Set origin of anatomical to anterior commissure (must provide
        % coordinates in section (4)).
        % example: sc1_imana('ANAT:centre_AC',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            img    = fullfile(studyDir{1},anatomicalDir,subj_name{sn(s)},['anatomical','.nii']);
            V               = spm_vol(img);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = loc_AC{sn(s)};
            spm_write_vol(V,dat);
            fprintf('Done for %s',subj_name{sn(s)})
        end
    case 'ANAT:segmentation'                 % STEP 1.4: Segmentation + Normalisation
        % STUDY 1 ONLY
        % example: sc1_imana('ANAT:segmentation',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=1:subjs,
            J.channel.vols = {fullfile(studyDir{1},anatomicalDir,subj_name{sn(s)},'anatomical.nii,1')};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 0];
            J.tissue(1).tpm = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
            J.warp.mrf = 1;
            J.warp.cleanup = 1;
            J.warp.reg = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg = 'mni';
            J.warp.fwhm = 0;
            J.warp.samp = 3;
            J.warp.write = [1 1];
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
            fprintf('Check segmentation results for %s\n', subj_name{sn(s)})
        end;
        
    case 'SURF:run_all'                      % STEP 2.1-2.5: enter sn - runs all surface constructions (freesurfer + caret format transformation)
        sn=varargin{1};
        
        subjs=length(sn) ;
        for s=1:subjs,
            sc1_sc2_imana('SURF:recon_all',sn(s));
            sc1_sc2_imana('SURF:xhemireg',sn(s));
            sc1_sc2_imana('SURF:map_ico',sn(s));
            sc1_sc2_imana('SURF:make_caret',sn(s));
        end
    case 'SURF:recon_all'                    % STEP 2.2: Calls recon_all
        % STUDY 1 ONLY
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        % example: sc1_imana('surf_freesurfer',1)
        sn=varargin{1}; % subjNum
        
        for i=sn
            freesurfer_reconall(fullfile(studyDir{1},freesurferDir),subj_name{i},fullfile(anatomicalDir,subj_name{i},['anatomical.nii']));
        end
    case 'SURF:xhemireg'                     % STEP 2.3: cross-register surfaces left / right hem
        % STUDY 1 ONLY
        % surface-based interhemispheric registration
        % example: sc1_imana('surf_xhemireg',1)
        s=varargin{1}; % subjNum
        
        for i=s
            freesurfer_registerXhem({subj_name{i}},fullfile(studyDir{1},freesurferDir),'hemisphere',[1 2]); % For debug... [1 2] orig
        end;
    case 'SURF:map_ico'                      % STEP 2.4: Align to the new atlas surface (map icosahedron)
        % STUDY 1 ONLY
        % Resampels a registered subject surface to a regular isocahedron
        % This allows things to happen in atlas space - each vertex number
        % corresponds exactly to an anatomical location
        % Makes a new folder, called ['x' subj] that contains the remapped subject
        % Uses function mri_surf2surf
        % mri_surf2surf: resamples one cortical surface onto another
        % example: sc1_sc2_imana('SURF:map_ico',1)
        sn=varargin{1}; % subjNum
        
        for i=sn
            freesurfer_mapicosahedron_xhem(subj_name{i},fullfile(studyDir{1},freesurferDir),'smoothing',1,'hemisphere',[1:2]);
        end;
    case 'SURF:make_caret'                   % STEP 2.5: Translate into caret format
        % STUDY 1 ONLY
        % Imports a surface reconstruction from Freesurfer automatically into Caret
        % Makes a new spec file and moves the coord-files to respond to World,
        % rather than to RAS_tkreg coordinates.
        % example: sc1_imana('surf_make_caret',1)
        sn=varargin{1}; % subjNum
        
        for i=sn
            caret_importfreesurfer(['x' subj_name{i}],fullfile(studyDir{1},freesurferDir),fullfile(studyDir{1},caretDir));
        end;
        
    case 'FUNC:func_dicom_import'            % STEP 3.1: Import func dicom
        % converts dicom to nifti files w/ spm_dicom_convert
        % example: sc1_sc2_imana('FUNC:func_dicom_import',1,1,1)
        sn=varargin{1}; % subjNum
        sess=varargin{2}; % sessNum
        study=varargin{1}; % studyNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            dircheck(fullfile(studyDir{study},dicomDir,[subj_name{sn(s)},sprintf('_%d',sess)]));
            cd(fullfile(studyDir{study},dicomDir,[subj_name{sn(s)},sprintf('_%d',sess)]));
            
            if sess==1,
                ss=(sn(s)*2)-1;
            else
                ss=(sn(s)*2);
            end
            
            for i=1:length(fscanNum{study}{ss})
                r=fscanNum{study}{ss}(i);
                DIR=dir(sprintf('%s.%4.4d.*.IMA',DicomName{study}{ss},r));      % Get DICOM FILE NAMES
                Names=vertcat(DIR.name);
                if (~isempty(Names))
                    HDR=spm_dicom_headers(Names,1);                             % Load dicom headers
                    dirname{r}=sprintf('series%2.2d',r);
                    if (~exist(dirname{r}))
                        mkdir(dirname{r});                                      % make dir for series{r} for .nii file output
                    end;
                    dircheck(fullfile(studyDir{study},dicomDir,[subj_name{sn(s)},sprintf('_%d',sess)],dirname{r}));
                    cd(fullfile(studyDir{study},dicomDir,[subj_name{sn(s)},sprintf('_%d',sess)],dirname{r}));
                    
                    spm_dicom_convert(HDR,'all','flat','nii');                  % Convert the data to nifti
                    cd ..
                end;
                display(sprintf('Series %d done \n',fscanNum{study}{ss}(i)))
            end;
            fprintf('Functional runs have been imported. Be sure to copy the unique .nii name for subj files and place into section (4).\n')
        end
    case 'FUNC:make_4dNifti'                 % STEP 3.2: Make 4dNifti
        % merges nifti files for each image into a 4-d nifti (time is 4th
        % dimension) w/ spm_file_merge
        % all imaging data (raw and realigned) for both sc1 and sc2 is
        % saved in sc1
        % example: sc1_sc2_imana('FUNC:make_4dNifti',1,1,1)
        sn=varargin{1}; % subjNum
        sess=varargin{2}; % sessNum
        study=varargin{3}; % studyNum
        
        subjs=length(sn);
        for s=1:subjs,
            dircheck(fullfile(studyDir{1},imagingDirRaw,subj_name{sn(s)}));
            if sess==1,
                ss=(sn(s)*2)-1;
            else
                ss=(sn(s)*2);
            end
            for i=1:length(fscanNum{study}{ss}),  % run number
                outfilename = fullfile(studyDir{study},imagingDirRaw,subj_name{sn(s)},sprintf('run_%2.2d.nii',runs{study}{sess}(i)));
                for j=1:numTRs-numDummys    % doesn't include dummy scans in .nii file
                    P{j}=fullfile(studyDir{study},dicomDir,[subj_name{sn(s)},sprintf('_%d',sess)],sprintf('series%2.2d',fscanNum{study}{ss}(i)),...
                        sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',NiiRawName{study}{ss},fscanNum{study}{ss}(i),j+numDummys,j+numDummys));
                end;
                spm_file_merge(char(P),outfilename);
                fprintf('Run %d done for %s \n',runs{study}{sess}(i),subj_name{sn(s)});
            end;
        end
    case 'FUNC:realign'                      % STEP 3.3: Realign functional images (both sessions)
        % SPM realigns first volume in each run to first volume of first
        % run, and then registers each image in that run to the first
        % volume of that run. Hence also why it's often better to run
        % anatomical before functional scans.
        
        % SPM does this with 4x4 affine transformation matrix in nifti
        % header (see function 'coords'). These matrices convert from voxel
        % space to world space(mm). If the first image has an affine
        % transformation matrix M1, and image two has one (M2), the mapping
        % from 1 to 2 is: M2/M1 (map image 1 to world space-mm - and then
        % mm to voxel space of image 2).
        
        % Registration determines the 6 parameters that determine the rigid
        % body transformation for each image (described above). Reslice
        % conducts these transformations; resampling each image according
        % to the transformation parameters. This is for functional only!
        
        % all of the imaging data for both sc1 and sc2 are saved sc1 folder
        % example: sc1_sc2_imana('FUNC:realign',1,[1:32])
        
        sn=varargin{1}; %subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(studyDir{1},imagingDirRaw,subj_name{sn(s)}));
            spm_jobman % run this step first
            
            data={};
            for i = 1:length(runs),
                for j=1:numTRs-numDummys;
                    data{i}{j,1}=sprintf('run_%2.2d.nii,%d',runs(i),j);
                end;
            end;
            spmj_realign(data);
            fprintf('runsrealigned for %s\n',subj_name{sn(s)});
        end
    case 'FUNC:move_data'                    % STEP 3.4: Move realigned data
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        % example: sc1_sc2_imana('FUNC:move_data',1,[1:32])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            dircheck(fullfile(studyDir{1},imagingDir,subj_name{sn(s)}))
            for r=1:length(runs);
                % move realigned data for each run
                source = fullfile(studyDir{1},imagingDirRaw,subj_name{sn(s)},sprintf('rrun_%2.2d.nii',runs(r)));
                dest = fullfile(studyDir{1},imagingDir,subj_name{sn(s)},sprintf('rrun_%2.2d.nii',runs(r)));
                copyfile(source,dest);
                
                % move realignment parameter files for each run
                source = fullfile(studyDir{1},imagingDirRaw,subj_name{sn(s)},sprintf('rp_run_%2.2d.txt',runs(r)));
                dest = fullfile(studyDir{1},imagingDir,subj_name{sn(s)},sprintf('rp_run_%2.2d.txt',runs(r)));
                copyfile(source,dest);
            end;
            % move mean_epis
            source = fullfile(studyDir{1},imagingDirRaw,subj_name{sn(s)},'meanrun_01.nii');
            dest = fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'meanepi.nii');
            copyfile(source,dest);
            
            fprintf('realigned epi''s moved for %s \n',subj_name{sn(s)})
        end
    case 'FUNC:coreg'                        % STEP 3:5  STUDY 1:Adjust meanepi to anatomical image REQUIRES USER INPUT
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        % - Coregistration is done for study 1 only
        % example: sc1_sc2_imana('FUNC:coreg',1)
        sn=varargin{1};% subjNum
        step=varargin{2}; % 'manual' or 'auto'
        
        cd(fullfile(studyDir{1},anatomicalDir,subj_name{sn}));
        
        switch step,
            case 'manual'
                coregtool;
                keyboard();
            case 'auto'
                % do nothing
        end
        
        % (2) Automatically co-register functional and anatomical images for study 1
        J.ref = {fullfile(studyDir{1},anatomicalDir,subj_name{sn},['anatomical.nii'])};
        J.source = {fullfile(studyDir{1},imagingDir,subj_name{sn},['rmeanepi.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Manually check again
        %         coregtool;
        %         keyboard();
        
        % NOTE:
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
    case 'FUNC:make_samealign'               % STEP 3.6: Align functional images to rmeanepi of study 1
        % Aligns all functional images from both sessions (each study done separately)
        % to rmeanepi of study 1
        % example: sc1_sc2_imana('FUNC:make_samealign',1,[1:32])
        sn=varargin{1}; % subjNum
        runs=varargin{2}; % runNum
        
        prefix='r';
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(studyDir{1},imagingDir,subj_name{sn(s)}));
            
            % Select image for reference
            % Always rmeanepi from study 1
            P{1} = fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmeanepi.nii');
            
            % Select images to be realigned
            Q={};
            for r=1:numel(runs)
                for i=1:numTRs-numDummys;
                    Q{end+1}    = fullfile(studyDir{1},imagingDir,subj_name{sn(s)},...
                        sprintf('%srun_%2.2d.nii,%d',prefix,runs(r),i));
                end;
            end;
            
            % Run spmj_makesamealign_nifti
            spmj_makesamealign_nifti(char(P),char(Q));
            fprintf('functional images realigned for %s',subj_name{sn(s)})
        end
        
        %spmj_checksamealign
    case 'FUNC:make_maskImage'               % STEP 3.7: STUDY 1:Make mask images (noskull and grey_only)
        % Make maskImage meanepi
        % example: sc1_sc2_imana('FUNC:make_maskImage',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            cd(fullfile(studyDir{1},imagingDir,subj_name{sn(s)}));
            
            nam{1}  = fullfile(studyDir{1},imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(studyDir{1},anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{3}  = fullfile(studyDir{1},anatomicalDir, subj_name{sn(s)}, 'c2anatomical.nii');
            nam{4}  = fullfile(studyDir{1},anatomicalDir, subj_name{sn(s)}, 'c3anatomical.nii');
            spm_imcalc(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            
            nam={};
            nam{1}  = fullfile(studyDir{1},imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            spm_imcalc(nam, 'rmask_gray.nii', 'i1>2 & i2>0.4')
            
            nam={};
            nam{1}  = fullfile(studyDir{1},imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{3}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c5anatomical.nii');
            spm_imcalc(nam, 'rmask_grayEyes.nii', 'i1>2400 & i2+i3>0.4')
            
            nam={};
            nam{1}  = fullfile(studyDir{1},imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c5anatomical.nii');
            nam{3}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{4}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c2anatomical.nii');
            nam{5}  = fullfile(studyDir{1}, anatomicalDir, subj_name{sn(s)}, 'c3anatomical.nii');
            spm_imcalc(nam, 'rmask_noskullEyes.nii', 'i1>2000 & (i2+i3+i4+i5)>0.2')
        end
        
    case 'RS:func_dicom_import'              % STEP 4.1: Import rs dicom
        % converts dicom to nifti files w/ spm_dicom_convert
        % example: sc1_sc2_imana('RS:func_dicom_import',2)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            % different saving format for subjs 23-31:
            if sn(s)>22,
                fileEnd='dcm';
            else
                fileEnd='IMA';
            end
            
            cd(fullfile(restDir,dicomDir,[subj_name{sn(s)}]));
            
            for i=1:length(rs_fscanNum{sn(s)})
                r=rs_fscanNum{sn(s)}(i);
                DIR=dir(sprintf('%s.%4.4d.*.%s',rs_DicomName{sn(s)},r,fileEnd));      % Get DICOM FILE NAMES
                Names=vertcat(DIR.name);
                if (~isempty(Names))
                    HDR=spm_dicom_headers(Names,1);                             % Load dicom headers
                    dirname{r}=sprintf('series%2.2d',r);
                    if (~exist(dirname{r}))
                        mkdir(dirname{r});                                      % make dir for series{r} for .nii file output
                    end;
                    dircheck(fullfile(restDir,dicomDir,[subj_name{sn(s)}],dirname{r}));
                    cd(fullfile(restDir,dicomDir,[subj_name{sn(s)}],dirname{r}));
                    
                    spm_dicom_convert(HDR,'all','flat','nii');                  % Convert the data to nifti
                    cd ..
                end;
                display(sprintf('Series %d done \n',rs_fscanNum{sn(s)}(i)))
            end;
            fprintf('Functional runs have been imported. Be sure to copy the unique .nii name for subj files and place into section (4).\n')
        end
    case 'RS:make_4dNifti'                   % STEP 4.2: Make 4dNifti
        % merges nifti files for each image into a 4-d nifti (time is 4th
        % dimension) w/ spm_file_merge
        % all imaging data (raw and realigned) for both sc1 and sc2 is
        % saved in sc1
        % example: sc1_sc2_imana('RS:make_4dNifti',2)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        for s=1:subjs,
            dircheck(fullfile(restDir,imagingDirRaw,subj_name{sn(s)}));
            for i=1:length(rs_fscanNum{sn(s)}),  % run number
                outfilename = fullfile(restDir,imagingDirRaw,subj_name{sn(s)},sprintf('run_%2.2d.nii',i));
                for j=1:numTRs-numDummys    % doesn't include dummy scans in .nii file
                    P{j}=fullfile(restDir,dicomDir,[subj_name{sn(s)}],sprintf('series%2.2d',rs_fscanNum{sn(s)}(i)),...
                        sprintf('f%s-%4.4d-%5.5d-%6.6d-01.nii',rs_NiiRawName{sn(s)},rs_fscanNum{sn(s)}(i),j+numDummys,j+numDummys));
                end;
                spm_file_merge(char(P),outfilename);
                fprintf('Run %d done for %s \n',(i),subj_name{sn(s)});
            end;
        end
    case 'RS:realign'                        % STEP 4.3: Realign functional images for rs
        % SPM realigns first volume in each run to first volume of first
        % run, and then registers each image in that run to the first
        % volume of that run. Hence also why it's often better to run
        % anatomical before functional scans.
        
        % SPM does this with 4x4 affine transformation matrix in nifti
        % header (see function 'coords'). These matrices convert from voxel
        % space to world space(mm). If the first image has an affine
        % transformation matrix M1, and image two has one (M2), the mapping
        % from 1 to 2 is: M2/M1 (map image 1 to world space-mm - and then
        % mm to voxel space of image 2).
        
        % Registration determines the 6 parameters that determine the rigid
        % body transformation for each image (described above). Reslice
        % conducts these transformations; resampling each image according
        % to the transformation parameters. This is for functional only!
        sn=varargin{1}; %subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(restDir,imagingDirRaw,subj_name{sn(s)}));
            spm_jobman % run this step first
            
            data={};
            for i = 1:length(rs_fscanNum{sn(s)}),
                for j=1:numTRs-numDummys;
                    data{i}{j,1}=sprintf('run_%2.2d.nii,%d',i,j);
                end;
            end;
            spmj_realign(data);
            fprintf('rs runs realigned for %s\n',subj_name{sn(s)});
        end
    case 'RS:move_data'                      % STEP 4.4: Move realigned data
        % Moves image data from imaging_dicom_raw into a "working dir":
        % imaging_dicom.
        % example: sc1_sc2_imana('RS:move_data',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            dircheck(fullfile(restDir,imagingDir,subj_name{sn(s)}))
            for r=1:length(rs_fscanNum{sn(s)});
                % move realigned data for each run
                source = fullfile(restDir,imagingDirRaw,subj_name{sn(s)},sprintf('rrun_%2.2d.nii',r));
                dest = fullfile(restDir,imagingDir,subj_name{sn(s)},sprintf('rrun_%2.2d.nii',r));
                copyfile(source,dest);
                
                % move realignment parameter files for each run
                source = fullfile(restDir,imagingDirRaw,subj_name{sn(s)},sprintf('rp_run_%2.2d.txt',r));
                dest = fullfile(restDir,imagingDir,subj_name{sn(s)},sprintf('rp_run_%2.2d.txt',r));
                copyfile(source,dest);
            end;
            % move mean_epis
            source = fullfile(restDir,imagingDirRaw,subj_name{sn(s)},'meanrun_01.nii');
            dest = fullfile(restDir,imagingDir,subj_name{sn(s)},'meanepi.nii');
            copyfile(source,dest);
            
            fprintf('realigned epi''s moved for %s \n',subj_name{sn(s)})
        end
    case 'RS:coreg'                          % STEP 4.5: Coreg to anat
        % (1) Manually seed the functional/anatomical registration
        % - Do "coregtool" on the matlab command window
        % - Select anatomical image and meanepi image to overlay
        % - Manually adjust meanepi image and save result as rmeanepi image
        % - Coregistration is done for study 1 only
        % example: sc1_sc2_imana('RS:coreg',1)
        sn=varargin{1};% subjNum
        step=varargin{2}; % 'manual' or 'auto'
        
        cd(fullfile(baseDir,'sc1',anatomicalDir,subj_name{sn}));
        
        switch step,
            case 'manual'
                coregtool;
                keyboard();
            case 'auto'
                % do nothing
        end
        
        % (2) Automatically co-register functional and anatomical images for study 1
        J.ref = {fullfile(baseDir,'sc1',anatomicalDir,subj_name{sn},['anatomical.nii'])};
        J.source = {fullfile(restDir,imagingDir,subj_name{sn},['rmeanepi.nii'])};
        J.other = {''};
        J.eoptions.cost_fun = 'nmi';
        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Manually check again
        %         coregtool;
        %         keyboard();
        
        % NOTE:
        % Overwrites meanepi, unless you update in step one, which saves it
        % as rmeanepi.
        % Each time you click "update" in coregtool, it saves current
        % alignment by appending the prefix 'r' to the current file
        % So if you continually update rmeanepi, you'll end up with a file
        % called r...rrrmeanepi.
    case 'RS:make_samealign'                 % STEP 4.6: Align rs func images to rmeanepi
        % Aligns all functional images from both sessions (each study done separately)
        % to rmeanepi of study 1
        % example: sc1_sc2_imana('RS:make_samealign',1)
        sn=varargin{1}; % subjNum
        
        prefix='r';
        subjs=length(sn);
        
        for s=1:subjs,
            
            cd(fullfile(restDir,imagingDir,subj_name{sn(s)}));
            
            % Select image for reference
            P{1} = fullfile(restDir,imagingDir,subj_name{sn(s)},'rmeanepi.nii');
            
            % Select images to be realigned
            Q={};
            for r=1:numel(runs)
                for i=1:numTRs-numDummys;
                    Q{end+1}    = fullfile(restDir,imagingDir,subj_name{sn(s)},...
                        sprintf('%srun_%2.2d.nii,%d',prefix,r,i));
                end;
            end;
            
            % Run spmj_makesamealign_nifti
            spmj_makesamealign_nifti(char(P),char(Q));
            fprintf('functional images realigned for %s',subj_name{sn(s)})
        end
    case 'RS:make_maskImage'
        % Make maskImage meanepi
        % example: sc1_sc2_imana('RS:make_maskImage',1)
        sn=varargin{1}; % subjNum
        
        subjs=length(sn);
        
        for s=1:subjs,
            cd(fullfile(restDir,imagingDir,subj_name{sn(s)}));
            
            nam{1}  = fullfile(restDir,imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(baseDir,'sc1',anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{3}  = fullfile(baseDir,'sc1',anatomicalDir, subj_name{sn(s)}, 'c2anatomical.nii');
            nam{4}  = fullfile(baseDir,'sc1',anatomicalDir, subj_name{sn(s)}, 'c3anatomical.nii');
            spm_imcalc(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            
            nam={};
            nam{1}  = fullfile(restDir,imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            spm_imcalc(nam, 'rmask_gray.nii', 'i1>2 & i2>0.4')
            
            nam={};
            nam{1}  = fullfile(restDir,imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{3}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c5anatomical.nii');
            spm_imcalc(nam, 'rmask_grayEyes.nii', 'i1>2400 & i2+i3>0.4')
            
            nam={};
            nam{1}  = fullfile(restDir,imagingDir,subj_name{sn(s)}, 'rmeanepi.nii');
            nam{2}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c5anatomical.nii');
            nam{3}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c1anatomical.nii');
            nam{4}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c2anatomical.nii');
            nam{5}  = fullfile(baseDir,'sc1', anatomicalDir, subj_name{sn(s)}, 'c3anatomical.nii');
            spm_imcalc(nam, 'rmask_noskullEyes.nii', 'i1>2000 & (i2+i3+i4+i5)>0.2')
        end
        
    case 'GLM:study1_glm4'                   % STEP 5.1c:FAST glm w/out hpf (complex:rest as baseline) - model one instruct period
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study1_glm4',[2:22],[1:16])
        sn=varargin{1};
        runs=varargin{2}; % [1:16]
        
        prefix='r';
        announceTime=5;
        glm=4;
        
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==1);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{1},'data', subj_name{sn(s)},['sc1_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{1}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            % annoying but reorder behavioural runs slightly for 2
            % subjects...
            if strcmp(subj_name{sn(s)},'s18')
                runB = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
            elseif strcmp(subj_name{sn(s)},'s21'),
                runB = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
            end
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDir,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames), % loop through trial-types (ex. congruent and incongruent)
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{1},behavDir, subj_name{sn(s)},['sc1_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        
                        % determine trialType (ugly -- but no other way)
                        if isfield(R,'trialType'),
                            tt=(R.trialType==Cc.trialType(c));
                        else
                            tt=Cc.trialType(c);
                        end
                        if strcmp(Cc.taskNames{c},'visualSearch'),
                            tt=(R.setSize==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'nBack') || strcmp(Cc.taskNames{c},'nBackPic')
                            tt=(R.respMade==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'motorImagery') || strcmp(Cc.taskNames{c},'ToM'),
                            tt=1;
                        end
                        
                        onset=[P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                    end
                    
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:study1_glm5'                   % STEP 5.2c:FAST glm w/out hpf (complex:rest as baseline) - model one instruct period - nonAggr!
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study1_glm5',[2:22],[1:16])
        sn=varargin{1};
        runs=varargin{2}; % [1:16];
        
        prefix='r';
        announceTime=5;
        glm=5;
        
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==1);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{1},'data', subj_name{sn(s)},['sc1_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{1}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            % annoying but reorder behavioural runs slightly for 2
            % subjects...
            if strcmp(subj_name{sn(s)},'s18')
                runB = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
            elseif strcmp(subj_name{sn(s)},'s21'),
                runB = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
            end
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDirNA,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])]; % nonAggr
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames), % loop through trial-types (ex. congruent and incongruent)
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{1},behavDir, subj_name{sn(s)},['sc1_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        
                        % determine trialType (ugly -- but no other way)
                        if isfield(R,'trialType'),
                            tt=(R.trialType==Cc.trialType(c));
                        else
                            tt=Cc.trialType(c);
                        end
                        if strcmp(Cc.taskNames{c},'visualSearch'),
                            tt=(R.setSize==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'nBack') || strcmp(Cc.taskNames{c},'nBackPic')
                            tt=(R.respMade==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'motorImagery') || strcmp(Cc.taskNames{c},'ToM'),
                            tt=1;
                        end
                        
                        onset=[P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                    end
                    
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:study1_glm6'                   % STEP 5.3c:FAST glm w/out hpf (complex:rest as baseline) - model one instruct period - Aggr!
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study1_glm6',[2:22],[1:16])
        sn=varargin{1};
        runs=varargin{2}; % [1:16]
        
        prefix='r';
        announceTime=5;
        glm=6;
        
        subjs=length(sn);
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==1);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{1},'data', subj_name{sn(s)},['sc1_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{1}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            % annoying but reorder behavioural runs slightly for 2
            % subjects...
            if strcmp(subj_name{sn(s)},'s18')
                runB = [51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,60];
            elseif strcmp(subj_name{sn(s)},'s21'),
                runB = [51,52,53,54,55,56,57,58,59,60,61,63,64,65,66,62];
            end
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDirA,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])]; % aggr
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames), % loop through trial-types (ex. congruent and incongruent)
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{1},behavDir, subj_name{sn(s)},['sc1_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        
                        % determine trialType (ugly -- but no other way)
                        if isfield(R,'trialType'),
                            tt=(R.trialType==Cc.trialType(c));
                        else
                            tt=Cc.trialType(c);
                        end
                        if strcmp(Cc.taskNames{c},'visualSearch'),
                            tt=(R.setSize==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'nBack') || strcmp(Cc.taskNames{c},'nBackPic')
                            tt=(R.respMade==Cc.trialType(c));
                        elseif strcmp(Cc.taskNames{c},'motorImagery') || strcmp(Cc.taskNames{c},'ToM'),
                            tt=1;
                        end
                        
                        onset=[P.realStartTime(ST)+R.startTimeReal(tt)+announceTime-(J.timing.RT*numDummys)];
                    end
                    
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:study2_glm4'                   % STEP 5.4: FAST glm w/out hpf (complex:rest as baseline) - model one instruct period
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study2_glm4',[2:22],[17:32])
        sn=varargin{1};
        runs=varargin{2}; % [17:32]
        
        subjs=length(sn);
        
        prefix='r';
        announceTime=5;
        glm=4;
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==2);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{2},behavDir,subj_name{sn(s)},['sc2_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{2}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDir,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames),
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{2},behavDir, subj_name{sn(s)},['sc2_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        onset=[P.realStartTime(ST)+R.startTimeReal(R.condition==Cc.trialType(c))+announceTime-(J.timing.RT*numDummys)];
                    end
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:study2_glm5'                   % STEP 5.5: FAST glm w/out hpf (complex:rest as baseline) - model one instruct period - nonAggr!
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study2_glm5',[2:22],[17:32])
        sn=varargin{1};
        runs=varargin{2}; % [17:32]
        
        subjs=length(sn);
        
        prefix='r';
        announceTime=5;
        glm=5;
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==2);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{2},behavDir,subj_name{sn(s)},['sc2_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{2}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDirNA,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames),
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{2},behavDir, subj_name{sn(s)},['sc2_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        onset=[P.realStartTime(ST)+R.startTimeReal(R.condition==Cc.trialType(c))+announceTime-(J.timing.RT*numDummys)];
                    end
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:study2_glm6'                   % STEP 5.6: FAST glm w/out hpf (complex:rest as baseline) - model one instruct period - Aggr!
        % GLM with FAST and no high pass filtering
        % 'spm_get_defaults' code modified to allow for -v7.3 switch (to save
        % >2MB FAST GLM struct)
        % Be aware: this switch (from -v6 to -v7.3) slows down the code!
        % EXAMPLE: sc1_sc2_imana('GLM:study2_glm6',[2:22],[17:32])
        sn=varargin{1};
        runs=varargin{2}; % [17:32]
        
        subjs=length(sn);
        
        prefix='r';
        announceTime=5;
        glm=6;
        
        % load in task information
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds_GLM.txt'));
        Cc=getrow(C,C.StudyNum==2);
        
        for s=1:subjs,
            T=[];
            A = dload(fullfile(studyDir{2},behavDir,subj_name{sn(s)},['sc2_',subj_name{sn(s)},'.dat']));
            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
            
            glmSubjDir =[studyDir{2}, sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];dircheck(glmSubjDir);
            
            J.dir = {glmSubjDir};
            J.timing.units = 'secs';
            J.timing.RT = 1.0;
            J.timing.fmri_t = 16;
            J.timing.fmri_t0 = 1;
            
            for r=1:numel(run) % loop through runs
                P=getrow(A,A.runNum==runB(r));
                for i=1:(numTRs-numDummys)
                    N{i} = [fullfile(studyDir{1},imagingDirA,subj_name{sn(s)},[prefix sprintf('run_%2.2d',runs(r)),'.nii,',num2str(i)])];
                end;
                J.sess(r).scans= N; % number of scans in run
                
                for c=1:length(Cc.condNames),
                    
                    % different onset time for instruct
                    if c==1,
                        onset=[P.realStartTime-J.timing.RT*numDummys];
                    else
                        D=dload(fullfile(studyDir{2},behavDir, subj_name{sn(s)},['sc2_',subj_name{sn(s)},'_',Cc.taskNames{c},'.dat']));
                        R=getrow(D,D.runNum==runB(r)); % functional runs
                        ST = find(strcmp(P.taskName,Cc.taskNames{c}));
                        onset=[P.realStartTime(ST)+R.startTimeReal(R.condition==Cc.trialType(c))+announceTime-(J.timing.RT*numDummys)];
                    end
                    % loop through trial-types (ex. congruent or incongruent)
                    J.sess(r).cond(c).name = Cc.condNames{c};
                    J.sess(r).cond(c).onset = onset; % correct start time for numDummys and announcetime included (not for instruct)
                    J.sess(r).cond(c).duration = Cc.duration(c);  % duration of trials (+ fixation cross) we are modeling
                    J.sess(r).cond(c).tmod = 0;
                    J.sess(r).cond(c).orth = 0;
                    J.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
                    S.SN    = sn(s);
                    S.run   = r;
                    S.task  = Cc.taskNum(c);
                    S.cond  = Cc.condNum(c);
                    S.TN    = {Cc.condNames{c}};
                    S.sess  = sess(r);
                    T=addstruct(T,S);
                end
                J.sess(r).multi = {''};
                J.sess(r).regress = struct('name', {}, 'val', {});
                J.sess(r).multi_reg = {''};
                J.sess(r).hpf = inf;                                        % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
            end
            J.fact = struct('name', {}, 'levels', {});
            J.bases.hrf.derivs = [0 0];
            J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
            J.volt = 1;
            J.global = 'None';
            J.mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_noskull.nii,1')};
            J.mthresh = 0.05;
            J.cvi_mask = {fullfile(studyDir{1},imagingDir,subj_name{sn(s)},'rmask_gray.nii')};
            J.cvi =  'fast';
            
            spm_rwls_run_fmri_spec(J);
            
            save(fullfile(J.dir{1},'SPM_info.mat'),'-struct','T');
            fprintf('glm_%d has been saved for %s \n',glm, subj_name{sn(s)});
        end
    case 'GLM:estimate_glm'                  % STEP 5.7: Enter subjNum & glmNum Takes approx 70 minutes!!
        % example: sc1_sc2_imana('estimate_glm',1,4,1)
        sn=varargin{1};
        glm=varargin{2};
        study=varargin{3};
        
        subjs=length(sn);
        
        for s=1:subjs,
            glmDir =[studyDir{study},sprintf('/GLM_firstlevel_%d/',glm),subj_name{sn(s)}];
            load(fullfile(glmDir,'SPM.mat'));
            SPM.swd=glmDir;
            spm_rwls_spm(SPM);
        end
    case 'GLM:contrast'                      % STEP 5.8: Define linear contrasts
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM)
        % example: sc1_sc2_imana('contrast',1,4,1)
        s=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        study=varargin{3}; % studyNum
        
        glmSubjDir =[studyDir{study},sprintf('/GLM_firstlevel_%d/%s',glm,subj_name{s})];
        cd(glmSubjDir);
        load SPM;
        SPM=rmfield(SPM,'xCon');
        T=load('SPM_info.mat');
        nrun=numel(SPM.nscan);
        
        % t contrast for tasks vs. rest
        condNum=unique(T.cond);
        idx=1;
        for tt=1:length(condNum), % 0 is "instruct" regressor
            con=zeros(1,size(SPM.xX.X,2));
            con(:,T.cond==condNum(tt))=1; % contrast against rest
            con=con/abs(sum(con));
            name=sprintf('%s-rest',char(unique(T.TN(T.cond==condNum(tt)))));
            SPM.xCon(idx)=spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
            idx=idx+1;
        end
        
        SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
        save('SPM.mat','SPM','-v7.3');
        SPM=rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
        save(fullfile(glmSubjDir,'SPM_light.mat'),'SPM');
        
        % rename contrast images and spmT images
        conName = {'con','spmT'};
        for i=1:length(SPM.xCon),
            for n=1:numel(conName),
                oldName{i} = fullfile(glmSubjDir,sprintf('%s_%2.4d.nii',conName{n},i));
                newName{i} = fullfile(glmSubjDir,sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                movefile(oldName{i},newName{i});
            end
        end
        
    case 'SUIT:run_all'                      % STEP 9.1-8.5
        % example: sc1_sc2_imana('SUIT:run_all',2,4)
        sn = varargin{1}; % subjNum
        glm = varargin{2}; % glmNum
        spm fmri
        % STUDY 1 ONLY EXCEPT FOR 'SUIT:reslice'
        
        for s=sn,
            sc1_sc2_imana('SUIT:isolate_segment',s)
            sc1_sc2_imana('SUIT:make_mask_cortex',s)
            sc1_sc2_imana('SUIT:corr_cereb_cortex_mask',s)
            sc1_sc2_imana('SUIT:normalise_dartel',s,'grey');
            %             sc1_sc2_imana('suit_normalise_dentate',s,'grey');
            sc1_sc2_imana('SUIT:make_mask',s,glm,'grey');
            %             sc1_sc2_imana('SUIT:reslice',s,glm,'contrast','cereb_prob_corr_grey');
            sc1_sc2_imana('SUIT:reslice',s,1,glm,'ResMS','cereb_prob_corr_grey');
            sc1_sc2_imana('SUIT:reslice',s,2,glm,'ResMS','cereb_prob_corr_grey');
            sc1_sc2_imana('SUIT:reslice',s,1,glm,'betas','cereb_prob_corr_grey');
            sc1_sc2_imana('SUIT:reslice',s,2,glm,'betas','cereb_prob_corr_grey');
            fprintf('suit data processed for %s',subj_name{s})
        end
    case 'SUIT:isolate_segment'              % STEP 9.2:Segment cerebellum into grey and white matter
        % STUDY 1 ONLY
        sn=varargin{1};
        %         spm fmri
        for s=sn,
            suitSubjDir = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{s});dircheck(suitSubjDir);
            source=fullfile(studyDir{1},anatomicalDir,subj_name{s},'anatomical.nii');
            dest=fullfile(suitSubjDir,'anatomical.nii');
            copyfile(source,dest);
            cd(fullfile(suitSubjDir));
            suit_isolate_seg({fullfile(suitSubjDir,'anatomical.nii')},'keeptempfiles',1);
        end
    case 'SUIT:make_mask_cortex'             % STEP 9.3:
        % STUDY 1
        sn=varargin{1};
        
        subjs=length(sn);
        for s=1:subjs,
            glmSubjDir =[studyDir{1} '/GLM_firstlevel_4/' subj_name{sn(s)}];
            
            for h=regSide,
                C=caret_load(fullfile(studyDir{1},caretDir,atlasname,hemName{h},[hem{h} '.cerebral_cortex.paint'])); % freesurfer
                caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                file=fullfile(glmSubjDir,'mask.nii');
                r{h}.type='surf_nodes';
                r{h}.location=find(C.data(:,1)~=6); % we don't want the medial wall
                r{h}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                r{h}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                r{h}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                r{h}.linedef=[5,0,1];
                r{h}.image=file;
                r{h}.file=file;
                r{h}.name=[hem{h}];
            end;
            r=region_calcregions(r);
            r{1}.location=vertcat(r{1}.location,r{2}.location);
            r{1}.linvoxidxs=vertcat(r{1}.linvoxidxs,r{2}.linvoxidxs);
            r{1}.data=vertcat(r{1}.data,r{2}.data);
            r{1}.weight=vertcat(r{1}.weight,r{2}.weight);
            r{1}.name='cortical_mask_grey';
            R{1}=r{1};
            dircheck(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)}));
            cd(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)}))
            region_saveasimg(R{1},R{1}.file);
        end
    case 'SUIT:corr_cereb_cortex_mask'       % STEP 9.4:
        sn=varargin{1};
        % STUDY 1
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            cortexGrey= fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey.nii'); % cerebellar mask grey (corrected)
            cerebGrey = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'c1anatomical.nii'); % was 'cereb_prob_corr_grey.nii'
            bufferVox = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'buffer_voxels.nii');
            
            % isolate overlapping voxels
            spm_imcalc({cortexGrey,cerebGrey},bufferVox,'(i1.*i2)')
            
            % mask buffer
            spm_imcalc({bufferVox},bufferVox,'i1>0')
            
            cerebGrey2 = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'cereb_prob_corr_grey.nii');
            cortexGrey2= fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
            
            % remove buffer from cerebellum
            spm_imcalc({cerebGrey,bufferVox},cerebGrey2,'i1-i2')
            
            % remove buffer from cortex
            spm_imcalc({cortexGrey,bufferVox},cortexGrey2,'i1-i2')
        end
    case 'SUIT:normalise_dartel'             % STEP 9.5: Normalise the cerebellum into the SUIT template.
        % STUDY 1
        % Normalise an individual cerebellum into the SUIT atlas template
        % Dartel normalises the tissue segmentation maps produced by suit_isolate
        % to the SUIT template
        % !! Make sure that you're choosing the correct isolation mask
        % (corr OR corr1 OR corr2 etc)!!
        % if you are running multiple subjs - change to 'job.subjND(s)."'
        % example: sc1_sc2_imana('SUIT:normalise_dartel',1,'grey')
        sn=varargin{1}; %subjNum
        type=varargin{2}; % 'grey' or 'whole' cerebellar mask
        
        subjs=length(sn);
        for s=1:subjs,
            cd(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)}));
            job.subjND.gray      = {'a_c_anatomical_seg1.nii'};
            job.subjND.white     = {'a_c_anatomical_seg2.nii'};
            switch type,
                case 'grey'
                    job.subjND.isolation= {'cereb_prob_corr_grey.nii'};
                case 'whole'
                    job.subjND.isolation= {'cereb_prob_corr.nii'};
            end
            suit_normalize_dartel(job);
        end
        
        % 'spm_dartel_warp' code was changed to look in the working
        % directory for 'u_a_anatomical_segment1.nii' file - previously it
        % was giving a 'file2mat' error because it mistakenly believed that
        % this file had been created
    case 'SUIT:normalise_dentate'            % STEP 9.6: Uses an ROI from the dentate nucleus to improve the overlap of the DCN
        % STUDY 1
        sn=varargin{1}; %subjNum
        type=varargin{2}; % 'grey' or 'whole'
        % example: 'sc1_sc2_imana('SUIT:normalise_dentate',2,'grey'
        
        cd(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn}));
        job.subjND.gray       = {'c_anatomical_seg1.nii'};
        job.subjND.white      = {'c_anatomical_seg2.nii'};
        job.subjND.dentateROI = {fullfile(studyDir{1},suitDir,'glm4',subj_name{sn},'dentate_mask.nii')};
        switch type,
            case 'grey'
                job.subjND.isolation  = {'cereb_prob_corr_grey.nii'};
            case 'whole'
                job.subjND.isolation  = {'cereb_prob_corr.nii'};
        end
        
        suit_normalize_dentate(job);
    case 'SUIT:make_mask'                    % STEP 9.7: Make cerebellar mask using SUIT
        % STUDY 1 ONLY
        % example: sc1_sc2_imana('SUIT:make_mask',1,4,'grey')
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        type=varargin{3}; % 'grey' or 'whole'
        
        subjs=length(sn);
        
        for s=1:subjs,
            glmSubjDir = fullfile(studyDir{1},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
            mask       = fullfile(glmSubjDir,'mask.nii'); % mask for functional image
            switch type
                case 'grey'
                    suit  = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'cereb_prob_corr_grey.nii'); % cerebellar mask grey (corrected)
                    omask = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'maskbrainSUITGrey.nii'); % output mask image - grey matter
                case 'whole'
                    suit  = fullfile(studyDir{1},suitDir,'anatomicals', subj_name{sn(s)},'cereb_prob_corr.nii'); % cerebellar mask (corrected)
                    omask = fullfile(studyDir{1},suitDir, sprintf('glm%d',glm),subj_name{sn(s)},'maskbrainSUIT.nii'); % output mask image
            end
            dircheck(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)}));
            cd(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)}));
            spm_imcalc({mask,suit},omask,'i1>0 & i2>0.7',{});
        end
    case 'SUIT:reslice'                      % STEP 9.8: Reslice the contrast images from first-level GLM
        % Reslices the functional data (betas, contrast images or ResMS)
        % from the first-level GLM using deformation from
        % 'suit_normalise_dartel'.
        % example: sc1_sc2_imana('SUIT:reslice',1,1,4,'betas','cereb_prob_corr_grey')
        % make sure that you reslice into 2mm^3 resolution
        sn=varargin{1}; % subjNum
        study=varargin{2}; % studyNum
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'betas' or 'contrast' or 'ResMS' or 'cerebellarGrey'
        mask=varargin{5}; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        
        subjs=length(sn);
        
        for s=1:subjs,
            switch type
                case 'betas'
                    glmSubjDir = fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='beta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'contrast'
                    glmSubjDir = fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='con';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'ResMS'
                    glmSubjDir = fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)});
                    outDir=fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    images='ResMS';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'cerebellarGrey'
                    source=dir(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'c1anatomical.nii')); % image to be resliced
                    cd(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)}));
            end
            job.subj.affineTr = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},sprintf('%s.nii',mask))};
            job.vox           = [2 2 2];
            suit_reslice_dartel(job);
            if ~strcmp(type,'cerebellarGrey'),
                source=fullfile(glmSubjDir,'*wd*');
                dircheck(fullfile(outDir));
                destination=fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)});
                movefile(source,destination);
            end
            fprintf('%s have been resliced into suit space for %s \n\n',type,glm,subj_name{sn(s)})
        end
        
    case 'HOUSEKEEPING:renameSPM'            % HOUSEKEEPING: subjs 23-31 were processed by chernandez so SPM dir is hardcoded as such
        sn=varargin{1};
        study=varargin{2};
        
        subjs=length(sn);
        
        newDir=fullfile(studyDir{1},imagingDir);
        
        for s=1:subjs,
            load(fullfile(studyDir{study},'GLM_firstlevel_4',subj_name{sn(s)},'SPM.mat'));
            newSubjDir=fullfile(newDir,subj_name{sn(s)}); dircheck(newSubjDir);
            SPM=spmj_move_rawdata(SPM,newSubjDir);
            SPM.swd=fullfile(studyDir{study},'GLM_firstlevel_4',subj_name{sn(s)});
            save(fullfile(studyDir{study},'GLM_firstlevel_4',subj_name{sn(s)},'SPM.mat'),'-v7.3','SPM');
        end
    case 'HOUSEKEEPING:renameROI'            % HOUSEKEEPING: subjs 23-31 were processed by chernandez so ROI file is hardcoded as such
        sn=varargin{1};
        type=varargin{2};
        
        subjs=length(sn);
        
        for s=1:subjs,
            load(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            R{1}.file=[];
        end
        
    case 'ROI:run_all'
        % STEP 9.1-9.6
        % example:
        % 'sc1_sc2_imana('ROI:run_all',[2:22],4,{'cerebellum_grey'},1)
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        type=varargin{3}; %'cortical_lobes','whole_brain','yeo','desikan','cerebellum_grey'
        study=varargin{4}; % studyNum
        
        subjs=length(sn);
        for s=1:subjs,
            for t=1:numel(type),
                %                 sc1_sc2_imana('ROI:define',sn(s),type{t})
                sc1_sc2_imana('ROI:betas',sn(s),study,glm,type{t})
                sc1_sc2_imana('ROI:stats',sn(s),study,glm,1,type{t}) % remove mean
                %             sc1_SC2_imana('ROI:RDM_stability',s,glm,type)
            end
        end
    case 'ROI:define'                        % STEP 10.1: Enter subjNum and glmNum. Defines ROIs that are referenced in (2) at start of sc1_sc2_imana
        % STUDY 1 ONLY
        
        % Run FREESURFER before this step!
        sn=varargin{1}; % subjNum
        type=varargin{2}; % 'cortical_lobes','yeo','desikan','cerebellum','cerebellum_grey'
        
        subjs=length(sn);
        idx=0;
        
        atlasName=fullfile(studyDir{1},caretDir,atlasname);
        
        for s=1:subjs,
            
            switch type
                case 'cortical_lobes'
                    for h=regSide,
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.cerebral_cortex.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                        if h==2,
                            idx=numel(C.paintnames)-1; % we don't want medial wall
                        end
                        % Get cortical regions
                        for i=1:numel(C.paintnames)-1, % we don't want the medial wall
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==i);
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i}];
                        end;
                    end;
                    r=region_calcregions(r);
                    for i=1:numel(C.paintnames)-1,
                        r{i}.location=vertcat(r{i}.location,r{i+idx}.location);
                        r{i}.linvoxidxs=vertcat(r{i}.linvoxidxs,r{i+idx}.linvoxidxs);
                        r{i}.data=vertcat(r{i}.data,r{i+idx}.data);
                        r{i}.weight=vertcat(r{i}.weight,r{i+idx}.weight);
                        R{i}=r{i};
                    end
                case 'yeo',
                    for h=regSide, % loop over hemispheres
                        % define Yeo networks as regions
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.Yeo17.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==i);
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i}];
                        end;
                    end
                    r=region_calcregions(r);
                    for i=1:numel(C.paintnames),
                        r{i}.location=vertcat(r{i}.location,r{i+idx}.location);
                        r{i}.linvoxidxs=vertcat(r{i}.linvoxidxs,r{i+idx}.linvoxidxs);
                        r{i}.data=vertcat(r{i}.data,r{i+idx}.data);
                        r{i}.weight=vertcat(r{i}.weight,r{i+idx}.weight);
                        R{i}=r{i};
                    end
                case 'yeo_hem'
                    for h=regSide,
                        % define desikan as regions
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.Yeo17.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==i);
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i} '_' hem{h}];
                        end;
                    end
                    R=region_calcregions(r);
                case 'desikan'
                    for h=regSide,
                        % define desikan as regions
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.desikan.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==i);
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i}];
                        end;
                    end
                    r=region_calcregions(r);
                    for i=1:numel(C.paintnames),
                        r{i}.location=vertcat(r{i}.location,r{i+idx}.location);
                        r{i}.linvoxidxs=vertcat(r{i}.linvoxidxs,r{i+idx}.linvoxidxs);
                        r{i}.data=vertcat(r{i}.data,r{i+idx}.data);
                        r{i}.weight=vertcat(r{i}.weight,r{i+idx}.weight);
                        R{i}=r{i};
                    end
                case 'desikan_hem'
                    for h=regSide,
                        % define desikan as regions
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.desikan.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==i);
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i} '_' hem{h}];
                        end;
                    end
                    R=region_calcregions(r);
                case 'cerebellum_grey'
                    % Get cerebellum
                    file = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'maskbrainSUITGrey.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['cerebellum_grey'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
                case '162_tessellation'
                    caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                    file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii');
                    for h=regSide,
                        C=caret_load(fullfile(atlasName,hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
                        tesselNames=unique(C.data);
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==tesselNames(i));
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i}];
                        end;
                        r{149}.name='medial_wall'; % tessel 149 is the medial wall
                    end
                    r=region_calcregions(r);
                    for i=1:numel(C.paintnames),
                        r{i}.location=vertcat(r{i}.location,r{i+idx}.location);
                        r{i}.linvoxidxs=vertcat(r{i}.linvoxidxs,r{i+idx}.linvoxidxs);
                        r{i}.data=vertcat(r{i}.data,r{i+idx}.data);
                        r{i}.weight=vertcat(r{i}.weight,r{i+idx}.weight);
                        R{i}=r{i};
                    end
                case '162_tessellation_hem'
                    caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                    file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii'); % used to be 'mask'
                    for h=regSide,
                        C=caret_load(fullfile(atlasName,hemName{h},sprintf('%s.tessel162.paint',hem{h}))); % freesurfer
                        tesselNames=unique(C.data);
                        if h==2,
                            idx=numel(C.paintnames);
                        end
                        for i=1:numel(C.paintnames),
                            r{i+idx}.type='surf_nodes';
                            r{i+idx}.location=find(C.data(:,1)==tesselNames(i));
                            r{i+idx}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                            r{i+idx}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                            r{i+idx}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                            r{i+idx}.linedef=[5,0,1];
                            r{i+idx}.image=file;
                            r{i+idx}.file=file;
                            r{i+idx}.name=[C.paintnames{i} '_' hem{h}];
                        end;
                        r{149}.name='medial_wall_lh'; % tessel 149 is the medial wall
                        r{149+158}.name='medial_wall_rh'; % tessel 149+158 is the other medial wall
                    end
                    R=region_calcregions(r);
                case 'cortex_grey'
                    for h=regSide,
                        C=caret_load(fullfile(atlasName,hemName{h},[hem{h} '.cerebral_cortex.paint'])); % freesurfer
                        caretSubjDir=fullfile(studyDir{1},caretDir,[atlasA subj_name{sn(s)}]);
                        file=fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},'cortical_mask_grey_corr.nii'); % used to be just 'mask'
                        r{h}.type='surf_nodes';
                        r{h}.location=find(C.data(:,1)~=6); % we don't want the medial wall
                        r{h}.white=fullfile(caretSubjDir,hemName{h},[hem{h} '.WHITE.coord']);
                        r{h}.pial=fullfile(caretSubjDir,hemName{h},[hem{h} '.PIAL.coord']);
                        r{h}.topo=fullfile(caretSubjDir,hemName{h},[hem{h} '.CLOSED.topo']);
                        r{h}.linedef=[5,0,1];
                        r{h}.image=file;
                        r{h}.file=file;
                        r{h}.name=[hem{h}];
                    end;
                    r=region_calcregions(r);
                    R=r;
                case 'dentate'
                    % Get dentate nucleus
                    file = fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'dentate_mask.nii');
                    R{1}.type = 'roi_image';
                    R{1}.file= file;
                    R{1}.name = ['dentate'];
                    R{1}.value = 1;
                    R=region_calcregions(R);
            end
            dircheck(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)}));
            %             dircheck(fullfile(regDir,'glm4',subj_name{sn(s)},sprintf('%s_masks',type)));
            %             cd(fullfile(regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('%s_masks',type)));
            %             for ii=1:length(R),
            %                 region_saveasimg(R{ii},R{ii}.file);
            %             end
            save(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type)),'R');
            fprintf('ROIs have been defined for %s for %s \n',type,subj_name{sn(s)})
        end
    case 'ROI:betas'                         % STEP 10.2: Extract betas and prewhiten (apply none, uni and multi noise normalisation)
        % Betas are not multivariately noise-normalised in this version of
        % the code. See 'sc1_imana_backUp' code for this option
        % Betas are computed for all regions (both hemispheres for cortical ROIs).
        % Regions are defined in section (2)
        sn=varargin{1}; % subjNum
        study=varargin{2}; % studyNum
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'cortical_lobes','whole_brain','yeo','desikan','cerebellum','yeo_cerebellum'
        
        B = [];
        glmDir =[studyDir{study} sprintf('/GLM_firstlevel_%d',glm)];dircheck(glmDir);
        subjs=length(sn);
        
        for s=1:subjs,
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            load(fullfile(glmDirSubj,'SPM.mat'));
            
            % load data
            load(fullfile(studyDir{1},regDir,'data',subj_name{sn(s)},sprintf('regions_%s.mat',type))); % 'regions' are defined in 'ROI_define'
            
            % Get the raw data files
            V=SPM.xY.VY;
            
            Y = region_getdata(V,R);  % Data is N x P
            
            for r=1:numel(R), % R is the output 'regions' structure from 'ROI_define'
                % Get betas (univariately prewhitened)
                [resMS,beta_hat]=rsa.spm.noiseNormalizeBeta(Y{r},SPM,'normmode','none');
                
                % option to include empty betasUW matrix (no betas for some
                % ROIs for some subjects)
                if isempty(beta_hat),
                    betasUW=zeros(size(beta_hat,1),1);
                else
                    betasUW = bsxfun(@rdivide,beta_hat,sqrt(resMS)); % univariate noise normalisation
                end
                
                B{r}.betasNW = beta_hat; % no noise normalisation
                B{r}.betasUW = betasUW; % univariate noise normalisation
                B{r}.region  = r;
                B{r}.regName = R{r}.name;
                B{r}.SN      = sn(s);
            end
            
            % Save output for each subject
            outfile = sprintf('betas_%s.mat',type);
            dircheck(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)}));
            save(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},outfile),'B');
            fprintf('betas computed and saved for %s (%s) for study%d %s \n',subj_name{sn(s)},sprintf('glm%d',glm),study,type);
        end
    case 'ROI:stats'                         % STEP 10.3: Calculate G/second-moment matrix,distance estimates (Mahalanobis),pattern consistencies for all regions
        sn=varargin{1}; % subjNum
        study=varargin{2}; % studyNum
        glm =varargin{3}; % glmNum
        remove_mean=varargin{4}; % better to remove mean for accurate pattern consistencies (input:1)
        type=varargin{5}; % 'all','yeo','desikan','cerebellum' 'cerebellum_grey'
        
        glmDir =[studyDir{study},sprintf('/GLM_firstlevel_%d',glm)];dircheck(glmDir);
        
        subjs=length(sn);
        
        for s=1:subjs,
            Ts=[];  % Seperated by tasks
            To=[];  % Overall
            
            glmDirSubj=fullfile(glmDir, subj_name{sn(s)});
            D=load(fullfile(glmDirSubj,'SPM_info.mat'));
            N=length(D.SN);
            
            % load activity patterns
            load(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('betas_%s.mat',type)));
            
            for r=1:numel(B), % loop over regions
                
                for w=1:numel(PW), % loop over different noise normalisation approaches
                    
                    betas = B{r}.(sprintf('%s',PW{w})); % get correct betas;
                    
                    % Pattern consistency and Split-Half Reliability
                    
                    % Split half reliability
                    split = [1;2];
                    
                    % loop over datasets (odd and even runs for distance measures and
                    % session 1 and session 2 for pattern consistency)
                    for i=1:numel(split),
                        % Get distances
                        hPart =D.run .* (mod(D.run+i-1,2)); % shouldn't trialType be specified here?? (D.run.*D.type==2.*(mod(Drun+i-1,2));
                        So{r}.(sprintf('RDMh_%d',i))  = rsa.distanceLDC(betas,hPart,D.cond); % LDC is the crossval Mahalanobis distance (+ multiv noise normalisation)
                        
                        % Get pattern consistencies within sessions
                        partition = (D.run.*(D.cond~=0)).*(D.sess==i); % partition is usually defined as runs
                        if remove_mean==0
                            So{r}.(sprintf('R2_Sess%d',split(i))) = rsa_patternConsistency(betas,partition,D.cond,'removeMean',0);
                        else
                            So{r}.(sprintf('R2_Sess%d',split(i))) = rsa_patternConsistency(betas,partition,D.cond); % better to remove mean
                        end
                    end
                    
                    % Get pattern consistencies across sessions
                    partition = D.run.*(D.cond~=0);
                    if remove_mean==0,
                        So{r}.R2_overall = rsa_patternConsistency(betas,partition,D.cond,'removeMean',0);
                    else
                        So{r}.R2_overall = rsa_patternConsistency(betas,partition,D.cond);
                    end
                    
                    % Estimate crossval second moment matrix
                    b{1}=betas(1:N,:);
                    condVec{1}=D.cond;
                    partVec{1}=partition;
                    
                    [G,Sig] = pcm_estGCrossval(b{1},partVec{1},condVec{1}); % all betas given here - D.tt indicates betas of interest
                    Ss{r}.IPM = rsa_vectorizeIPM(G); % vectorise second-moment matrix for later (multi-dimensional scaling)
                    Ss{r}.Sig = rsa_vectorizeIPM(Sig); % sigma is optional
                    
                    So{r}.SN         = sn(s); % subjNum
                    So{r}.region     = r; % regionNum
                    So{r}.regName    = {B{r}.regName}; % regionName
                    So{r}.GLM        = glm; % glmNum
                    So{r}.numvox     = size(betas,2); % number of voxels
                    So{r}.method     = PW(w); % noise normalisation approach (none or uni)
                    So{r}.method_num = w;
                    Ss{r}.SN         = sn(s);
                    Ss{r}.region     = r;
                    Ss{r}.regName    = {B{r}.regName};
                    Ss{r}.GLM        = glm;
                    Ss{r}.method     = PW(w);
                    Ss{r}.method_num = w;
                    Ss{r}.numvox     = size(betas,2);
                    Ts               = addstruct(Ts,Ss{r});
                    To               = addstruct(To,So{r});
                end
            end
            % Save output
            save(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('Toverall_%s.mat',type)),'To'); % used for correlation statistics in 'ROI_RDM_stability'
            save(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('Ttasks_%s.mat',type)),'Ts'); % used to visualise RDMs and MDS
            fprintf('GLM %d stats are done for study%d %s (%s). \n',glm,study,subj_name{sn(s)},type)
        end
    case 'ROI:RDM_stability'                 % STEP 10.4: Pearson's corr for split-half reliability of LDC distances
        % example: sc1_imana('ROI_RDM_stability',1,1)
        sn=varargin{1}; % subjNum
        glm=varargin{2}; % glmNum
        type=varargin{3}; % 'cortical_lobes','whole_brain','yeo','desikan','cerebellum'
        
        split = [1;2]; % runs - odd and even
        
        T = [];
        subjs=length(sn);
        for s=1:subjs,
            T=[];
            
            switch glm
                case 4
                    glmType=1;
                case 13
                    glmType=2;
            end
            
            % load overall statistics
            load(fullfile(regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('Toverall_%s.mat',type))); % load in output from 'ROI_stats'
            
            numReg=length(unique(To.region));
            
            for r=1:numReg, % loop over regions
                for w=1:numel(PW), % loop over methods
                    P = getrow(To,To.method_num==w & To.region==r);
                    for i=1:numel(split), % loop over splits (odd and even runs)
                        X(i,:)  = P.(sprintf('RDMh_%d',i)); % distance estimates for splits
                    end
                    [Cs,~] = corr(X','type','pearson');
                    C.region    = r;
                    C.regName   = P.regName;
                    C.corrName  = 'pearson';
                    C.method    = PW(w);
                    C.method_num= w;
                    C.SN        = sn(s);
                    C.corr      = Cs(2);
                    C.GLM       = glm;
                    C.glm_type  = glm_type(glmType);
                    
                    T = addstruct(T,C);
                    C=[];
                end
            end
            
            % Save output
            save(fullfile(regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('correlations_%s.mat',type)),'T');
            fprintf('Correlation stats are done for %s (%s) for %s \n',subj_name{sn(s)}, sprintf('glm%d',glm),type)
        end
        
    case 'PREP:run_indiv'                    % STEP 11.1:
        sn=varargin{1}; % subjNum. use 'returnSubjs' for the 'PREP:voxels_all' and 'PREP:avrg_betas' cases
        study=varargin{2};
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'cereb' or 'cortex'
        
        switch type,
            case 'cereb'
                sc1_sc2_imana('PREP:cereb:suit_betas',sn,study,glm)
                sc1_sc2_imana('PREP:cereb:indiv_voxels',sn,study,glm,'grey_nan')
            case 'cortex'
                sc1_sc2_imana('PREP:cortex:surface_betas',sn,study,glm,'beta')
                sc1_sc2_imana('PREP:cortex:surface_betas',sn,study,glm,'ResMS')
                sc1_sc2_imana('PREP:cortex:indiv_vertices',sn,study,glm)
        end
    case 'PREP:run_allSubj'                  % STEP 11.2:
        study=varargin{1};
        glm=varargin{2}; % glmNum
        type=varargin{3}; % 'cereb' or 'cortex'
        
        switch type,
            case 'cereb'
                sc1_sc2_imana('PREP:betas_all',returnSubjs,glm,'cerebellum',[1])
                sc1_sc2_imana('PREP:avrg_betas',returnSubjs,study,glm,'cereb',[1])
            case 'cortex'
                sc1_sc2_imana('PREP:betas_all',returnSubjs,glm,'cortex',[1:2])
                sc1_sc2_imana('PREP:avrg_betas',returnSubjs,study,glm,'cortex',[1:2])
        end
    case 'PREP:avrgMask_cereb'               % STEP 11.3:
        sn=varargin{1};
        
        subjs=length(sn);
        for s=1:subjs,
            nam{s}=fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'wdc1anatomical.nii');
        end
        opt.dmtx = 1;
        cd(fullfile(studyDir{1},suitDir,'anatomicals'));
        spm_imcalc(nam,'cerebellarGreySUIT.nii','mean(X)',opt);
        
        fprintf('averaged cerebellar grey mask in SUIT space has been computed \n')
    case 'PREP:cereb:suit_betas'             % STEP 11.4: Reslice univar pre-whitened betas into suit space
        sn=varargin{1};
        study=varargin{2};
        glm=varargin{3};
        
        subjs=length(sn);
        
        for s=1:subjs,
            
            %             load betas (grey) from cerebellum
            load(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},'betas_cerebellum_grey.mat'));
            
            %             load cerebellar mask in individual func space
            Vi=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'maskbrainSUITGrey.nii'));
            X=spm_read_vols(Vi);
            indx=find(X>0);
            
            %             make volume
            for b=1:size(B{1}.betasUW,1),
                Yy=zeros(1,Vi.dim(1)*Vi.dim(2)*Vi.dim(3));
                Yy(1,indx)=B{1}.betasUW(b,:);
                Yy=reshape(Yy,[Vi.dim(1),Vi.dim(2),Vi.dim(3)]);
                Yy(Yy==0)=NaN;
                Vi.fname=fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)},sprintf('temp_cereb_beta_%2.4d.nii',b));
                spm_write_vol(Vi,Yy);
                clear Yy
                filenames{b}=Vi.fname;
                fprintf('beta %d done \n',b)
            end
            %             reslice univar prewhitened betas into suit space
            job.subj.affineTr = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = filenames';
            job.subj.mask     = {fullfile(studyDir{1},suitDir,'anatomicals',subj_name{sn(s)},'cereb_prob_corr_grey.nii')};
            job.vox           = [2 2 2];
            job.outFile       = 'mat';
            D=suit_reslice_dartel(job);
            %             delete temporary files
            deleteFiles=dir(fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)},'*temp*'));
            for b=1:length(deleteFiles),
                delete(char(fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm),subj_name{sn(s)},deleteFiles(b).name)));
            end
            save(fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)},'wdBetas_UW.mat'),'D');
            fprintf('UW betas resliced into suit space for %s \n',subj_name{sn(s)});
        end
    case 'PREP:cortex:surface_betas'         % STEP 11.5: Map betas and ResMS (.nii) onto surface (.metric)
        % Run FREESURFER before this step!
        % map volume images to metric file and save them in individual
        % surface folder
        % example: sc1_sc2_imana('map_con_surf',2,1,4,'betas')
        sn   = varargin{1}; % subjNum
        study=varargin{2};
        glm  = varargin{3}; % glmNum
        contrast = varargin{4}; % 'beta'or 'ResMS'
        
        glmDir =[studyDir{study} sprintf('/GLM_firstlevel_%d',glm)];dircheck(glmDir);
        
        subjs=length(sn);
        
        vararginoptions({varargin{5:end}},{'atlas','regSide'});
        
        for s=1:subjs,
            glmSubjDir=fullfile(glmDir, subj_name{sn(s)});
            for h=regSide,
                caretSDir = fullfile(studyDir{study}, caretDir,[atlasA,subj_name{sn(s)}],hemName{h}); dircheck(caretSDir);
                white=fullfile(studyDir{1},caretDir,[atlasA,subj_name{sn(s)}],hemName{h},[hem{h} '.WHITE.coord']);
                pial=fullfile(studyDir{1},caretDir,[atlasA,subj_name{sn(s)}],hemName{h},[hem{h} '.PIAL.coord']);
                
                C1=caret_load(white);
                C2=caret_load(pial);
                fileList = [];
                outfile  = [];
                
                filenames = dir(fullfile(glmSubjDir,sprintf('*%s*',contrast)));
                outfile = sprintf('%s_glm%d_%s_cortex_%s.metric',subj_name{sn(s)},glm,contrast,hem{h});
                
                for t = 1:length(filenames),
                    fileList{t} = {filenames(t).name};
                end
                
                for f=1:length(fileList),
                    images(f)=spm_vol(fullfile(glmSubjDir,fileList{f}));
                end;
                M=caret_vol2surf_own(C1.data,C2.data,images);
                caret_save(fullfile(caretSDir,outfile),M);
                
                fprintf('%s map to surface for %s:%s \n',contrast,subj_name{sn(s)},hemName{h});
            end;
        end
    case 'PREP:cereb:indiv_voxels'           % STEP 11.6: Get UW cerebellar data (voxels)
        sn=varargin{1};
        study=varargin{2};
        glm=varargin{3};
        data=varargin{4}; % 'grey_white' or 'grey' or 'grey_nan'
        
        subjs=length(sn);
        numRuns=numel(run);
        restCond=[29;32]; % SC1: 29th, SC2: 32nd - REST
        
        for s=1:subjs,
            
            VresMS=spm_vol(fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)},'wdResMS.nii'));
            ResMS= spm_read_vols(VresMS);
            ResMS(ResMS==0)=NaN;
            
            % Load over all grey matter mask
            if strcmp(data,'grey'),
                V=spm_vol(fullfile(studyDir{1},'suit','anatomicals',subj_name{sn(s)},'wdc1anatomical.nii')); % call from sc1
            else
                V=spm_vol(fullfile(studyDir{1},'suit','anatomicals','cerebellarGreySUIT.nii')); % call from sc1
            end
            
            X=spm_read_vols(V);
            % Check if V.mat is the the same as wdResMS!!!
            grey_threshold = 0.1; % grey matter threshold
            indx=find(X>grey_threshold);
            [i,j,k]= ind2sub(size(X),indx');
            
            encodeSubjDir = fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),subj_name{sn(s)}); dircheck(encodeSubjDir);
            glmSubjDir =[studyDir{study} sprintf('/GLM_firstlevel_%d/',glm) subj_name{sn(s)}];
            Y=load(fullfile(glmSubjDir,'SPM_info.mat'));
            
            % reorganise
            r=ones(numRuns,1)*restCond(study); % rest is 29th (SC1) and 32nd (SC2) task cond (not zero)
            T=struct('run',[1:numRuns]','sess',kron([1:2]',ones(numRuns/2,1)),'cond',r,'SN',ones(numRuns,1)*sn(s));
            T.TN=repmat({'rest'},numRuns,1);
            Y=addstruct(Y,T);
            Y=rmfield(Y,'task');
            
            switch data,
                case 'grey'
                    % univariately pre-whiten cerebellar voxels
                    nam={};
                    for b=1:length(Y.SN),
                        nam{1}=fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('wdbeta_%2.4d.nii',b));
                        V=spm_vol(nam{1});
                        B1(b,:)=spm_sample_vol(V,i,j,k,0);
                        B1(b,:)=bsxfun(@rdivide,B1(b,:),sqrt(ResMS(indx)')); % univariate noise normalisation
                    end
                    for b=1:length(Y.SN),
                        Yy=zeros(1,V.dim(1)*V.dim(2)*V.dim(3));
                        Yy(1,indx)=B1(b,:);
                        Yy=reshape(Yy,[V.dim(1),V.dim(2),V.dim(3)]);
                        Yy(Yy==0)=NaN;
                        idx=find(~isnan(Yy));
                        Yy=Yy(:);
                        Bb(b,:)=Yy(idx,:);
                    end
                    clear B1 indx
                    B1=Bb;
                    indx=idx;
                case 'grey_nan'
                    load(fullfile(studyDir{study},suitDir,sprintf('glm%d',glm),subj_name{sn(s)},'wdBetas_UW.mat'));
                    for b=1:size(D,1),
                        dat=D(b,:,:,:);
                        Vi.dat=reshape(dat,[V.dim(1),V.dim(2),V.dim(3)]);
                        B1(b,:)=spm_sample_vol(Vi,i,j,k,0);
                    end
            end
            
            % write out new structure ('Y_info')
            Y.data=B1;
            Y.data(end-numRuns+1:end,:)=0; % make intercept (last 16 regressors) equal zero
            Y.identity=indicatorMatrix('identity_p',Y.cond);
            Y.nonZeroInd=repmat(indx',size(B1,1),1);
            
            outName= fullfile(encodeSubjDir,sprintf('Y_info_glm%d_%s.mat',glm,data));
            save(outName,'Y','-v7.3');
            fprintf('cerebellar voxels (%s) computed for %s \n',data,subj_name{sn(s)});
            clear B1 idx Bb indx
        end
    case 'PREP:cortex:indiv_vertices'        % STEP 11.7: Get UW cortical data( vertices)
        sn=varargin{1};
        study=varargin{2};
        glm=varargin{3};
        
        subjs=length(sn);
        numRuns=numel(run);
        
        restCond=[29;32]; % SC1: 29th, SC2: 32nd - REST
        
        for s=1:subjs,
            
            for h=1:2,
                
                encodeSubjDir = fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),subj_name{sn(s)}); dircheck(encodeSubjDir);
                glmSubjDir =[studyDir{study} sprintf('/GLM_firstlevel_%d/',glm) subj_name{sn(s)}];
                Y=load(fullfile(glmSubjDir,'SPM_info.mat'));
                
                % modify existing 'spm_info' structure (include zero intercept)
                r=ones(numRuns,1)*restCond(study); % rest is 29th (SC1) and 32nd (SC2) task cond (not zero)
                T=struct('run',[1:numRuns]','sess',kron([1:2]',ones(numRuns/2,1)),'cond',r,'SN',ones(numRuns,1)*sn(s));
                T.TN=repmat({'rest'},numRuns,1);
                Y=addstruct(Y,T);
                Y=rmfield(Y,'task');
                
                B=caret_load(fullfile(studyDir{study},caretDir,sprintf('x%s',subj_name{sn(s)}),hemName{h},sprintf('%s_glm%d_beta_cortex_%s.metric',subj_name{sn(s)},glm,hem{h})));
                R=caret_load(fullfile(studyDir{study},caretDir,sprintf('x%s',subj_name{sn(s)}),hemName{h},sprintf('%s_glm%d_ResMS_cortex_%s.metric',subj_name{sn(s)},glm,hem{h})));
                
                % write out new structure ('Y_info')
                Y.data=bsxfun(@rdivide,B.data',sqrt(R.data')); % UW betas
                Y.data(end-numRuns+1:end,:)=0; % make intercept (last 16 regressors) equal zero
                Y.nonZeroInd=B.index';
                
                outName=fullfile(encodeSubjDir,sprintf('Y_info_glm%d_cortex_%s.mat',glm,hem{h}));
                save(outName,'Y','-v7.3');
                fprintf('cortical vertices: (Y data) computed for %s \n',subj_name{sn(s)});
                clear B R Y
            end
        end
    case 'PREP:betas_all'                    % STEP 11.8
        sn=varargin{1};
        glm=varargin{2};
        type=varargin{3}; % 'cortex' or 'cerebellum'
        hemN=varargin{4}; % [2] for cortex and [1] for cerebellum
        
        subjs=length(sn);
        
        C=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        for s=1:subjs,
            for study=1:2,
                for h=hemN,
                    numTasks=length(C.condNames(C.StudyNum==study));
                    condNum{study} = [kron(ones(numel(run),1),[1:numTasks]');ones(numel(run),1)*numTasks+1];
                    partNum{study} = [kron([1:numel(run)]',ones(numTasks,1));[1:numel(run)]'];
                    encodeSubjDir=fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),subj_name{sn(s)});
                    switch type,
                        case 'cortex'
                            inFile=sprintf('Y_info_glm%d_%s_%s.mat',glm,type,hem{h});
                            outFile=sprintf('allVert_sc1_sc2_sess_%s.mat',type,hem{h});
                        case 'cerebellum'
                            inFile=sprintf('Y_info_glm%d_grey_nan.mat',glm);
                            outFile=sprintf('allVox_sc1_sc2_sess_%s.mat',type);
                    end
                    load(fullfile(encodeSubjDir,inFile))
                    Yy{s}{study}=Y.data;
                    fprintf('subj%d study%d done \n',sn(s),study)
                end
            end
        end
        save(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),outFile),'Yy','condNum','partNum','-v7.3');
    case 'PREP:IPM_all'                      % FIX THIS TO READ IN 'allVox_sc1_sc2_sess_cerebellum'
        sn=varargin{1};
        glm=varargin{2};
        type=varargin{3}; % 'cortex' or 'cerebellum'
        step=varargin{4}; % 1- 1 IPM (sc1+sc2); 2 - 4 sIPM's (sc1+sc2)
        
        subjs=length(sn);
        
        % prep inputs for PCM modelling functions
        for s=1:subjs,
            V=[];
            load(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),sprintf('allVox_sc1_sc2_sess_%s.mat',type))); % region stats (T)
            for study=1:2,
                D = load(fullfile(studyDir{study},sprintf('GLM_firstlevel_%d',glm), subj_name{sn(s)}, 'SPM_info.mat'));   % load subject's trial structure
                betaW              = Yy{s}{study};
                % get subject's partitions and second moment matrix
                N                  = length(D.run);
                numConds(study)    = length(D.cond(D.cond~=1))/numel(run); % remove instruct
                betaW              = betaW(1:N,:);
                betaW              = [betaW(D.cond~=1,:);zeros(numel(run),size(betaW,2))]; % remove instruct & add rest
                B.partVec          = [D.run(D.cond~=1);[1:numel(run)]'];   % remove instruct & add rest
                B.sessNum          = [D.sess(D.cond~=1);kron([1:2],ones(1,numel(run)/2))'];   % remove instruct & add rest
                % subtract block mean from run
                I = indicatorMatrix('identity',B.partVec);
                R  = eye(N)-I*pinv(I);
                B.Y = R*betaW;            % Subtract block mean
                if study==1,
                    B.condVec      = [repmat([1:numConds(study)],1,numel(run))';repmat(numConds(study)+1,numel(run),1)];  % remove instruct % add rest
                else
                    B.condVec      = [repmat([numConds(2)-1:numConds(1)+numConds(2)+1],1,numel(run))';repmat(numConds(1)+numConds(2)+2,numel(run),1)];
                end
                B.studyNum         = repmat(study,size(betaW,1),1);
                V=addstruct(V,B);
                clear B
            end
            
            % Overall IPM or 4-session IPM's
            if step==1,
                G_hat(:,:,s)    = pcm_estGCrossval(V.Y,V.partVec,V.condVec);  % get IPM
            else
                for study=1:2,
                    for sess=1:2,
                        sessStudyIdx=V.studyNum==study & V.sessNum==sess;
                        if study==1,
                            G_hat_sc1(:,:,sess,s)=pcm_estGCrossval(V.Y(sessStudyIdx,:),V.partVec(sessStudyIdx),V.condVec(sessStudyIdx));
                        else
                            G_hat_sc2(:,:,sess,s)=pcm_estGCrossval(V.Y(sessStudyIdx,:),V.partVec(sessStudyIdx),V.condVec(sessStudyIdx));
                        end
                    end
                end
            end
            
            fprintf('IPM calculated for subj%d \n',sn(s));
            clear Y partVec condVec
        end
        % Save overall IPM or 4-session IPM's
        if step==1,
            save(fullfile(studyDir{2},regDir,sprintf('glm%d',glm),sprintf('G_hat_sc1_sc2_%s.mat',type)),'G_hat');
        else
            save(fullfile(studyDir{2},regDir,sprintf('glm%d',glm),sprintf('G_hat_sc1_sc2_sess_%s.mat',type)),'G_hat_sc1','G_hat_sc2');
        end
    case 'PREP:avrg_betas'                   % STEP 11.9
        sn=varargin{1};
        study=varargin{2}; % studyNum
        glm=varargin{3}; % glmNum
        type=varargin{4}; % 'cortex' or 'cereb'
        hemN=varargin{5}; % [1:2] and [1]
        
        % instructions are removed in this case!! (but rest remains)
        restCond=[29,32];
        
        subjs=length(sn);
        for h=hemN,
            T=[];
            for s=1:subjs,
                switch type,
                    case 'cortex'
                        inFile=sprintf('Y_info_glm%d_cortex_%s.mat',glm,hemN{h});
                        outName=sprintf('%s_%s_avrgDataStruct_vert.mat',type,hem{h});
                    case 'cereb'
                        inFile=sprintf('Y_info_glm%d_grey_nan.mat',glm);
                        outName=sprintf('%s_avrgDataStruct.mat',type);
                        V=spm_vol(fullfile(studyDir{1},suitDir,'anatomicals','cerebellarGreySUIT.nii'));
                end
                load(fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),subj_name{sn(s)},inFile));
                numTasks=size(Y.identity,2);
                %                 Y.cond = Y.cond-1;
                %                 Y.cond(Y.cond==-1)=numTasks;
                numTasks=numTasks-1;
                for sess=1:2,
                    if sess==1,r=1;else r=9;end;
                    v=ones(numTasks,1);
                    S.SN = v*sn(s);
                    S.sess = v*sess;
                    S.study= v*study;
                    S.cond = [1:numTasks]';
                    indx = (Y.sess==sess & Y.cond~=0 & Y.cond~=restCond(study));
                    S.TN   = Y.TN(indx & Y.run==r); % get condNames
                    X=indicatorMatrix('identity',Y.cond(indx,:));
                    S.data = pinv(X)*Y.data(indx,:);
                    S.data = bsxfun(@minus,S.data,nanmean(S.data)); % remove mean within each session
                    R=struct('sess',sess,'cond',restCond(study),'SN',sn(s),'study',study,'TN','rest','data',zeros(1,size(S.data,2))); % add rest
                    S=addstruct(S,R);
                    T=addstruct(T,S);
                end;
                if strcmp(type,'cereb'),
                    volIndx=Y.nonZeroInd(1,:);
                end
                clear Y;
                fprintf('subj%d done...\n',sn(s))
            end;
            % save out struct
            if strcmp(type,'cereb'),
                save(fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),outName),'T','volIndx','V');
            else
                save(fullfile(studyDir{study},encodeDir,sprintf('glm%d',glm),outName),'T');
            end
        end
        
    case 'CHECK:run_times'                   % CHECK: Ensure that start-times (real and predicted) match
        % Check that the startimes match TR times
        % example: sc1_sc2_imana('CHECK:run_times',2,1)
        sn=varargin{1}; % subjNum
        study=varargin{2}; % studyNum
        
        cd(fullfile(studyDir{study},behavDir,subj_name{sn}));
        
        D = dload(sprintf('sc%d_%s.dat',study,subj_name{sn}));
        D = getrow(D,D.runNum>=funcRunNum(1)); % gets behavioural data for functional runs
        %         lineplot(D.realStartTime(D.runNum>=funcRunNum(1)),D.startTime(D.runNum>=funcRunNum(1)));
        pivottable(D.runNum,D.startTime,D.realStartTime,'(x)'); % visually inspect that the start-times match
        taskNames=unique(D.taskName);
        for tt=1:length(taskNames),
            table(D.taskFile(strcmpi(D.taskName,taskNames(tt))),D.runNum(strcmpi(D.taskName,taskNames(tt))));
        end
        fprintf('Visually inspect that the startimes match across all runs \n')
    case 'CHECK:behavioural'                 % CHECK: Check behavioural performance
        % Get behavioural results
        % example: sc1_sc2_imana('CHECK:behavioural',2,1,'scanning','stroop')
        sn=varargin{1};  % subjNum
        study=varargin{2}; % studyNum
        sess=varargin{3}; % 'behavioural' or 'scanning'
        type=varargin{4}; % taskType
        
        switch type
            case 'stroop'
                A=[];
                for s=sn,
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_stroop.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot([A.runNum], A.rt, 'split',A.trialType,'leg',{'incongruent','congruent'},'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot([A.runNum], A.numCorr, 'split',A.trialType,'leg',{'incongruent','congruent'},'subset', A.respMade==1);
                xlabel('Run')
                ylabel('Percent correct')
            case 'nBack'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_nBack.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.possCorr,'leg', {'No Match','Match'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.possCorr,'leg', {'No Match','Match'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'visualSearch'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_visualSearch.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(4,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'absent','present'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(4,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'absent','present'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                
                subplot(4,1,3)
                lineplot(A.runNum, A.rt, 'split', A.setSize,'leg', {'4','8','12'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(4,1,4)
                lineplot(A.runNum, A.numCorr, 'split', A.setSize,'leg', {'4','8','12'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'GoNoGo'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_GoNoGo.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'negative','positive'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'negative','positive'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'nBackPic'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_nBackPic.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.possCorr,'leg', {'No Match','Match'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.possCorr,'leg', {'No Match','Match'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'affective'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_affective.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'unpleasant','pleasant'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'unpleasant','pleasant'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'emotional'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_emotional.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'sad','happy'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'sad','happy'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'ToM'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_ToM.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt,'split',A.condition,'leg',{'false','true'},'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr,'split',A.condition,'leg',{'false','true'},'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'arithmetic'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_arithmetic.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'equations','control'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'equations','control'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'intervalTiming'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_intervalTiming.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                figure(s)
                subplot(2,1,1);
                lineplot(A.runNum,A.rt,'split',A.trialType,'leg', {'long', 'short'},'subset',A.respMade==1);
                xlabel('Run');
                ylabel('Reaction Time');
                
                subplot(2,1,2);
                lineplot(A.runNum,A.numCorr,'split',A.trialType,'leg',{'long','short'},'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
            case 'motorSequence'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc1_%s_motorSequence.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                figure(s)
                subplot(2,1,1);
                lineplot(A.runNum,A.rt,'split',A.trialType,'leg',{'control','sequence'});
                xlabel('Run');
                ylabel('Reaction Time');
                
                subplot(2,1,2);
                lineplot(A.runNum,A.numCorr,'split', A.trialType,'leg', {'control','sequence'});
                xlabel('Run')
                ylabel('Percent Correct')
            case 'CPRO'
                A=[];
                for s=sn,
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_CPRO.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot([A.runNum], A.rt,'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot([A.runNum], A.numCorr,'subset', A.respMade==1);
                xlabel('Run')
                ylabel('Percent correct')
                title(main)
            case 'prediction'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_prediction.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'meaningful','not meaningful'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'meaningful','not meaningful'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'spatialMap'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_spatialMap.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'mentalRotation'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_mentalRotation.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'emotionProcess'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_emotionProcess.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'intact','scrambled'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'intact','scrambled'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'respAlt'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_respAlt.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'easy','med','diff'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'visualSearch2'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_visualSearch2.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(4,1,1)
                lineplot(A.runNum, A.rt, 'split', A.trialType,'leg', {'absent','present'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(4,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.trialType,'leg', {'absent','present'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
                
                subplot(4,1,3)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'4','8','12'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(4,1,4)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'4','8','12'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'nBackPic2'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_nBackPic2.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt, 'split', A.condition,'leg', {'response','no response'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr, 'split', A.condition,'leg', {'response','no response'}, 'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'ToM2'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_ToM2.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1)
                lineplot(A.runNum, A.rt,'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Reaction Time')
                title(main)
                
                subplot(2,1,2)
                lineplot(A.runNum, A.numCorr,'subset',A.respMade==1);
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
            case 'motorSequence2'
                A=[];
                for s=sn,
                    
                    D = dload(fullfile(studyDir{study},behavDir,subj_name{s},sprintf('sc2_%s_motorSequence2.dat',subj_name{s})));
                    A = addstruct(A,D);
                    switch sess,
                        case 'behavioural'
                            A = getrow(A,A.runNum>=1 & A.runNum<=51);
                        case 'scanning'
                            A = getrow(A,A.runNum>=funcRunNum(1) & A.runNum<=funcRunNum(2));
                    end
                end
                
                if numel(sn)>1,
                    main=sprintf('all subjects-%s',type);
                else
                    main=sprintf('%s-%s',subj_name{s},type);
                end
                
                figure(s)
                subplot(2,1,1);
                lineplot(A.runNum,A.rt,'split',A.condition,'leg',{'control','sequence'});
                xlabel('Run');
                ylabel('Reaction Time');
                title(main)
                
                subplot(2,1,2);
                lineplot(A.runNum,A.numCorr,'split', A.condition,'leg', {'control','sequence'});
                xlabel('Run')
                ylabel('Percent Correct')
                title(main)
        end
    case 'CHECK:designMatrix'                % CHECK: Visually inspect design matrix
        % example: sc1_sc2_imana('CHECK:designMatrix',2,1,4)
        % run 'make_glm#','estimate_glm' and 'contrast' before this step
        sn=varargin{1};
        study=varargin{2};
        glm=varargin{3};
        
        glmDir =[studyDir{study} sprintf('/GLM_firstlevel_%d/',glm) subj_name{sn}];
        
        A=load(fullfile(glmDir,'SPM.mat'));
        
        X=A.SPM.xX.X(A.SPM.Sess(1).row,A.SPM.Sess(1).col);
        
        figure(sn);
        imagesc(X);
        title('Design Matrix for run 1')
        
        figure(sn+1);
        indx=A.SPM.Sess(1).col;
        imagesc(A.SPM.xX.Bcov(indx,indx))
        title('variance of predicted betas')
        
        if length(varargin)==3,
            reg=varargin{3}; % specify [10,11] for example
            figure(3);
            plot([1:size(X,1)],X(:,reg),'b'); % which regressors are we plotting?
        end
    case 'CHECK:movement'                    % CHECK: Check residual SD
        % example: sc1_sc2_imana('CHECK:movement',2,1,4)
        s=varargin{1}; % subjNum
        study=varargin{2}; % studyNum
        glm=varargin{3}; % glmNum
        
        glmSubjDir =[studyDir{study} sprintf('/GLM_firstlevel_%d/',glm)  subj_name{s}];dircheck(glmSubjDir);
        load(fullfile(glmSubjDir,'SPM.mat'));
        spm_rwls_resstats(SPM);
        
        % 'spm_rwls_resstats' was edited to include check for FAST GLM
        % ('spm_light') - no 'xVi' field
    case 'CHECK:SNR'                         % CHECK: Test pattern consistencies (across subs,across sessions,across regions)
        % example: sc1_sc2_imana('CHECK:SNR',[1,2,3],1,4,{'cerebellum_grey'},[1],'subjs')
        % default noise normalisation is univariate. For NW or MW add 1 or
        % 3 respectively.
        % For example: sc1_sc2_imana('CHECK:SNR',2,1,4,{'cerebellum_grey'},{[1]},'subjs')
        sn=varargin{1};
        study=varargin{2}; % studyNum
        glm=varargin{3}; % glmNum
        types=varargin{4}; % {'cortex_grey','cerebellum_grey'}
        reg=varargin{5}; % which region of the cortex (lh or rh) or cortical_lobes (1,2,3,4) or cerebellar lobules (1-10)?
        step=varargin{6};
        
        subjs=length(sn);
        
        A=[];
        idx=1;
        for s=1:subjs,
            idx2=1;
            for i=1:numel(types),
                for r=1:numel(reg{i}),
                    % load statistics for all subjects and all GLMs
                    load(fullfile(studyDir{study},regDir,sprintf('glm%d',glm),subj_name{sn(s)},sprintf('Toverall_%s.mat',types{i})));
                    A=getrow(To,To.method_num==2 & To.region==reg{i}(r));
                    sess1(idx)=A.R2_Sess1;
                    sess2(idx)=A.R2_Sess2;
                    overall(idx)=A.R2_overall;
                    S(idx)=sn(s);
                    R(idx)=idx2;
                    regName{idx}=sprintf('%d.%s',idx,char(A.regName));
                    idx=idx+1;
                    idx2=idx2+1;
                end
            end
        end
        
        switch step,
            case 'subjs'
                figure()
                barplot(S',overall','split',R','leg',regName')
                ylabel('R2')
            case 'regions'
                figure()
                barplot(R',overall','split',R','leg',regName')
                ylabel('R2')
            case 'sess1'
                figure()
                barplot(R',sess1','split',R','leg',regName')
                ylabel('R2')
            case 'sess2'
                figure()
                barplot(R',sess2','split',R','leg',regName')
                ylabel('R2')
        end
    case 'CHECK:DIST'                        % CHECK: Within and across study reliability of distances
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or 'cortex'
        % example 'sc1_sc2_imana('CHECK:DIST',4,'cerebellum')
        
        load(fullfile(regDir,sprintf('glm%d',glm),sprintf('G_hat_sc1_sc2_sess_%s.mat',type)));
        D=dload('sc1_sc2_taskConds.txt');
        D1=getrow(D,D.StudyNum==1);
        D2=getrow(D,D.StudyNum==2);
        
        % Look at the shared conditions only
        i1 = find(D1.overlap==1);
        i2 = find(D2.overlap==1);
        [~,b] = sort(D2.condNumUni(i2));     % Bring the indices for sc2 into the right order.
        i2=i2(b);
        numCond = length(i1);
        numSubj = size(G_hat_sc1,4);
        numSess = 2;
        
        C=indicatorMatrix('allpairs',[1:numCond]);
        for i=1:numSubj
            for j=1:numSess
                dist(:,j  ,i)  = ssqrt(diag(C*G_hat_sc1(i1,i1,j,i)*C'));
                dist(:,j+2,i)  = ssqrt(diag(C*G_hat_sc2(i2,i2,j,i)*C'));
            end;
            CORR(:,:,i)    = corr(dist(:,:,i));
            T.SN(i,1)      = i;
            T.within1(i,1) = CORR(1,2,i);
            T.within2(i,1) = CORR(3,4,i);
            T.across(i,1)  = mean(mean(CORR(1:2,3:4,i)));
        end;
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_tukey');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
        % group reliability
        X=nanmean(dist,3);
        groupCorr=corr(X);
        fprintf('group reliability for study1 is %2.2f \n',groupCorr(1,2));
        fprintf('group reliability for study2 is %2.2f \n',groupCorr(3,4));
        fprintf('group reliability for shared tasks is %2.2f \n',mean(mean(groupCorr(1:2,3:4))));
        varargout={T};
    case 'CHECK:PATTERN'                     % CHECK: Within and across study reliability of patterns
        glm=varargin{1};
        type=varargin{2}; % 'cerebellum' or 'cortex'
        
        D=dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        
        load(fullfile(studyDir{2},encodeDir,sprintf('glm%d',glm),sprintf('allVox_sc1_sc2_sess_%s.mat',type)));
        Y=Yy;clear Yy;
        
        numSubj=length(Y);
        
        for subj=1:numSubj,
            idx=1;
            for study=1:2,
                
                D1=getrow(D,D.StudyNum==study);
                
                cN=condNum{study}-1;  % Important: In the allVox file, instruction is still included!
                pN=partNum{study};    % Partition Numner
                sN=(pN>8)+1;          % Sessions Number
                for se=1:2,
                    X1=indicatorMatrix('identity_p',cN.*(sN==se));  % This one is the matrix that related trials-> condition numbers
                    X2=indicatorMatrix('identity_p',D1.condNumUni.*D1.overlap); % THis goes from condNum to shared condNumUni
                    Yf(:,:,idx,subj)=pinv(X1*X2)*Y{subj}{study};
                    Yf(:,:,idx,subj)=bsxfun(@minus,Yf(:,:,idx,subj),mean(Yf(:,:,idx,subj)));
                    idx=idx+1;
                end;
            end;
            CORR(:,:,subj)=interSess_corr(Yf(:,:,:,subj));
            T.SN(subj,1)      = subj;
            T.within1(subj,1) = CORR(1,2,subj);
            T.within2(subj,1) = CORR(3,4,subj);
            T.across(subj,1)  = mean(mean(CORR(1:2,3:4,subj)));
        end
        
        % within & between-subj reliability
        myboxplot([],[T.within1 T.within2 T.across],'style_tukey');
        drawline(0,'dir','horz');
        set(gca,'XTickLabel',{'SC1','SC2','across'});
        set(gcf,'PaperPosition',[2 2 4 3]);
        wysiwyg;
        ttest(sqrt(T.within1.*T.within2),T.across,2,'paired');
        
        % group reliability
        X=nanmean(CORR,3);
        fprintf('group reliability for study1 is %2.2f \n',X(1,2));
        fprintf('group reliability for study2 is %2.2f \n',X(3,4));
        fprintf('group reliability for shared tasks is %2.2f \n',mean(mean(X(1:2,3:4))));
        varargout={T};
    case 'CHECK:saccades'                    % CHECK: Pivottable (average saccades across tasks)
        sn=varargin{1};
        glm=varargin{2};
        Rr=[];
        
        switch glm,
            case 13
                taskNames = taskNames;
            case 4
                taskNames=condNames;
        end
        
        for s=sn,
            eyeSubjDir=fullfile(eyeDir,subj_name{s});
            load(fullfile(eyeSubjDir,sprintf('glm%d_saccades.mat',glm)));  % load saccade info
            Rr=addstruct(Rr,R);
        end
        
        lineplot(Rr.condNum,Rr.avgSacc,'split',Rr.taskName,'leg',taskNames,'markersize',10);
        xlabel('Tasks')
        ylabel('Number of Saccades')
        varargout={Rr};
end

% Local functions
function dircheck(dir)
if ~exist(dir,'dir');
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);
end

% InterSubj Corr
function C=interSess_corr(Y)
numSess=size(Y,3);
for i=1:numSess
    for j=1:numSess
        C(i,j)=nansum(nansum(Y(:,:,i).*Y(:,:,j)))/...
            sqrt(nansum(nansum(Y(:,:,i).*Y(:,:,i)))*...
            nansum(nansum(Y(:,:,j).*Y(:,:,j))));
    end;
end



