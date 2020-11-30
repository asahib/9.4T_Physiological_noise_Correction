datadir=('/p/ashish/MultiBandEPI/EPILEPSY_STUDY');
subj={'PAT_01','PAT_02','PAT_03','PAT_04','PAT_05','PAT_06','PAT_07',};
restdir_sw=('/p/ashish/MultiBandEPI/EPILEPSY_STUDY/PAT_01/DY_BC_3d');
restdir_ROI=('/p/ashish/MultiBandEPI/EPILEPSY_STUDY/PAT_01/op_tc');

MaskFilename='/p/ashish/MultiBandEPI/EPILEPSY_STUDY/PAT_01/T1/wc1s04TT890518MRMB-0002-00001-000176-01.nii';

ROIDef_wm={'/p/ashish/MultiBandEPI/EPILEPSY_STUDY/PAT_01/T1/wc2s04TT890518MRMB-0002-00001-000176-01.nii'};

ROIDef_csf={'/p/ashish/MultiBandEPI/EPILEPSY_STUDY/PAT_01/T1/wc3s04TT890518MRMB-0002-00001-000176-01.nii'};

rp=load(spm_select('FPList',fullfile(datadir,subj{1},'topup'),'^rp_cT.*txt$'));
rp_dt=detrend(rp,'constant');


[theROITimeCourses_wm] = rest_ExtractROITC(restdir_sw, ROIDef_wm,restdir_ROI);

[theROITimeCourses_csf] = rest_ExtractROITC(restdir_sw, ROIDef_csf,restdir_ROI);


SCovar=[theROITimeCourses_wm-mean(theROITimeCourses_wm) theROITimeCourses_csf-mean(theROITimeCourses_csf) rp(1:end,:)];
% SCovar=[theROITimeCourses_wm-mean(theROITimeCourses_wm) theROITimeCourses_csf-mean(theROITimeCourses_csf) rp(1:77,:)];

Subject_Covariables.ort_file=fullfile(restdir_ROI, 'wcT_WM_rp_csfSignals.txt');

DirPostfix=['rp_wm_csf' '_removed'];
Subject_Covariables.polort=-1;
save(Subject_Covariables.ort_file,'SCovar', '-ASCII', '-DOUBLE','-TABS');
rest_RegressOutCovariates(restdir_sw,Subject_Covariables,DirPostfix,MaskFilename);


% disp(datestr(datenum(0,0,0,0,0,f),'HH: MM:SS'))