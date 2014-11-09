mgg[100,180];


mgg_bkg_8TeV_slope1_cat0[9,-10.0, 20.0];
mgg_bkg_8TeV_slope1_cat1[9,-10.0, 20.0];

mgg_sig_m0_cat0[124.2, 123, 125];
mgg_sig_sigma_cat0[2.0, 1.0, 3.0];
mgg_sig_alpha_cat0[1.0, 0.0, 2.5];
mgg_sig_n_cat0[2.0, 1.0, 5.0];
mgg_sig_gsigma_cat0[1.2, 0.8, 1.6];
mgg_sig_frac_cat0[0.6, 0.4, 1.0];

mgg_sig_m0_cat1[124.2, 123, 125];
mgg_sig_sigma_cat1[2.0, 1.0, 3.0];
mgg_sig_alpha_cat1[1.0, 0.0, 2.5];
mgg_sig_n_cat1[2.0, 1.5, 10];
mgg_sig_gsigma_cat1[1.2, 0.8, 1.6];
mgg_sig_frac_cat1[0.6, 0.4, 1.0];

mggGaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
mggCBSig_cat0 = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
mggSig_cat0 = AddPdf(mggGaussSig_cat0, mggCBSig_cat0, mgg_sig_frac_cat0);

mggGaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
mggCBSig_cat1 = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
mggSig_cat1 = AddPdf(mggGaussSig_cat1, mggCBSig_cat1, mgg_sig_frac_cat1);

mgg_hig_m0_0_cat0[124.2, 123, 125];
mgg_hig_sigma_0_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat0[1.0, 0.0, 2.5];
mgg_hig_n_0_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_0_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_0_cat1[124.2, 123, 125];
mgg_hig_sigma_0_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_0_cat1[1.0, 0.0, 2.5];
mgg_hig_n_0_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_0_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_0_cat1[0.6, 0.4, 1.0];

mggGaussHig_0_cat0 = Gaussian(mgg, mgg_hig_m0_0_cat0, mgg_hig_gsigma_0_cat0);
mggCBHig_0_cat0 = CBShape(mgg, mgg_hig_m0_0_cat0, mgg_hig_sigma_0_cat0, mgg_hig_alpha_0_cat0, mgg_hig_n_0_cat0);
mggHig_0_cat0 = AddPdf(mggGaussHig_0_cat0, mggCBHig_0_cat0, mgg_hig_frac_0_cat0);

mggGaussHig_0_cat1 = Gaussian(mgg, mgg_hig_m0_0_cat1, mgg_hig_gsigma_0_cat1);
mggCBHig_0_cat1 = CBShape(mgg, mgg_hig_m0_0_cat1, mgg_hig_sigma_0_cat1, mgg_hig_alpha_0_cat1, mgg_hig_n_0_cat1);
mggHig_0_cat1 = AddPdf(mggGaussHig_0_cat1, mggCBHig_0_cat1, mgg_hig_frac_0_cat1);

mgg_hig_m0_1_cat0[124.2, 123, 125];
mgg_hig_sigma_1_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat0[1.0, 0.0, 2.5];
mgg_hig_n_1_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_1_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_1_cat1[124.2, 123, 125];
mgg_hig_sigma_1_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_1_cat1[1.0, 0.0, 2.5];
mgg_hig_n_1_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_1_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_1_cat1[0.6, 0.4, 1.0];

mggGaussHig_1_cat0 = Gaussian(mgg, mgg_hig_m0_1_cat0, mgg_hig_gsigma_1_cat0);
mggCBHig_1_cat0 = CBShape(mgg, mgg_hig_m0_1_cat0, mgg_hig_sigma_1_cat0, mgg_hig_alpha_1_cat0, mgg_hig_n_1_cat0);
mggHig_1_cat0 = AddPdf(mggGaussHig_1_cat0, mggCBHig_1_cat0, mgg_hig_frac_1_cat0);

mggGaussHig_1_cat1 = Gaussian(mgg, mgg_hig_m0_1_cat1, mgg_hig_gsigma_1_cat1);
mggCBHig_1_cat1 = CBShape(mgg, mgg_hig_m0_1_cat1, mgg_hig_sigma_1_cat1, mgg_hig_alpha_1_cat1, mgg_hig_n_1_cat1);
mggHig_1_cat1 = AddPdf(mggGaussHig_1_cat1, mggCBHig_1_cat1, mgg_hig_frac_1_cat1);

mgg_hig_m0_2_cat0[124.2, 123, 125];
mgg_hig_sigma_2_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat0[1.0, 0.0, 2.5];
mgg_hig_n_2_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_2_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_2_cat1[124.2, 123, 125];
mgg_hig_sigma_2_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_2_cat1[1.0, 0.0, 2.5];
mgg_hig_n_2_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_2_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_2_cat1[0.6, 0.4, 1.0];

mggGaussHig_2_cat0 = Gaussian(mgg, mgg_hig_m0_2_cat0, mgg_hig_gsigma_2_cat0);
mggCBHig_2_cat0 = CBShape(mgg, mgg_hig_m0_2_cat0, mgg_hig_sigma_2_cat0, mgg_hig_alpha_2_cat0, mgg_hig_n_2_cat0);
mggHig_2_cat0 = AddPdf(mggGaussHig_2_cat0, mggCBHig_2_cat0, mgg_hig_frac_2_cat0);

mggGaussHig_2_cat1 = Gaussian(mgg, mgg_hig_m0_2_cat1, mgg_hig_gsigma_2_cat1);
mggCBHig_2_cat1 = CBShape(mgg, mgg_hig_m0_2_cat1, mgg_hig_sigma_2_cat1, mgg_hig_alpha_2_cat1, mgg_hig_n_2_cat1);
mggHig_2_cat1 = AddPdf(mggGaussHig_2_cat1, mggCBHig_2_cat1, mgg_hig_frac_2_cat1);

mgg_hig_m0_3_cat0[124.2, 123, 125];
mgg_hig_sigma_3_cat0[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat0[1.0, 0.0, 2.5];
mgg_hig_n_3_cat0[2.0, 1.0, 5.0];
mgg_hig_gsigma_3_cat0[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat0[0.6, 0.4, 1.0];

mgg_hig_m0_3_cat1[124.2, 123, 125];
mgg_hig_sigma_3_cat1[2.0, 1.0, 3.0];
mgg_hig_alpha_3_cat1[1.0, 0.0, 2.5];
mgg_hig_n_3_cat1[2.0, 1.5, 10];
mgg_hig_gsigma_3_cat1[1.2, 0.8, 1.6];
mgg_hig_frac_3_cat1[0.6, 0.4, 1.0];

mggGaussHig_3_cat0 = Gaussian(mgg, mgg_hig_m0_3_cat0, mgg_hig_gsigma_3_cat0);
mggCBHig_3_cat0 = CBShape(mgg, mgg_hig_m0_3_cat0, mgg_hig_sigma_3_cat0, mgg_hig_alpha_3_cat0, mgg_hig_n_3_cat0);
mggHig_3_cat0 = AddPdf(mggGaussHig_3_cat0, mggCBHig_3_cat0, mgg_hig_frac_3_cat0);

mggGaussHig_3_cat1 = Gaussian(mgg, mgg_hig_m0_3_cat1, mgg_hig_gsigma_3_cat1);
mggCBHig_3_cat1 = CBShape(mgg, mgg_hig_m0_3_cat1, mgg_hig_sigma_3_cat1, mgg_hig_alpha_3_cat1, mgg_hig_n_3_cat1);
mggHig_3_cat1 = AddPdf(mggGaussHig_3_cat1, mggCBHig_3_cat1, mgg_hig_frac_3_cat1);


mjj[60,180];
mjj_sig_m0[110.0, 70, 160];
mjj_sig_sigma[10.0, 5.0, 20.0];
mjj_sig_alpha[1.0, 0.5, 3]; 
mjj_sig_n[4.0, 0.5, 10]; 
mjj_sig_gsigma[25, 10.0, 50.0];
mjj_sig_frac[0.3, 0, 0.5];

mjjGaussSig = Gaussian(mjj, mjj_sig_m0, mjj_sig_gsigma);
mjjCBSig    = CBShape(mjj, mjj_sig_m0, mjj_sig_sigma, mjj_sig_alpha, mjj_sig_n);
mjjSig      = AddPdf(mjjGaussSig, mjjCBSig, mjj_sig_frac);

mjj_sig_m0_cat0[110.0, 105, 155];
mjj_sig_sigma_cat0[10.0, 5.0, 20.0];
mjj_sig_alpha_cat0[2.0, 1.0, 2.5]; 
mjj_sig_n_cat0[2.0, 1.0, 5.0]; 
mjj_sig_gsigma_cat0[25.0, 10.0, 50.0];
mjj_sig_frac_cat0[0.1, 0, 0.5];


mjjGaussSig_cat0 = Gaussian(mjj, mjj_sig_m0_cat0, mjj_sig_gsigma_cat0);
mjjCBSig_cat0    = CBShape(mjj, mjj_sig_m0_cat0, mjj_sig_sigma_cat0, mjj_sig_alpha_cat0, mjj_sig_n_cat0);
mjjSig_cat0      = AddPdf(mjjGaussSig_cat0, mjjCBSig_cat0, mjj_sig_frac_cat0);

mjj_sig_m0_cat1[110.0, 70, 160];
mjj_sig_sigma_cat1[10.0, 5.0, 20.0];
mjj_sig_alpha_cat1[2.0, 1.2, 5]; 
mjj_sig_n_cat1[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat1[25.0, 10.0, 50.0];
mjj_sig_frac_cat1[0.1, 0, 0.5];

mjjGaussSig_cat1 = Gaussian(mjj, mjj_sig_m0_cat1, mjj_sig_gsigma_cat1);
mjjCBSig_cat1    = CBShape(mjj, mjj_sig_m0_cat1, mjj_sig_sigma_cat1, mjj_sig_alpha_cat1, mjj_sig_n_cat1);
mjjSig_cat1      = AddPdf(mjjGaussSig_cat1, mjjCBSig_cat1, mjj_sig_frac_cat1);

mjj_sig_m0_cat2[110.0, 70, 160];
mjj_sig_sigma_cat2[10.0, 5.0, 20.0];
mjj_sig_alpha_cat2[1.0, 1.0, 3.5]; 
mjj_sig_n_cat2[2.0, 1.5, 10]; 
mjj_sig_gsigma_cat2[25.0, 10.0, 25.0];
mjj_sig_frac_cat2[0.1, 0, 0.3];

mjjGaussSig_cat2 = Gaussian(mjj, mjj_sig_m0_cat2, mjj_sig_gsigma_cat2);
mjjCBSig_cat2    = CBShape(mjj, mjj_sig_m0_cat2, mjj_sig_sigma_cat2, mjj_sig_alpha_cat2, mjj_sig_n_cat2);
mjjSig_cat2      = AddPdf(mjjGaussSig_cat2, mjjCBSig_cat2, mjj_sig_frac_cat2);

mjj_bkg_8TeV_slope[1.0,0, 1];
mjj_bkg_8TeV_slope1[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope2[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope3[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope4[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope5[0.1,-10.0, 10.0];

mjj_bkg_8TeV_slope_cat0[1.0,0, 1];
mjj_bkg_8TeV_slope1_cat0[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope2_cat0[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope3_cat0[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope4_cat0[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope5_cat0[0.1,-10.0, 10.0];

mjj_bkg_8TeV_slope_cat1[1.0,0, 1];
mjj_bkg_8TeV_slope1_cat1[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope2_cat1[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope3_cat1[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope4_cat1[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope5_cat1[0.1,-10.0, 10.0];

mjj_bkg_8TeV_slope_cat2[1.0,0, 1];
mjj_bkg_8TeV_slope1_cat2[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope2_cat2[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope3_cat2[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope4_cat2[0.1,-10.0, 10.0];
mjj_bkg_8TeV_slope5_cat2[0.1,-10.0, 10.0];

