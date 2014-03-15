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

GaussSig_cat0 = Gaussian(mgg, mgg_sig_m0_cat0, mgg_sig_gsigma_cat0);
CBSig_cat0 = CBShape(mgg, mgg_sig_m0_cat0, mgg_sig_sigma_cat0, mgg_sig_alpha_cat0, mgg_sig_n_cat0);
mggSig_cat0 = AddPdf(GaussSig_cat0, CBSig_cat0, mgg_sig_frac_cat0);

GaussSig_cat1 = Gaussian(mgg, mgg_sig_m0_cat1, mgg_sig_gsigma_cat1);
CBSig_cat1 = CBShape(mgg, mgg_sig_m0_cat1, mgg_sig_sigma_cat1, mgg_sig_alpha_cat1, mgg_sig_n_cat1);
mggSig_cat1 = AddPdf(GaussSig_cat1, CBSig_cat1, mgg_sig_frac_cat1);

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

GaussHig_0_cat0 = Gaussian(mgg, mgg_hig_m0_0_cat0, mgg_hig_gsigma_0_cat0);
CBHig_0_cat0 = CBShape(mgg, mgg_hig_m0_0_cat0, mgg_hig_sigma_0_cat0, mgg_hig_alpha_0_cat0, mgg_hig_n_0_cat0);
mggHig_0_cat0 = AddPdf(GaussHig_0_cat0, CBHig_0_cat0, mgg_hig_frac_0_cat0);

GaussHig_0_cat1 = Gaussian(mgg, mgg_hig_m0_0_cat1, mgg_hig_gsigma_0_cat1);
CBHig_0_cat1 = CBShape(mgg, mgg_hig_m0_0_cat1, mgg_hig_sigma_0_cat1, mgg_hig_alpha_0_cat1, mgg_hig_n_0_cat1);
mggHig_0_cat1 = AddPdf(GaussHig_0_cat1, CBHig_0_cat1, mgg_hig_frac_0_cat1);

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

GaussHig_1_cat0 = Gaussian(mgg, mgg_hig_m0_1_cat0, mgg_hig_gsigma_1_cat0);
CBHig_1_cat0 = CBShape(mgg, mgg_hig_m0_1_cat0, mgg_hig_sigma_1_cat0, mgg_hig_alpha_1_cat0, mgg_hig_n_1_cat0);
mggHig_1_cat0 = AddPdf(GaussHig_1_cat0, CBHig_1_cat0, mgg_hig_frac_1_cat0);

GaussHig_1_cat1 = Gaussian(mgg, mgg_hig_m0_1_cat1, mgg_hig_gsigma_1_cat1);
CBHig_1_cat1 = CBShape(mgg, mgg_hig_m0_1_cat1, mgg_hig_sigma_1_cat1, mgg_hig_alpha_1_cat1, mgg_hig_n_1_cat1);
mggHig_1_cat1 = AddPdf(GaussHig_1_cat1, CBHig_1_cat1, mgg_hig_frac_1_cat1);

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

GaussHig_2_cat0 = Gaussian(mgg, mgg_hig_m0_2_cat0, mgg_hig_gsigma_2_cat0);
CBHig_2_cat0 = CBShape(mgg, mgg_hig_m0_2_cat0, mgg_hig_sigma_2_cat0, mgg_hig_alpha_2_cat0, mgg_hig_n_2_cat0);
mggHig_2_cat0 = AddPdf(GaussHig_2_cat0, CBHig_2_cat0, mgg_hig_frac_2_cat0);

GaussHig_2_cat1 = Gaussian(mgg, mgg_hig_m0_2_cat1, mgg_hig_gsigma_2_cat1);
CBHig_2_cat1 = CBShape(mgg, mgg_hig_m0_2_cat1, mgg_hig_sigma_2_cat1, mgg_hig_alpha_2_cat1, mgg_hig_n_2_cat1);
mggHig_2_cat1 = AddPdf(GaussHig_2_cat1, CBHig_2_cat1, mgg_hig_frac_2_cat1);

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

GaussHig_3_cat0 = Gaussian(mgg, mgg_hig_m0_3_cat0, mgg_hig_gsigma_3_cat0);
CBHig_3_cat0 = CBShape(mgg, mgg_hig_m0_3_cat0, mgg_hig_sigma_3_cat0, mgg_hig_alpha_3_cat0, mgg_hig_n_3_cat0);
mggHig_3_cat0 = AddPdf(GaussHig_3_cat0, CBHig_3_cat0, mgg_hig_frac_3_cat0);

GaussHig_3_cat1 = Gaussian(mgg, mgg_hig_m0_3_cat1, mgg_hig_gsigma_3_cat1);
CBHig_3_cat1 = CBShape(mgg, mgg_hig_m0_3_cat1, mgg_hig_sigma_3_cat1, mgg_hig_alpha_3_cat1, mgg_hig_n_3_cat1);
mggHig_3_cat1 = AddPdf(GaussHig_3_cat1, CBHig_3_cat1, mgg_hig_frac_3_cat1);


