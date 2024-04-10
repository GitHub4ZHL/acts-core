import ROOT

# Get root files
file_path_wot = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/wot_multiplicity_20/trackstates_ckf.root"
file_path_wt = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/wt_multiplicity_20/trackstates_ckf.root"
file_wot = ROOT.TFile.Open(file_path_wot, "READ")
file_wt = ROOT.TFile.Open(file_path_wt, "READ")

# Make Canvas
c1 = ROOT.TCanvas("c1", "Canvas for DUT Resolution loc1", 800, 600)

# Get histogram
trackstates_wot= file_wot.Get("trackstates")
trackstates_wot.Draw("res_eLOC1_ubs>>hist_DUT_Res_wot(150, -0.2, 0.2)", "layer_id==8 && nMeasurementsExcluded>=4")
hist_DUT_Res_wot = ROOT.gDirectory.Get("hist_DUT_Res_wot")
#hist_DUT_Res_wot.SetStats(False)
hist_DUT_Res_wot.GetXaxis().SetTitle("Resolution")
hist_DUT_Res_wot.GetYaxis().SetTitle("Entries")

trackstates_wt= file_wt.Get("trackstates")
trackstates_wt.Draw("res_eLOC1_ubs>>hist_DUT_Res_wt(150, -0.2, 0.2)", "layer_id==8 && nMeasurementsExcluded>=4")
hist_DUT_Res_wt = ROOT.gDirectory.Get("hist_DUT_Res_wt")

hist_DUT_Res_wot.SetTitle("DUT loc1 resolution w/wo time")

# Gauss fit
'''
gaussFit_wot = ROOT.TF1("gaussFit", "gaus", -0.2, 0.2)
gaussFit_wot.SetParLimits(0, 1600, 2000)
gaussFit_wot.SetLineColor(ROOT.kBlue)
gaussFit_wot.GetXaxis().SetTitle("Resolution")
gaussFit_wot.GetYaxis().SetTitle("Entries")

hist_DUT_Res_wot.Fit(gaussFit_wot)
mean_wot = gaussFit_wot.GetParameter(1)
sigma_wot = gaussFit_wot.GetParameter(2)

gaussFit_wt = ROOT.TF1("gaussFit", "gaus", -0.2, 0.2)
gaussFit_wt.SetParLimits(0, 1050, 1300)
gaussFit_wt.SetLineColor(ROOT.kRed)

hist_DUT_Res_wt.Fit(gaussFit_wt)
mean_wt = gaussFit_wt.GetParameter(1)
sigma_wt = gaussFit_wt.GetParameter(2)
gaussFit_wt.GetXaxis().SetTitle("Resolution")
gaussFit_wt.GetYaxis().SetTitle("Entries")
'''
hist_DUT_Res_wot.Draw()
hist_DUT_Res_wt.Draw("same")

# Legend
'''
legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
legend.SetTextSize(0.025)
legend.SetFillStyle(1) # Transparent

legend.AddEntry(hist_DUT_Res_wot, f"mean:{mean_wot:.0f} stddev: {sigma_wot:.3f} mm", "L")
legend.AddEntry(hist_DUT_Res_wt, f"mean: {mean_wt:.0f} stddev: {sigma_wt:.3f} mm", "L")
legend.Draw()
'''
c1.Update()
#c1.SaveAs("DUT-Res-eLOC1-digi.pdf")
c1.SaveAs("DUT-Res-eLOC1.png")

file_wot.Close()
file_wt.Close()

input()
