import ROOT

file_path_wot = "/home/ZhangHaolin/ACTS/acts/RunSpace/wot_time_stddev_2/trackstates_ckf.root"
file_path_wt = "/home/ZhangHaolin/ACTS/acts/RunSpace/wt_time_stddev_2/trackstates_ckf.root"
file_wot = ROOT.TFile.Open(file_path_wot, "READ")
file_wt = ROOT.TFile.Open(file_path_wt, "READ")

c1 = ROOT.TCanvas("c1", "Canvas for DUT Resolution loc0", 800, 600)

trackstates_wot= file_wot.Get("trackstates")
trackstates_wot.Draw("res_eLOC0_ubs>>hist_DUT_Res_wot(150, -0.2, 0.2)", "layer_id==8")
hist_DUT_Res_wot = ROOT.gDirectory.Get("hist_DUT_Res_wot")

trackstates_wt= file_wt.Get("trackstates")
trackstates_wt.Draw("res_eLOC0_ubs>>hist_DUT_Res_wt(150, -0.2, 0.2)", "layer_id==8", "same")
hist_DUT_Res_wt = ROOT.gDirectory.Get("hist_DUT_Res_wt")

#hist_DUT_Res.SetTitle("DUT Resolution-loc0 with Time")

if hist_DUT_Res_wot:
    gaussFit = ROOT.TF1("gaussFit", "gaus", -0.2, 0.2)
    gaussFit.SetParLimits(0, 2400, 3000)
    #gaussFit.SetParLimits(1, -0.1, 0.1)
    #gaussFit.SetParLimits(2, 0, 0.1)
    hist_DUT_Res_wot.Fit(gaussFit)
c1.Update()
    #c1.SaveAs("DUT-Res-eLOC0-wt.pdf")
    #c1.SaveAs("DUT-Res-eLOC0-wt.png")

if hist_DUT_Res_wt:
    gaussFit = ROOT.TF1("gaussFit", "gaus", -0.2, 0.2)
    gaussFit.SetParLimits(0, 1500, 2000)
    #gaussFit.SetParLimits(1, -0.1, 0.1)
    #gaussFit.SetParLimits(2, 0, 0.1)
    hist_DUT_Res_wt.Fit(gaussFit)
c1.Update()
    #c1.SaveAs("DUT-Res-eLOC0-wt.pdf")
    #c1.SaveAs("DUT-Res-eLOC0-wt.png")

file_wot.Close()
file_wt.Close()
input()
