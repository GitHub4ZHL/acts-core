import ROOT

file_path = "/home/ZhangHaolin/ACTS/acts/RunSpace/wot_time_stddev_2/trackstates_ckf.root"
file = ROOT.TFile.Open(file_path, "READ")

c1 = ROOT.TCanvas("c1", "Canvas for DUT Resolution loc0", 800, 600)

trackstates= file.Get("trackstates")
trackstates.Draw("res_eLOC0_smt>>hist_DUT_Res(150, -0.5, 0.5)", "layer_id==8")

hist_DUT_Res = ROOT.gDirectory.Get("hist_DUT_Res")
hist_DUT_Res.SetTitle("DUT Resolution-loc0 without Time")

if hist_DUT_Res:

    gaussFit = ROOT.TF1("gaussFit", "gaus", -0.5, 0.5)
    gaussFit.SetParLimits(0, 6400, 8000)
    #gaussFit.SetParLimits(1, -0.1, 0.1)
    #gaussFit.SetParLimits(2, 0, 0.1)
    hist_DUT_Res.Fit(gaussFit)

    #polyFit = ROOT.TF1("polyFit", "pol10", -0.5, 0.5)
    #hist_DUT_Res.Fit(polyFit, "S")

    c1.Update()
    c1.SaveAs("DUT-Res-eLOC0-wot.pdf")
    c1.SaveAs("DUT-Res-eLOC0-wot.png")

file.Close()

input()
