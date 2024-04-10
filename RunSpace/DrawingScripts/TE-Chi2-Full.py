import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_chi = 100
chi2_range = range(10, max_chi+5, 5) #(begin number, final number + 1, step length)
file_paths_1_wot = [f"{base_path}wot_1_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]
file_paths_15_wot = [f"{base_path}wot_15_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]
file_paths_1_wt = [f"{base_path}wt_1_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]
file_paths_15_wt = [f"{base_path}wt_15_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

te_1_wot = []
te_15_wot = []
TEErr_1_wot = []
TEErr_15_wot = []

te_1_wt = []
te_15_wt = []
TEErr_1_wt = []
TEErr_15_wt = []

# Build WO-Time TGraph

for file_path in file_paths_1_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    #print(value)
    teerr_1_wot = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_1_wot.append(value)
    TEErr_1_wot.append(teerr_1_wot)
    file.Close()

n_points = len(chi2_range)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_1_wot)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_1_wot)):
    graph1.SetPointError(i, 0, TEErr_1_wot[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 110)
graph1.GetYaxis().SetRangeUser(0.94, 1.02)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Tracking Efficiency Compare about Chi2")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_15_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    #print(value)
    teerr_15_wot = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_15_wot.append(value)
    TEErr_15_wot.append(teerr_15_wot)
    file.Close()
n_points = len(chi2_range)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_15_wot)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_15_wot)):
    graph2.SetPointError(i, 0, TEErr_15_wot[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 110)
graph2.GetYaxis().SetRangeUser(0.94, 1.02)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

for file_path in file_paths_1_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    #print(value)
    teerr_1_wt = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_1_wt.append(value)
    TEErr_1_wt.append(teerr_1_wt)
    file.Close()

n_points = len(chi2_range)

graph3 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_1_wt)):
    graph3.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_1_wt)):
    graph3.SetPointError(i, 0, TEErr_1_wt[i-1])

graph3.SetMarkerStyle(ROOT.kFullCircle)
graph3.SetMarkerSize(0.6)
graph3.SetMarkerColor(ROOT.kGreen)
graph3.GetXaxis().SetRangeUser(0, 110)
graph3.GetYaxis().SetRangeUser(0.94, 1.02)
graph3.SetLineStyle(1)
graph3.SetLineColor(ROOT.kGreen)
graph3.SetTitle("Tracking Efficiency Compare about Chi2")
graph3.Draw("PL same")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_15_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    #print(value)
    teerr_15_wt = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_15_wt.append(value)
    TEErr_15_wt.append(teerr_15_wt)
    file.Close()
n_points = len(chi2_range)

graph4 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_15_wt)):
    graph4.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_15_wt)):
    graph4.SetPointError(i, 0, TEErr_15_wt[i-1])

graph4.SetMarkerStyle(ROOT.kFullSquare)
graph4.SetMarkerSize(0.6)
graph4.SetMarkerColor(ROOT.kYellow)
graph4.GetXaxis().SetRangeUser(0, 110)
graph4.GetYaxis().SetRangeUser(0.94, 1.02)
graph4.SetLineStyle(1)
graph4.SetLineColor(ROOT.kYellow)
graph4.Draw("PL same")


# Axis Title
graph1.GetXaxis().SetTitle("Chi2")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "mul=1 wot", "p")
legend.AddEntry(graph2, "mul=15 wot", "p")
legend.AddEntry(graph3, "mul=1 wt", "p")
legend.AddEntry(graph4, "mul=15 wt", "p")

legend.Draw()

canvas.Update()
canvas.Draw()

# Title
canvas.SetTitle("Tracking Efficiency Compare about Chi2")
canvas.SaveAs("TE-Chi2-4.pdf")
canvas.SaveAs("TE-Chi2-4.png")
