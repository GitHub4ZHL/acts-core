import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_chi = 100
chi2_range = range(10, max_chi+5, 5) #(begin number, final number + 1, step length)
file_paths_1 = [f"{base_path}wt_1_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]
file_paths_15 = [f"{base_path}wt_15_chi2_{chi2}/performance_ckf.root" for chi2 in chi2_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

te_1 = []
te_15 = []
TEErr_1 = []
TEErr_15 = []

# Build WO-Time TGraph

for file_path in file_paths_1:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    print(value)
    teerr_1 = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_1.append(value)
    TEErr_1.append(teerr_1)
    file.Close()

n_points = len(chi2_range)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_1)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_1)):
    graph1.SetPointError(i, 0, TEErr_1[i-1])

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

for file_path in file_paths_15:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    print(value)
    teerr_15 = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_15.append(value)
    TEErr_15.append(teerr_15)
    file.Close()
n_points = len(chi2_range)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(chi2_range, te_15)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_15)):
    graph2.SetPointError(i, 0, TEErr_15[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 110)
graph2.GetYaxis().SetRangeUser(0.94, 1.02)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Chi2")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "multiplicity=1", "p")
legend.AddEntry(graph2, "multiplicity=15", "p")
legend.Draw()

canvas.Update()
canvas.Draw()

# Title
canvas.SetTitle("Tracking Efficiency Compare about Chi2")
canvas.SaveAs("TE-Chi2-wt.pdf")
canvas.SaveAs("TE-Chi2-wt.png")
