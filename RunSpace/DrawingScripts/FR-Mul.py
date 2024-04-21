import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_mul = 20
multiplicity_range = range(1, max_mul+1, 1) #(begin number, final number + 1, step length)
file_paths_wot = [f"{base_path}wot_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_30ps = [f"{base_path}wt_multiplicity_{mult}_30ps/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_50ps = [f"{base_path}wt_multiplicity_{mult}_50ps/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_200ps = [f"{base_path}wt_multiplicity_{mult}_200ps/performance_ckf.root" for mult in multiplicity_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

# Efficiency
fr_wot = []
fr_30ps = []
fr_50ps = []
fr_200ps = []

# Error
frErr_wot = []
frErr_30ps = []
frErr_50ps = []
frErr_200ps = []

# Get efficiency and error
for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    fr_wot_value = eff.GetEfficiency(2)
    frErr_wot_value = (eff.GetEfficiencyErrorLow(2)+eff.GetEfficiencyErrorUp(2))/2
    fr_wot.append(fr_wot_value)
    frErr_wot.append(frErr_wot_value)
    file.Close()

for file_path in file_paths_wt_30ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    fr_30ps_value = eff.GetEfficiency(2)
    frErr_30ps_value = (eff.GetEfficiencyErrorLow(2)+eff.GetEfficiencyErrorUp(2))/2
    fr_30ps.append(fr_30ps_value)
    frErr_30ps.append(frErr_30ps_value)
    file.Close()

for file_path in file_paths_wt_50ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    fr_50ps_value = eff.GetEfficiency(2)
    frErr_50ps_value = (eff.GetEfficiencyErrorLow(2)+eff.GetEfficiencyErrorUp(2))/2
    fr_50ps.append(fr_50ps_value)
    frErr_50ps.append(frErr_50ps_value)
    file.Close()

for file_path in file_paths_wt_200ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    fr_200ps_value = eff.GetEfficiency(2)
    frErr_200ps_value = (eff.GetEfficiencyErrorLow(2)+eff.GetEfficiencyErrorUp(2))/2
    fr_200ps.append(fr_200ps_value)
    frErr_200ps.append(frErr_200ps_value)
    file.Close()

# Draw efficiency
n_points = len(multiplicity_range)

# wot
fr_wot_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, fr_wot)):
    fr_wot_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(frErr_wot)):
    fr_wot_graph.SetPointError(i, 0, frErr_wot[i-1])

fr_wot_graph.SetMarkerStyle(ROOT.kOpenCircle)
fr_wot_graph.SetMarkerSize(0.8)
fr_wot_graph.SetMarkerColor(ROOT.kBlue)
fr_wot_graph.GetXaxis().SetRangeUser(0, 21)
fr_wot_graph.GetYaxis().SetRangeUser(0, 0.002)
fr_wot_graph.SetLineStyle(1)
fr_wot_graph.SetLineColor(ROOT.kBlue)
fr_wot_graph.SetTitle("Fake Rate vs. Multiplicity")
fr_wot_graph.Draw("AP")  # "A" as Axis，"P" as Plot

# 30ps
fr_30ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, fr_30ps)):
    fr_30ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(frErr_30ps)):
    fr_30ps_graph.SetPointError(i, 0, frErr_30ps[i-1])

fr_30ps_graph.SetMarkerStyle(ROOT.kFullCircle)
fr_30ps_graph.SetMarkerSize(0.8)
fr_30ps_graph.SetMarkerColor(ROOT.kRed)
fr_30ps_graph.SetLineStyle(1)
fr_30ps_graph.SetLineColor(ROOT.kRed)
fr_30ps_graph.Draw("P same")

# 50ps
fr_50ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, fr_50ps)):
    fr_50ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(frErr_50ps)):
    fr_50ps_graph.SetPointError(i, 0, frErr_30ps[i-1])

fr_50ps_graph.SetMarkerStyle(ROOT.kFullSquare)
fr_50ps_graph.SetMarkerSize(0.8)
fr_50ps_graph.SetMarkerColor(ROOT.kGreen)
fr_50ps_graph.SetLineStyle(1)
fr_50ps_graph.SetLineColor(ROOT.kGreen)
fr_50ps_graph.Draw("P same")

# 200ps
fr_200ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, fr_200ps)):
    fr_200ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(frErr_200ps)):
    fr_200ps_graph.SetPointError(i, 0, frErr_200ps[i-1])

fr_200ps_graph.SetMarkerStyle(ROOT.kFullTriangleUp)
fr_200ps_graph.SetMarkerSize(0.8)
fr_200ps_graph.SetMarkerColor(ROOT.kViolet)
fr_200ps_graph.SetLineStyle(1)
fr_200ps_graph.SetLineColor(ROOT.kViolet)
fr_200ps_graph.Draw("P same")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(fr_wot_graph, "wot", "pl")
legend.AddEntry(fr_30ps_graph, "30ps", "pl")
legend.AddEntry(fr_50ps_graph, "50ps", "pl")
legend.AddEntry(fr_200ps_graph, "200ps", "pl")
legend.Draw()

# Axis Title
fr_30ps_graph.GetXaxis().SetTitle("Multiplicity")
fr_30ps_graph.GetYaxis().SetTitle("Fake Rate")

# Draw canvas
canvas.SetTitle("Fake Rate vs. Multiplicity")
canvas.Update()
canvas.Draw()
canvas.SaveAs("fr-mul.png")
#canvas.SaveAs("fr-mul.pdf")
input()

