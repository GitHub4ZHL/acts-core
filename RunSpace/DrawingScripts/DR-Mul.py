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
dr_wot = []
dr_30ps = []
dr_50ps = []
dr_200ps = []

# Error
drErr_wot = []
drErr_30ps = []
drErr_50ps = []
drErr_200ps = []

# Get efficiency and error
for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_wot_value = eff.GetEfficiency(3)
    drErr_wot_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_wot.append(dr_wot_value)
    drErr_wot.append(drErr_wot_value)
    file.Close()

for file_path in file_paths_wt_30ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_30ps_value = eff.GetEfficiency(3)
    drErr_30ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_30ps.append(dr_30ps_value)
    drErr_30ps.append(drErr_30ps_value)
    file.Close()

for file_path in file_paths_wt_50ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_50ps_value = eff.GetEfficiency(3)
    drErr_50ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_50ps.append(dr_50ps_value)
    drErr_50ps.append(drErr_50ps_value)
    file.Close()

for file_path in file_paths_wt_200ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_200ps_value = eff.GetEfficiency(3)
    drErr_200ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_200ps.append(dr_200ps_value)
    drErr_200ps.append(drErr_200ps_value)
    file.Close()

# Draw efficiency
n_points = len(multiplicity_range)

# wot
dr_wot_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, dr_wot)):
    dr_wot_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(drErr_wot)):
    dr_wot_graph.SetPointError(i, 0, drErr_wot[i-1])

dr_wot_graph.SetMarkerStyle(ROOT.kOpenCircle)
dr_wot_graph.SetMarkerSize(0.8)
dr_wot_graph.SetMarkerColor(ROOT.kBlue)
dr_wot_graph.GetXaxis().SetRangeUser(0, 21)
dr_wot_graph.GetYaxis().SetRangeUser(0, 0.018)
dr_wot_graph.SetLineStyle(1)
dr_wot_graph.SetLineColor(ROOT.kBlue)
dr_wot_graph.SetTitle("Duplicated Rate vs. Multiplicity")
dr_wot_graph.Draw("AP")  # "A" as Axis，"P" as Plot

# 30ps
dr_30ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, dr_30ps)):
    dr_30ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(drErr_30ps)):
    dr_30ps_graph.SetPointError(i, 0, drErr_30ps[i-1])

dr_30ps_graph.SetMarkerStyle(ROOT.kFullCircle)
dr_30ps_graph.SetMarkerSize(0.8)
dr_30ps_graph.SetMarkerColor(ROOT.kRed)
dr_30ps_graph.SetLineStyle(1)
dr_30ps_graph.SetLineColor(ROOT.kRed)
dr_30ps_graph.Draw("P same")

# 50ps
dr_50ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, dr_50ps)):
    dr_50ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(drErr_50ps)):
    dr_50ps_graph.SetPointError(i, 0, drErr_30ps[i-1])

dr_50ps_graph.SetMarkerStyle(ROOT.kFullSquare)
dr_50ps_graph.SetMarkerSize(0.8)
dr_50ps_graph.SetMarkerColor(ROOT.kGreen)
dr_50ps_graph.SetLineStyle(1)
dr_50ps_graph.SetLineColor(ROOT.kGreen)
dr_50ps_graph.Draw("P same")

# 200ps
dr_200ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, dr_200ps)):
    dr_200ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(drErr_200ps)):
    dr_200ps_graph.SetPointError(i, 0, drErr_200ps[i-1])

dr_200ps_graph.SetMarkerStyle(ROOT.kFullTriangleUp)
dr_200ps_graph.SetMarkerSize(0.8)
dr_200ps_graph.SetMarkerColor(ROOT.kViolet)
dr_200ps_graph.SetLineStyle(1)
dr_200ps_graph.SetLineColor(ROOT.kViolet)
dr_200ps_graph.Draw("P same")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(dr_wot_graph, "wot", "pl")
legend.AddEntry(dr_30ps_graph, "30ps", "pl")
legend.AddEntry(dr_50ps_graph, "50ps", "pl")
legend.AddEntry(dr_200ps_graph, "200ps", "pl")
legend.Draw()

# Axis Title
dr_30ps_graph.GetXaxis().SetTitle("Multiplicity")
dr_30ps_graph.GetYaxis().SetTitle("Duplicated Rate")

# Draw canvas
canvas.SetTitle("Duplicated Rate vs. Multiplicity")
canvas.Update()
canvas.Draw()
canvas.SaveAs("dr-mul.png")
#canvas.SaveAs("dr-mul.pdf")
input()

