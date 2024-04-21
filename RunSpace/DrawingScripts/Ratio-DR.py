import ROOT
import math

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_mul = 20
multiplicity_range = range(1, max_mul+1, 1) #(begin number, final number + 1, step length)
file_paths_wot = [f"{base_path}wot_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_30ps = [f"{base_path}wt_multiplicity_{mult}_30ps/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_50ps = [f"{base_path}wt_multiplicity_{mult}_50ps/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt_200ps = [f"{base_path}wt_multiplicity_{mult}_200ps/performance_ckf.root" for mult in multiplicity_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Ratio", 800, 600)

# Efficiency
dr_wot = []
dr_wt_30ps = [] 
dr_wt_50ps = []
dr_wt_200ps = []

# Error
drErr_wot = []
drErr_wt_30ps = []
drErr_wt_50ps = []
drErr_wt_200ps = []

# Ratio
rat_30ps = []
rat_50ps = []
rat_200ps = []

# Ratio Error
ratErr_30ps = []
ratErr_50ps = []
ratErr_200ps = []

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
    dr_wt_30ps_value = eff.GetEfficiency(3)
    drErr_wt_30ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_wt_30ps.append(dr_wt_30ps_value)
    drErr_wt_30ps.append(drErr_wt_30ps_value)
    file.Close()

for file_path in file_paths_wt_50ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_wt_50ps_value = eff.GetEfficiency(3)
    drErr_wt_50ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_wt_50ps.append(dr_wt_50ps_value)
    drErr_wt_50ps.append(drErr_wt_50ps_value)
    file.Close()

for file_path in file_paths_wt_200ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    dr_wt_200ps_value = eff.GetEfficiency(3)
    drErr_wt_200ps_value = (eff.GetEfficiencyErrorLow(3)+eff.GetEfficiencyErrorUp(3))/2
    dr_wt_200ps.append(dr_wt_200ps_value)
    drErr_wt_200ps.append(drErr_wt_200ps_value)
    file.Close()

# Get ratio value
for i in range(0, 20):
  rat_30ps_value = dr_wt_30ps[i]/dr_wot[i]
  rat_30ps.append(rat_30ps_value)
  ratErr_30ps_value = rat_30ps[i]*math.hypot(drErr_wot[i]/dr_wot[i], drErr_wt_30ps[i]/dr_wt_30ps[i])
  ratErr_30ps.append(ratErr_30ps_value)

  rat_50ps_value = dr_wt_50ps[i]/dr_wot[i]
  rat_50ps.append(rat_50ps_value)
  ratErr_50ps_value = rat_50ps[i]*math.hypot(drErr_wot[i]/dr_wot[i], drErr_wt_50ps[i]/dr_wt_50ps[i])
  ratErr_50ps.append(ratErr_50ps_value)

  rat_200ps_value = dr_wt_200ps[i]/dr_wot[i]
  rat_200ps.append(rat_200ps_value)
  ratErr_200ps_value = rat_200ps[i]*math.hypot(drErr_wot[i]/dr_wot[i], drErr_wt_200ps[i]/dr_wt_200ps[i])
  ratErr_200ps.append(ratErr_200ps_value)

# Draw ratio for 30ps
n_points = len(multiplicity_range)
rat_30ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_30ps)):
    rat_30ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_30ps)):
    rat_30ps_graph.SetPointError(i, 0, ratErr_30ps[i-1])

rat_30ps_graph.SetMarkerStyle(ROOT.kFullSquare)
rat_30ps_graph.SetMarkerSize(0.8)
rat_30ps_graph.SetMarkerColor(ROOT.kRed)
rat_30ps_graph.SetLineStyle(1)
rat_30ps_graph.SetLineColor(ROOT.kRed)
rat_30ps_graph.Draw("AP")
rat_30ps_graph.GetXaxis().SetRangeUser(0, 21)
rat_30ps_graph.GetYaxis().SetRangeUser(0.4, 1.3)
rat_30ps_graph.SetTitle("Duplicated Rate Ratio vs. Multiplicity")

# Draw ratio for 50ps
n_points = len(multiplicity_range)
rat_50ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_50ps)):
    rat_50ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_50ps)):
    rat_50ps_graph.SetPointError(i, 0, ratErr_50ps[i-1])

rat_50ps_graph.SetMarkerStyle(ROOT.kFullSquare)
rat_50ps_graph.SetMarkerSize(0.8)
rat_50ps_graph.SetMarkerColor(ROOT.kBlue)
rat_50ps_graph.SetLineStyle(1)
rat_50ps_graph.SetLineColor(ROOT.kBlue)
rat_50ps_graph.Draw("P same")

# Draw ratio for 200ps
n_points = len(multiplicity_range)
rat_200ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_200ps)):
    rat_200ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_200ps)):
    rat_200ps_graph.SetPointError(i, 0, ratErr_200ps[i-1])

rat_200ps_graph.SetMarkerStyle(ROOT.kFullSquare)
rat_200ps_graph.SetMarkerSize(0.8)
rat_200ps_graph.SetMarkerColor(ROOT.kGreen)
rat_200ps_graph.SetLineStyle(1)
rat_200ps_graph.SetLineColor(ROOT.kGreen)
rat_200ps_graph.Draw("P same")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(rat_30ps_graph, "30ps", "pl")
legend.AddEntry(rat_50ps_graph, "50ps", "pl")
legend.AddEntry(rat_200ps_graph, "200ps", "pl")
legend.Draw()

# Axis Title
rat_30ps_graph.GetXaxis().SetTitle("Multiplicity")
rat_30ps_graph.GetYaxis().SetTitle("Efficiency Ratio")

# Draw canvas
canvas.SetTitle("Duplicated Rate Ratio vs. Multiplicity")
canvas.Update()
canvas.Draw()
#canvas.SaveAs("ratio-dr.png")
canvas.SaveAs("ratio-dr.pdf")
#input()
