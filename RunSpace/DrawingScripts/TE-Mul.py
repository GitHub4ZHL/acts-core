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
canvas = ROOT.TCanvas("canvas", "", 1000, 800)
canvas.Divide(1, 2)

# Efficiency
te_wot = []
te_30ps = []
te_50ps = []
te_200ps = []

# Error
teErr_wot = []
teErr_30ps = []
teErr_50ps = []
teErr_200ps = []

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
    te_wot_value = eff.GetEfficiency(1)
    teErr_wot_value = (eff.GetEfficiencyErrorLow(1)+eff.GetEfficiencyErrorUp(1))/2
    te_wot.append(te_wot_value)
    teErr_wot.append(teErr_wot_value)
    file.Close()

for file_path in file_paths_wt_30ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    te_30ps_value = eff.GetEfficiency(1)
    teErr_30ps_value = (eff.GetEfficiencyErrorLow(1)+eff.GetEfficiencyErrorUp(1))/2
    te_30ps.append(te_30ps_value)
    teErr_30ps.append(teErr_30ps_value)
    file.Close()

for file_path in file_paths_wt_50ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    te_50ps_value = eff.GetEfficiency(1)
    teErr_50ps_value = (eff.GetEfficiencyErrorLow(1)+eff.GetEfficiencyErrorUp(1))/2
    te_50ps.append(te_50ps_value)
    teErr_50ps.append(teErr_50ps_value)
    file.Close()

for file_path in file_paths_wt_200ps:
    file = ROOT.TFile(file_path, "READ")
    eff_name = "perfSummary"
    eff = file.Get(eff_name)
    te_200ps_value = eff.GetEfficiency(1)
    teErr_200ps_value = (eff.GetEfficiencyErrorLow(1)+eff.GetEfficiencyErrorUp(1))/2
    te_200ps.append(te_200ps_value)
    teErr_200ps.append(teErr_200ps_value)
    file.Close()

# Get ratio value
for i in range(0, 20):
    rat_30ps_value = te_30ps[i]/te_wot[i]
    rat_30ps.append(rat_30ps_value)
    ratErr_30ps_value = rat_30ps[i]*math.hypot(teErr_wot[i]/te_wot[i], teErr_30ps[i]/te_30ps[i])
    ratErr_30ps.append(ratErr_30ps_value)

    rat_50ps_value = te_50ps[i]/te_wot[i]
    rat_50ps.append(rat_50ps_value)
    ratErr_50ps_value = rat_50ps[i]*math.hypot(teErr_wot[i]/te_wot[i], teErr_50ps[i]/te_50ps[i])
    ratErr_50ps.append(ratErr_50ps_value)

    rat_200ps_value = te_200ps[i]/te_wot[i]
    rat_200ps.append(rat_200ps_value)
    ratErr_200ps_value = rat_200ps[i]*math.hypot(teErr_wot[i]/te_wot[i], teErr_200ps[i]/te_200ps[i])
    ratErr_200ps.append(ratErr_200ps_value)

# Draw efficiency
eff_pad = canvas.cd(1)
eff_pad.SetTopMargin(0.1)
eff_pad.SetBottomMargin(0.1)

n_points = len(multiplicity_range)

# wot
te_wot_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, te_wot)):
    te_wot_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(teErr_wot)):
    te_wot_graph.SetPointError(i, 0, teErr_wot[i-1])

te_wot_graph.SetMarkerStyle(ROOT.kFullCircle)
te_wot_graph.SetMarkerSize(0.9)
te_wot_graph.SetMarkerColor(ROOT.kBlue)
te_wot_graph.GetXaxis().SetRangeUser(0, 21)
te_wot_graph.GetYaxis().SetRangeUser(0.955, 1.02)
te_wot_graph.SetLineStyle(1)
te_wot_graph.SetLineColor(ROOT.kBlue)
te_wot_graph.SetTitle("Efficiency vs. Multiplicity")
# Axis Title
te_wot_graph.GetXaxis().SetTitle("Multiplicity")
te_wot_graph.GetYaxis().SetTitle("Tracking Efficiency")
te_wot_graph.Draw("AP")  # "A" as Axis，"P" as Plot

# 30ps
te_30ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, te_30ps)):
    te_30ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(teErr_30ps)):
    te_30ps_graph.SetPointError(i, 0, teErr_30ps[i-1])

te_30ps_graph.SetMarkerStyle(ROOT.kOpenCircle)
te_30ps_graph.SetMarkerSize(0.9)
te_30ps_graph.SetMarkerColor(ROOT.kRed)
te_30ps_graph.SetLineStyle(1)
te_30ps_graph.SetLineColor(ROOT.kRed)
te_30ps_graph.Draw("P same")

# 50ps
te_50ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, te_50ps)):
    te_50ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(teErr_50ps)):
    te_50ps_graph.SetPointError(i, 0, teErr_30ps[i-1])

te_50ps_graph.SetMarkerStyle(ROOT.kOpenSquare)
te_50ps_graph.SetMarkerSize(0.9)
te_50ps_graph.SetMarkerColor(ROOT.kGreen)
te_50ps_graph.SetLineStyle(1)
te_50ps_graph.SetLineColor(ROOT.kGreen)
te_50ps_graph.Draw("P same")

# 200ps
te_200ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, te_200ps)):
    te_200ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(teErr_200ps)):
    te_200ps_graph.SetPointError(i, 0, teErr_200ps[i-1])

te_200ps_graph.SetMarkerStyle(ROOT.kOpenTriangleUp)
te_200ps_graph.SetMarkerSize(0.9)
te_200ps_graph.SetMarkerColor(ROOT.kViolet)
te_200ps_graph.SetLineStyle(1)
te_200ps_graph.SetLineColor(ROOT.kViolet)
te_200ps_graph.Draw("P same")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(te_wot_graph, "wot", "pl")
legend.AddEntry(te_30ps_graph, "30ps", "pl")
legend.AddEntry(te_50ps_graph, "50ps", "pl")
legend.AddEntry(te_200ps_graph, "200ps", "pl")
legend.Draw()



# Draw ratio
rat_pad=canvas.cd(2)
rat_pad.SetTopMargin(0)
rat_pad.SetBottomMargin(0.7)
# Draw ratio for 30ps
n_points = len(multiplicity_range)
rat_30ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_30ps)):
    rat_30ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_30ps)):
    rat_30ps_graph.SetPointError(i, 0, ratErr_30ps[i-1])

rat_30ps_graph.SetMarkerStyle(ROOT.kOpenCircle)
rat_30ps_graph.SetMarkerSize(0.9)
rat_30ps_graph.SetMarkerColor(ROOT.kRed)
rat_30ps_graph.SetLineStyle(1)
rat_30ps_graph.SetLineColor(ROOT.kRed)
rat_30ps_graph.GetXaxis().SetTitle("Multiplicity")
rat_30ps_graph.GetYaxis().SetTitle("W/Wo")
rat_30ps_graph.Draw("AP")
rat_30ps_graph.GetXaxis().SetRangeUser(0, 21)
rat_30ps_graph.GetYaxis().SetRangeUser(0.99, 1.042)
rat_30ps_graph.SetTitle("")

# Draw ratio for 50ps
rat_50ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_50ps)):
    rat_50ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_50ps)):
    rat_50ps_graph.SetPointError(i, 0, ratErr_50ps[i-1])

rat_50ps_graph.SetMarkerStyle(ROOT.kOpenSquare)
rat_50ps_graph.SetMarkerSize(0.9)
rat_50ps_graph.SetMarkerColor(ROOT.kGreen)
rat_50ps_graph.SetLineStyle(1)
rat_50ps_graph.SetLineColor(ROOT.kGreen)
rat_50ps_graph.Draw("P same")

# Draw ratio for 200ps
rat_200ps_graph = ROOT.TGraphErrors(n_points)
for i, (x, y) in enumerate(zip(multiplicity_range, rat_200ps)):
    rat_200ps_graph.SetPoint(i, x, y)
for i, (e) in enumerate(zip(ratErr_200ps)):
    rat_200ps_graph.SetPointError(i, 0, ratErr_200ps[i-1])

rat_200ps_graph.SetMarkerStyle(ROOT.kOpenTriangleUp)
rat_200ps_graph.SetMarkerSize(0.9)
rat_200ps_graph.SetMarkerColor(ROOT.kViolet)
rat_200ps_graph.SetLineStyle(1)
rat_200ps_graph.SetLineColor(ROOT.kViolet)
rat_200ps_graph.Draw("P same")

# Dotted line
dotted_line = ROOT.TLine(rat_30ps_graph.GetXaxis().GetXmin(), 1, 21, 1)
dotted_line.SetLineStyle(2)
dotted_line.Draw("same")

# Draw canvas
canvas.SetTitle("Efficiency vs. Multiplicity")
canvas.Update()
canvas.Draw()
canvas.SaveAs("te-mul.png")
#canvas.SaveAs("te-mul.pdf")
input()
