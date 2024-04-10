import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_mul = 20
multiplicity_range = range(1, max_mul+1, 1) #(begin number, final number + 1, step length)
file_paths_wot = [f"{base_path}wot_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt = [f"{base_path}wt_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]


# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

multiplicity = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
dr_wotime = []
dr_wtime = []
DRErr_WoT = []
DRErr_WT = []

# Build WO-Time TGraph

for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(3)
    print(value)
    drerr_wot = (tefficiency.GetEfficiencyErrorLow(3)+tefficiency.GetEfficiencyErrorUp(3))/2
    dr_wotime.append(value)
    DRErr_WoT.append(drerr_wot)
    file.Close()

n_points = len(multiplicity_range)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity_range, dr_wotime)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(DRErr_WoT)):
    graph1.SetPointError(i, 0, DRErr_WoT[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 21)
graph1.GetYaxis().SetRangeUser(0, 0.025)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Duplicated Rate Compare about Multiplicity")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
for file_path in file_paths_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(3)
    print(value)
    drerr_wt = (tefficiency.GetEfficiencyErrorLow(3)+tefficiency.GetEfficiencyErrorUp(3))/2
    dr_wtime.append(value)
    DRErr_WT.append(drerr_wt)
    file.Close()
n_points = len(multiplicity_range)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity_range, dr_wtime)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(DRErr_WT)):
    graph2.SetPointError(i, 0, DRErr_WT[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 21)
graph2.GetYaxis().SetRangeUser(0, 0.025)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Multiplicity")
graph1.GetYaxis().SetTitle("Duplicated Rate")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

canvas.Update()
canvas.Draw()
# Title
canvas.SetTitle("Duplicated Rate Compare about Multiplicity")

canvas.SaveAs("DR-Multiplicity.pdf")
canvas.SaveAs("DR-Multiplicity.png")

