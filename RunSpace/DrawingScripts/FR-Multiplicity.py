import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_mul = 20
multiplicity_range = range(1, max_mul+1, 1) #(begin number, final number + 1, step length)
file_paths_wot = [f"{base_path}wot_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]
file_paths_wt = [f"{base_path}wt_multiplicity_{mult}/performance_ckf.root" for mult in multiplicity_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

fr_wotime = []
fr_wtime = []
FRErr_WoT = []
FRErr_WT = []

# Build WO-Time TGraph

for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    print(value)
    frerr_wot = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_wotime.append(value)
    FRErr_WoT.append(frerr_wot)
    file.Close()

n_points = len(multiplicity_range)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity_range, fr_wotime)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_WoT)):
    graph1.SetPointError(i, 0, FRErr_WoT[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 21)
graph1.GetYaxis().SetRangeUser(0, 0.4)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Fake Rate Compare about Multiplicity")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    print(value)
    frerr_wt = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_wtime.append(value)
    FRErr_WT.append(frerr_wt)
    file.Close()
n_points = len(multiplicity_range)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity_range, fr_wtime)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_WT)):
    graph2.SetPointError(i, 0, FRErr_WT[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 21)
graph2.GetYaxis().SetRangeUser(0, 0.4)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Multiplicity")
graph1.GetYaxis().SetTitle("Fake Rate")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

canvas.Update()
canvas.Draw()
# Title
canvas.SetTitle("Fake Rate Compare about Multiplicity")

canvas.SaveAs("FR-Multiplicity.pdf")

