import ROOT

# File Path
base_path = "/home/haolin/HEP-Software/ACTS/acts/RunSpace/"
max_nmea = 10
nmea_range = range(1, max_nmea+1, 1) #(begin number, final number + 1, step length)
file_paths_1_wot = [f"{base_path}wot_1_nmea_{nmea}/performance_ckf.root" for nmea in nmea_range]
file_paths_15_wot = [f"{base_path}wot_15_nmea_{nmea}/performance_ckf.root" for nmea in nmea_range]
file_paths_1_wt = [f"{base_path}wt_1_nmea_{nmea}/performance_ckf.root" for nmea in nmea_range]
file_paths_15_wt = [f"{base_path}wt_15_nmea_{nmea}/performance_ckf.root" for nmea in nmea_range]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

fr_1_wot = []
fr_15_wot = []
FRErr_1_wot = []
FRErr_15_wot = []

fr_1_wt = []
fr_15_wt = []
FRErr_1_wt = []
FRErr_15_wt = []

# Build WO-Time TGraph

for file_path in file_paths_1_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    print(value)
    frerr_1_wot = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_1_wot.append(value)
    FRErr_1_wot.append(frerr_1_wot)
    file.Close()

n_points = len(nmea_range)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(nmea_range, fr_1_wot)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_1_wot)):
    graph1.SetPointError(i, 0, FRErr_1_wot[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 110)
graph1.GetYaxis().SetRangeUser(0, 0.02)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Fake Rate Compare about numMeasurementsCut")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_15_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    #print(value)
    frerr_15_wot = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_15_wot.append(value)
    FRErr_15_wot.append(frerr_15_wot)
    file.Close()
n_points = len(nmea_range)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(nmea_range, fr_15_wot)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_15_wot)):
    graph2.SetPointError(i, 0, FRErr_15_wot[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 110)
graph2.GetYaxis().SetRangeUser(0, 0.02)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

for file_path in file_paths_1_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    #print(value)
    frerr_1_wt = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_1_wt.append(value)
    FRErr_1_wt.append(frerr_1_wt)
    file.Close()

n_points = len(nmea_range)

graph3 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(nmea_range, fr_1_wt)):
    graph3.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_1_wt)):
    graph3.SetPointError(i, 0, FRErr_1_wt[i-1])

graph3.SetMarkerStyle(ROOT.kFullCircle)
graph3.SetMarkerSize(0.6)
graph3.SetMarkerColor(ROOT.kGreen)
graph3.GetXaxis().SetRangeUser(0, 110)
graph3.GetYaxis().SetRangeUser(0, 0.02)
graph3.SetLineStyle(1)
graph3.SetLineColor(ROOT.kGreen)
graph3.SetTitle("Fake Rate Compare about numMeasurementsCut")
graph3.Draw("PL same")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_15_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    #print(value)
    frerr_15_wt = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    fr_15_wt.append(value)
    FRErr_15_wt.append(frerr_15_wt)
    file.Close()
n_points = len(nmea_range)

graph4 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(nmea_range, fr_15_wt)):
    graph4.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_15_wt)):
    graph4.SetPointError(i, 0, FRErr_15_wt[i-1])

graph4.SetMarkerStyle(ROOT.kFullSquare)
graph4.SetMarkerSize(0.6)
graph4.SetMarkerColor(ROOT.kViolet)
graph4.GetXaxis().SetRangeUser(0, 110)
graph4.GetYaxis().SetRangeUser(0, 0.02)
graph4.SetLineStyle(1)
graph4.SetLineColor(ROOT.kViolet)
graph4.Draw("PL same")


# Axis Title
graph1.GetXaxis().SetTitle("numMeasurementsCut")
graph1.GetYaxis().SetTitle("Fake Rate")

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
canvas.SetTitle("Fake Rate Compare about numMeasurementsCut")
canvas.SaveAs("FR-nMea-4.pdf")
#canvas.SaveAs("FR-nMea-4.png")
