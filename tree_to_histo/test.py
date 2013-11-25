import ROOT
h = ROOT.TH1F("test","test",10,-2.5,2.5)
h.Fill(-2.5)
h.Fill(0)
h.Fill(1)
outfile = ROOT.TFile("test.root", "recreate")

h.Write()

outfile.Close()
