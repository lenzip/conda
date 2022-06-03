import ROOT

def getWorkspace(ww, top, dytt, vv, sig, label):

  # get x axis histogram from data (assume it is the same for all)
  # x is a variable running over the bins that have been measured
  x = ROOT.RooRealVar("x", "x", ww.GetXaxis().GetXmin(),  ww.GetXaxis().GetXmax())

  # transform the normalizations of each background in a RooRealVar, fixed
  nu_ww   = ROOT.RooRealVar("nu_ww",   "nu_ww",   ww.Integral())
  nu_top  = ROOT.RooRealVar("nu_top",  "nu_top",  top.Integral())
  nu_dytt = ROOT.RooRealVar("nu_dytt", "nu_dytt", dytt.Integral())
  nu_vv   = ROOT.RooRealVar("nu_vv",   "nu_vv",   vv.Integral())
  # same for the signal
  nu_s    = ROOT.RooRealVar("nu_s",    "nu_s",    sig.Integral())

  # Define a multiplier for each normalization, including the signal
  mu = ROOT.RooRealVar("mu",      "mu",      0., -10., 10.) # this is the POI
  # the following are nuisances
  mu_ww   = ROOT.RooRealVar("mu_ww",   "mu_ww",   1, 0.01, 5)
  mu_top  = ROOT.RooRealVar("mu_top",  "mu_top",  1, 0.01, 5) 
  mu_dytt = ROOT.RooRealVar("mu_dytt", "mu_dytt", 1, 0.01, 5) 
  mu_vv   = ROOT.RooRealVar("mu_vv",   "mu_vv",   1, 0.01, 5) 

  # we now write the normalization of each background and the signal as
  # the relevant mu_ times teh relevant nu_
  norm_s    = ROOT.RooFormulaVar("norm_s",    "mu*nu_s",         ROOT.RooArgSet(mu, nu_s))
  norm_ww   = ROOT.RooFormulaVar("norm_ww",   "mu_ww*nu_ww",     ROOT.RooArgSet(mu_ww, nu_ww))
  norm_top  = ROOT.RooFormulaVar("norm_top",  "mu_top*nu_top",   ROOT.RooArgSet(mu_top, nu_top))
  norm_dytt = ROOT.RooFormulaVar("norm_dytt", "mu_dytt*nu_dytt", ROOT.RooArgSet(mu_dytt, nu_dytt))
  norm_vv   = ROOT.RooFormulaVar("norm_vv",   "mu_vv*nu_vv",     ROOT.RooArgSet(mu_vv, nu_vv))

  # Define PDFs in x for each background and for the signal:
  # to do so we need to transform each input histogram
  # to a RooDataHist and then to a RooHistPdf
  h_s = ROOT.RooDataHist("h_s",    "h_s",    ROOT.RooArgSet(x), sig)
  h_ww   = ROOT.RooDataHist("h_ww",   "h_ww",   ROOT.RooArgSet(x), ww)
  h_top  = ROOT.RooDataHist("h_top",  "h_top",  ROOT.RooArgSet(x), top)
  h_dytt = ROOT.RooDataHist("h_dytt", "h_dytt", ROOT.RooArgSet(x), dytt)
  h_vv   = ROOT.RooDataHist("h_vv",   "h_vv",   ROOT.RooArgSet(x), vv)

  pdf_s    =  ROOT.RooHistPdf("pdf_s",    "pdf_s",    ROOT.RooArgSet(x), h_s)
  pdf_ww   =  ROOT.RooHistPdf  ("pdf_ww",   "pdf_ww", ROOT.RooArgSet(x), h_ww)
  pdf_top  =  ROOT.RooHistPdf ("pdf_top",  "pdf_top", ROOT.RooArgSet(x), h_top)
  pdf_dytt =  ROOT.RooHistPdf("pdf_dytt", "pdf_dytt", ROOT.RooArgSet(x), h_dytt)
  pdf_vv   =  ROOT.RooHistPdf("pdf_vv",   "pdf_vv",   ROOT.RooArgSet(x), h_vv)

  # we now build a model.
  # data are represented by the sum of each of the above PDFs, 
  # each one multiplied by its integral 
  model = ROOT.RooAddPdf("model", "model", ROOT.RooArgSet(pdf_s,  pdf_ww,  pdf_top,  pdf_dytt,  pdf_vv),
                                           ROOT.RooArgSet(norm_s, norm_ww, norm_top, norm_dytt, norm_vv))

  # we now impose LogNormal constraints on the nuisances
  constraint_ww   = ROOT.RooLognormal("constraint_ww",   "constraint_ww",   mu_ww,   ROOT.RooFit.RooConst(1), ROOT.RooFit.RooConst(1.1))
  constraint_top  = ROOT.RooLognormal("constraint_top",  "constraint_top",  mu_top,  ROOT.RooFit.RooConst(1), ROOT.RooFit.RooConst(1.1))
  constraint_dytt = ROOT.RooLognormal("constraint_dytt", "constraint_dytt", mu_dytt, ROOT.RooFit.RooConst(1), ROOT.RooFit.RooConst(1.1))
  constraint_vv   = ROOT.RooLognormal("constraint_vv",   "constraint_vv",   mu_vv,   ROOT.RooFit.RooConst(1), ROOT.RooFit.RooConst(1.1))

  constraints = ROOT.RooArgSet(constraint_ww, constraint_top, constraint_dytt, constraint_vv)
  nuisances   = ROOT.RooArgSet(mu_ww, mu_top, mu_dytt, mu_vv)

  constrained_model = ROOT.RooProdPdf("constrained_model", "constrained_model", ROOT.RooArgSet(model,constraints))

  # save everything to a workspace for later use
  w = ROOT.RooWorkspace("w"+label, "w"+label)
  getattr(w,'import')(constrained_model)
  # also give a name to this set, so we can reues it later
  w.defineSet("constraints", constraints);
  w.defineSet("nuisances", nuisances);
  w.saveSnapshot("bonly", ROOT.RooArgSet(mu, nuisances))

  return w

