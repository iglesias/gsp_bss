def crunch_hists(success_hist, hist):
  common_keys = set(success_hist._dict.keys()).intersection(hist._dict.keys())
  success_ratio = [1.0*success_hist._dict[k] / hist._dict[k] for k in common_keys]
  return success_ratio
