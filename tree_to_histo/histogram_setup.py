
def get_eta_intervals(max_eta = -2.5, min_eta = 2.5, n_eta = 100):
    """
    get a vector of eta intervals for efficiency and fake-rate plots
    """
    
    etaintervals = []
    eta_step = (max_eta - min_eta)/n_eta
    
    for k in range (0,n_eta):
        d = min_eta + k*eta_step;
        etaintervals.append(d)
    etaintervals.append(max_eta)

    #print etaintervals
    return etaintervals

