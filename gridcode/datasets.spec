[__many__]
    container = string
    grl = string(default='')
    tree = string(default=tauPerf)
    type = option(DATA, MC, default=MC)
    label = option(TAU, ELEC, MUON, JET, default=TAU)
    weight = float(default=1.)
