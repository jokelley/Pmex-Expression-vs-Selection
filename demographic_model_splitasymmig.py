import numpy
import dadi

# Demographic model with two populations (one population split into two populations at time T)
# Asymmetric migration rates between the populations (m12 and m21)
# Effective population sizes are represented by the parameters nu1 and nu2

def split_asym_mig(params, ns, pts):
    nu1,nu2,T,m12,m21 = params
    xx = yy = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D (xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs
