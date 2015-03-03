import math
import array
import ROOT
from rootpy.vector import LorentzVector


def is_MET_bisecting(dphi_taus, dphi_tau1_MET, dphi_tau2_MET):
    """
    check whether the MET is bisecting the
    two taus in the transverse plane

    dphi_taus = between the 2 taus
    dphi_tau1_MET = between tau 1 and MET
    dphi_tau2_MET = between tau 2 and MET
    """
    return (((max(dphi_tau1_MET, dphi_tau2_MET) <= dphi_taus) and
            (dphi_tau1_MET + dphi_tau2_MET <= math.pi)))


def mass(tau1, tau2, METpx, METpy):
    """
    Calculate and return the collinear mass and momentum fractions
    of tau1 and tau2

    TODO: set visible mass of taus. 1.2 GeV for 3p and 0.8 GeV for 1p
    """
    recTau1 = LorentzVector()
    recTau2 = LorentzVector()
    print tau1, tau2
    # tau 4-vector; synchronize for MMC calculation
    if tau1.nTracks() < 3:
        recTau1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), 800.) # MeV
    else:
        recTau1.SetPtEtaPhiM(tau1.pt(), tau1.eta(), tau1.phi(), 1200.) # MeV

    if tau2.nTracks() < 3:
        recTau2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), 800.) # MeV
    else:
        recTau2.SetPtEtaPhiM(tau2.pt(), tau2.eta(), tau2.phi(), 1200.) # MeV

    K = ROOT.TMatrixD(2, 2)
    K.InsertRow(0, 0, array.array('d', [recTau1.Px(), recTau2.Px()]));
    K.InsertRow(1, 0, array.array('d', [recTau1.Py(), recTau2.Py()]));

    # deprecated in ROOT6
    # K[0][0] = recTau1.Px()
    # K[0][1] = recTau2.Px()
    # K[1][0] = recTau1.Py(); 
    # K[1][1] = recTau2.Py()
    print K.Determinant()
    if K.Determinant() == 0:
        return -1., -1111., -1111.

    M = ROOT.TMatrixD(2, 1)
    M.InsertRow(0, 0, array.array('d', [METpx]))
    M.InsertRow(1, 0, array.array('d', [METpy]))
    print
    print METpx, METpy
    print
    print M[0][0], M[1][0]
    # M[0][0] = METpx
    # M[1][0] = METpy
    print 
    print recTau1.Px(), recTau2.Px()
    print recTau1.Py(), recTau2.Py()
    print
    print K[0][0], K[0][1]
    print K[1][0], K[1][1]
    print 'Something wrong after inversion'
    Kinv = K.Invert()
    print Kinv[0][0], Kinv[0][1]
    print Kinv[1][0], Kinv[1][1]
    print
    print Kinv.Determinant()
    print (K * Kinv)[0][0], (K * Kinv)[0][1]
    print (K * Kinv)[1][0], (K * Kinv)[1][1]

    X = Kinv * M

    X1 = X(0, 0)
    X2 = X(1, 0)

    x1 = 1./(1. + X1)
    x2 = 1./(1. + X2)

    p1 = recTau1 * (1. / x1)
    p2 = recTau2 * (1. / x2)
    m_col = (p1 + p2).M()
    m_vis = (recTau1 + recTau2).M()

    return m_vis, m_col, x1, x2
