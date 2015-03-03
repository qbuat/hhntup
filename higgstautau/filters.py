import ROOT

from rootpy.tree.filtering import *
from rootpy.extern.hep import pdg
from rootpy import stl
VectorTLorentzVector = stl.vector("TLorentzVector")
Vector = stl.vector('float')
from itertools import ifilter
from math import *
from array import array as carray

from xaod import TOOLS
from . import datasets
# from .corrections import reweight_ggf
from .units import GeV
from .tautools import TauDecay
from . import utils
from . import store_helper
from . import log; log = log[__name__]

from goodruns import GRL


BCH_TOOLS = []


class GRLFilter(EventFilter):

    def __init__(self, grl, **kwargs):
        super(GRLFilter, self).__init__(**kwargs)
        if isinstance(grl, GRL):
            self.grl = grl
        else:
            self.grl = GRL(grl)
    def passes(self, event):
        if not self.grl:
            return False
        return (int(event.EventInfo.runNumber()), int(event.EventInfo.lumiBlock())) in self.grl

class GRLFilterOfficial(EventFilter):
    # ALTERNATIVE USING OFFICIAL TOOL
    def __init__(self, xml, **kwargs):
        super(GRLFilterOfficial, self).__init__(**kwargs)
        from ROOT import Root
        reader = Root.TGoodRunsListReader(xml)
        reader.Interpret()
        self.grl = reader.GetMergedGRLCollection()
        self.grl.Summary()
    def passes(self, event):
        return self.grl.HasRunLumiBlock(
            int(event.EventInfo.runNumber()), 
            int(event.EventInfo.lumiBlock()))


def primary_vertex_selection(vxp):
    return vxp.vertexType() == 1 and vxp.nTrackParticles() >= 4


def pileup_vertex_selection(vxp):
    return vxp.vertexType() == 3 and vxp.nTrackParticles() >= 2


def vertex_selection(vxp):
    """ Does the full primary and pileup vertex selection """
    return primary_vertex_selection(vxp) or pileup_vertex_selection(vxp)


class PriVertex(EventFilter):

    def passes(self, event):
        event.vertices.select(vertex_selection)
        return any(ifilter(primary_vertex_selection, event.vertices))


class CoreFlags(EventFilter):

    def passes(self, event):
        Core = event.EventInfo.Core
        return event.EventInfo.errorState(Core) == 0


class NvtxJets(EventFilter):

    def __init__(self, tree, **kwargs):
        super(NvtxJets, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        # Check for a good primary vertex
        # This is needed for jet and soft term systematics
        goodPV = False
        nvtxsoftmet = 0
        nvtxjets = 0
        for vertex in event.vertices:
            if vertex.vertexType() == 1 and vertex.nTrackParticles() > 2 and abs(vertex.z()) < 200:
                goodPV = True
        if goodPV:
            for vertex in event.vertices:
                if vertex.nTrackParticles() > 2:
                    nvtxsoftmet += 1
                if vertex.nTrackParticles() > 1:
                    nvtxjets += 1
        self.tree.nvtxsoftmet = nvtxsoftmet
        self.tree.nvtxjets = nvtxjets
        return True


class BCHCleaning(EventFilter):
    # NOT CONVERTED TO XAOD YET
    """
    https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BCHCleaningTool
    """
    def __init__(self, tree, passthrough, datatype, **kwargs):
        if not passthrough:
            # from externaltools import TileTripReader
            # from externaltools import BCHCleaningTool
            from ROOT import Root
            from ROOT import BCHTool
            self.tree = tree
            self.datatype = datatype
            self.tiletool = Root.TTileTripReader()
            self.tiletool.setTripFile(TileTripReader.get_resource("CompleteTripList_2011-2012.root"))
            self.bchtool_data = BCHTool.BCHCleaningToolRoot()
            self.bchtool_mc = BCHTool.BCHCleaningToolRoot()
            self.bchtool_data.InitializeTool(
                True, self.tiletool, BCHCleaningTool.get_resource("FractionsRejectedJetsMC.root"))
            self.bchtool_mc.InitializeTool(
                False, self.tiletool, BCHCleaningTool.get_resource("FractionsRejectedJetsMC.root"))
            BCH_TOOLS.append(self.bchtool_data)
            BCH_TOOLS.append(self.bchtool_mc)
        super(BCHCleaning, self).__init__(passthrough=passthrough, **kwargs)

    def passes(self, event):
        if self.datatype in (datasets.DATA, datasets.MC, datasets.MCEMBED):
            if self.datatype == datasets.DATA:
                jet_tool = self.bchtool_data
                runnumber = event.EventInfo.runNumber()
                lbn = event.EventInfo.lumiBlock()
            else:
                jet_tool = self.bchtool_mc
                runnumber = self.tree.RunNumber
                lbn = self.tree.lbn
                #jet_tool.SetSeed(314159 + event.EventNumber * 2)
            for jet in event.jets:
                jet.BCHMedium = jet_tool.IsBadMediumBCH(
                    runnumber, lbn, jet.eta(), jet.phi(), 
                    jet.auxdataConst('float')('BchCorrCell'), 
                    jet.auxdataConst('float')('EMFrac'), jet.pt())
                jet.BCHTight = jet_tool.IsBadTightBCH(
                    runnumber, lbn, jet.eta(), jet.phi(), 
                    jet.auxdataConst('float')('BchCorrCell'), 
                    jet.auxdataConst('float')('EMFrac'), jet.pt())
            for tau in event.taus:
                tau.BCHMedium = jet_tool.IsBadMediumBCH(
                    runnumber, lbn, tau.jet().eta(), tau.jet().phi(), 
                    tau.jet().auxdataConst('float')('BchCorrCell'),
                    tau.jet().auxdataConst('float')('EMFrac'), tau.jet().pt())
                tau.BCHTight = jet_tool.IsBadTightBCH(
                    runnumber, lbn, tau.jet().eta(), tau.jet().phi(), 
                    tau.jet().auxdataConst('float')('BchCorrCell'),
                    tau.jet().auxdataConst('float')('EMFrac'), tau.jet().pt())

        elif self.datatype == datasets.EMBED:
            # Do truth-matching to find out if MC taus
            #self.bchtool_data.SetSeed(314159 + event.EventNumber * 2)
            #self.bchtool_mc.SetSeed(314159 + event.EventNumber * 3)
            runnumber = event.RunNumber
            lbn = event.lbn
            for jet in event.jets:
                if jet.matched:
                    jet_tool = self.bchtool_mc
                else:
                    jet_tool = self.bchtool_data
                jet.BCHMedium = jet_tool.IsBadMediumBCH(
                    runnumber, lbn, jet.eta(), jet.phi(), 
                    jet.auxdataConst('float')('BchCorrCell'), 
                    jet.auxdataConst('float')('EMFrac'), jet.pt())
                jet.BCHTight = jet_tool.IsBadTightBCH(
                    runnumber, lbn, jet.eta(), jet.phi(), 
                    jet.auxdataConst('float')('BchCorrCell'), 
                    jet.auxdataConst('float')('EMFrac'), jet.pt())
            for tau in event.taus:
                if tau.matched:
                    jet_tool = self.bchtool_mc
                else:
                    jet_tool = self.bchtool_data
                tau.BCHMedium = jet_tool.IsBadMediumBCH(
                    runnumber, lbn, tau.jet().eta(), tau.jet().phi(), 
                    tau.jet().auxdataConst('float')('BchCorrCell'),
                    tau.jet().auxdataConst('float')('EMFrac'), tau.jet().pt())
                tau.BCHTight = jet_tool.IsBadTightBCH(
                    runnumber, lbn, tau.jet().eta(), tau.jet().phi(), 
                    tau.jet().auxdataConst('float')('BchCorrCell'),
                    tau.jet().auxdataConst('float')('EMFrac'), tau.jet().pt())

        return True


class TileTrips(EventFilter):
    """
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/DataPreparationCheckListForPhysicsAnalysis#Rejection_of_bad_corrupted_event
    """
    def __init__(self, passthrough=False, **kwargs):
        if not passthrough:
            #from externaltools import TileTripReader
            from ROOT import Root
            self.tool = Root.TTileTripReader()
        super(TileTrips, self).__init__(passthrough=passthrough, **kwargs)

    def passes(self, event):
        return self.tool.checkEvent(
            event.EventInfo.runNumber(),
            event.EventInfo.lumiBlock(),
            event.EventInfo.eventNumber())


class JetCleaning(EventFilter):

    BAD_TILE = [
        202660, 202668, 202712, 202740, 202965, 202987, 202991, 203027, 203169
    ]

    def __init__(self,
                 datatype,
                 year,
                 pt_thresh=20*GeV,
                 eta_max=4.5,
                 **kwargs):
        super(JetCleaning, self).__init__(**kwargs)
        self.year = year
        self.datatype = datatype
        self.pt_thresh = pt_thresh
        self.eta_max = eta_max
        from ROOT import JetCleaningTool
        self.tool = JetCleaningTool(JetCleaningTool.LooseBad)

    def passes(self, event):
        # using LC jets
        for jet in event.jets:
            if jet.pt() <= self.pt_thresh or abs(jet.eta()) >= self.eta_max:
                continue
            if self.tool.accept(jet):
                return False

        if (self.datatype in (datasets.DATA, datasets.EMBED)) and self.year == 2012:
            # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HowToCleanJets2012
            # Hot Tile calorimeter in period B1 and B2
            if event.EventInfo.runNumber() in JetCleaning.BAD_TILE:
                # recommendation is to use EM jets
                for jet in event.jets_EM:
                    _etaphi28 = (
                        -0.2 < jet.eta() < -0.1 and
                        2.65 < jet.phi() < 2.75)
                    FracSamplingMax = jet.auxdataConst('float')('FracSamplingMax')
                    SamplingMax = jet.auxdataConst('int')('FracSamplingMaxIndex')
                    if FracSamplingMax > 0.6 and SamplingMax == 13 and _etaphi28:
                        return False

        return True


class LArError(EventFilter):

    def passes(self, event):
        LAr = event.EventInfo.LAr
        return event.EventInfo.errorState(LAr) == 0


class TileError(EventFilter):

    def passes(self, event):
        Tile = event.EventInfo.Tile
        return event.EventInfo.errorState(Tile) == 0


class JetCrackVeto(EventFilter):

    def passes(self, event):
        for jet in event.jets:
            if jet.pt() <= 20 * GeV: 
                continue
            if 1.3 < abs(jet.eta()) < 1.7: 
                return False
        return True


class TauElectronVeto(EventFilter):

    def __init__(self, min_taus, **kwargs):
        super(TauElectronVeto, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):
        #https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauRecommendationsWinterConf2013#Electron_veto
        # only apply eveto on 1p taus with cluster and track eta less than 2.47
        # Eta selection already applied by TauEta filter
        event.taus.select(lambda tau:
            tau.nTracks() > 1 or
            tau.isTau(ROOT.xAOD.TauJetParameters.EleBDTLoose) == 0)
        return len(event.taus) >= self.min_taus


class TauMuonVeto(EventFilter):

    def __init__(self, min_taus, **kwargs):
        super(TauMuonVeto, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):
        event.taus.select(lambda tau: tau.isTau(ROOT.xAOD.TauJetParameters.MuonVeto) == 0)
        return len(event.taus) >= self.min_taus


class TauHasTrack(EventFilter):

    def __init__(self, min_taus, **kwargs):
        super(TauHasTrack, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):
        event.taus.select(lambda tau: tau.nTracks() > 0)
        return len(event.taus) >= self.min_taus


class TauPT(EventFilter):

    def __init__(self, min_taus, thresh=20 * GeV, **kwargs):
        self.min_taus = min_taus
        self.thresh = thresh
        super(TauPT, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.pt() > self.thresh)
        return len(event.taus) >= self.min_taus


class TauEta(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauEta, self).__init__(**kwargs)

    def passes(self, event):
        # both calo and leading track eta within 2.47
        event.taus.select(lambda tau:
            abs(tau.eta()) < 2.47 and
            abs(tau.track(0).eta()) < 2.47)
        return len(event.taus) >= self.min_taus

def jvf_selection(tau):
    if abs(tau.track(0).eta()) > 2.1:
        return True
    else:
        jvf = 0
        jvf_vec = tau.jet().auxdataConst('std::vector<float, std::allocator<float> >')('JVF')
        if not jvf_vec.empty():
            jvf = jvf_vec[0]
        return jvf > .5

class TauJVF(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        self.filter_func = jvf_selection
        super(TauJVF, self).__init__(**kwargs)
        
    def passes(self, event):
        event.taus.select(self.filter_func)
        return len(event.taus) >= self.min_taus


class Tau1Track3Track(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(Tau1Track3Track, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.nTracks() in (1, 3))
        return len(event.taus) >= self.min_taus


class Tau1P3P(EventFilter):
    """
    Only keep 1P + 3P and 3P + 3P
    """
    def passes(self, event):
        assert len(event.taus) == 2
        tau1, tau2 = event.taus
        # 1P + 3P
        if (tau1.nTracks() == 1 and tau2.nTracks() == 3) or \
           (tau1.nTracks() == 3 and tau2.nTracks() == 1):
            return True
        # 3P + 3P
        if tau1.nTracks() == 3 and tau2.nTracks() == 3:
            return True
        return False


class TauCharge(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauCharge, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: abs(tau.charge()) == 1)
        return len(event.taus) >= self.min_taus


class TauIDLoose(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauIDLoose, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.isTau(ROOT.xAOD.TauJetParameters.JetBDTSigLoose) == 1)
        return len(event.taus) >= self.min_taus


class TauIDMedium(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauIDMedium, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.isTau(ROOT.xAOD.TauJetParameters.JetBDTSigMedium) == 1)
        return len(event.taus) >= self.min_taus


class TauCrack(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauCrack, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(
            lambda tau: not (
                1.37 <= abs(tau.track(0).eta()) <= 1.52))
        return len(event.taus) >= self.min_taus

class TrueTauSelection(EventFilter):
    """
    True tau selection from the truth particle container
    using the official tool (does it work for all the generators ?)
    """
    def __init__(self, passthrough=False, **kwargs):
        super(TrueTauSelection, self).__init__(
            passthrough=passthrough, **kwargs)
        if not passthrough:
            from ROOT.TauAnalysisTools import TauTruthMatchingTool
            self.tau_truth_tool = TauTruthMatchingTool('tau_truth_tool')
            # Should add an argument for the sample type
            self.tau_truth_tool.initialize()
            truth_matching_tool = self.tau_truth_tool

    def passes(self, event):
        self.tau_truth_tool.setTruthParticleContainer(event.truetaus.collection)
        self.tau_truth_tool.createTruthTauContainer()
        truth_taus = self.tau_truth_tool.getTruthTauContainer()
        truth_taus_aux = self.tau_truth_tool.getTruthTauAuxContainer()
        truth_taus.setNonConstStore(truth_taus_aux)
        event.truetaus.collection = truth_taus
        # OLD METHOD using the edm itself
        # event.truetaus.select(lambda p: p.isTau() and p.status() in (2,))
        return True

class TruthMatching(EventFilter):

    def __init__(self, passthrough=False, **kwargs):
        super(TruthMatching, self).__init__(
            passthrough=passthrough, **kwargs)
        if not passthrough:
            # This implies that this filter is ALWAYS applied after
            # the TrueTauSelection
            self.tau_truth_tool = TOOLS.get('tau_truth_tool')

    def passes(self, event):
        for tau in event.taus:
            tau.matched_dr = 1111.
            tau.matched_object = None
            self.tau_truth_tool.setTruthParticleContainer(event.mc.collection)
            true_tau = self.tau_truth_tool.applyTruthMatch(tau)            
            tau.matched = tau.auxdataConst('bool')('IsTruthMatched')
            if tau.matched:
                tau.matched_object = true_tau
                tau.matched_dr = utils.dR(
                    tau.eta(), tau.phi(),
                    true_tau.auxdataConst('double')('eta_vis'),
                    true_tau.auxdataConst('double')('phi_vis'))
        return True


class RecoJetTrueTauMatching(EventFilter):
    """
    To use after the TrueTauSelection filter
    """

    def passes(self, event):
        for jet in event.jets:
            jet.matched = False
            jet.matched_dr = 1111.
            jet.matched_object = None
            for p in event.truetaus:
                dr = utils.dR(
                    jet.eta(), jet.phi(), 
                    p.auxdataConst('double')('eta_vis'),
                    p.auxdataConst('double')('phi_vis'))
                if dr < 0.2:
                    # TODO: handle possible collision!
                    jet.matched = True
                    jet.matched_dr = dr
                    jet.matched_object = p
                    break
        return True


class TauCalibration(EventFilter):
    """
    Apply Energy shift in data and 
    systematic variation in MC (Not yet)
    """
    def __init__(self, datatype, **kwargs):
        super(TauCalibration, self).__init__(**kwargs)
        self.datatype = datatype

        from ROOT.TauAnalysisTools import TauSmearingTool
        self.tool = TauSmearingTool('tau_smearing_tool')
        self.tool.setProperty('bool')(
            'IsData', self.datatype == datasets.DATA)

    def passes(self, event):
        taus_copy = store_helper.shallowCopyTauJetContainer(event.taus.collection)
        for tau in taus_copy:
            self.tool.applyCorrection(tau)
        event.taus.collection = taus_copy
        return True


class TauJetOverlapRemoval(EventFilter):
    """
    Precendence: taus > jets
    """
    def __init__(self, dr=.2, **kwargs):
        super(TauJetOverlapRemoval, self).__init__(**kwargs)
        self.dr = dr

    def passes(self, event):
        # remove overlap with taus
        event.jets.select(lambda jet:
                not any([tau for tau in event.taus if
                (utils.dR(jet.eta(), jet.phi(), tau.eta(), tau.phi()) < self.dr)]))
        return True

class NumJets25(EventFilter):

    def __init__(self, tree, **kwargs):
        super(NumJets25, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        self.tree.numJets25 = len([j for j in event.jets if
            j.pt() > 25 * GeV and abs(j.eta()) < 4.5])
        return True


class NonIsolatedJet(EventFilter):
    """
    https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=200403
    """
    def __init__(self, tree, **kwargs):
        super(NonIsolatedJet, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        # only write flag instead of vetoing the event so this
        # can be turned on and off after
        self.tree.nonisolatedjet = False
        for tau in event.taus:
            for jet in event.jets:
                if 0.4 < utils.dR(tau.eta(), tau.phi(), jet.eta(), jet.phi()) < 1.0:
                    self.tree.nonisolatedjet = True
        return True


def jet_selection_2011(jet):
    """ Finalizes the jet selection """

    if not (jet.pt() > 25 * GeV):
        return False

    if not (abs(jet.eta()) < 4.5):
        return False

    # suppress forward jets
    if (abs(jet.eta()) > 2.4) and not (jet.pt() > 30 * GeV):
        return False

    # JVF cut on central jets
    #if (abs(jet.eta) < 2.4) and not (abs(jet.jvtxf) > 0.75):
    #    return False
    # NO JVFUncertaintyTool for 2011!

    return True


def jet_selection_2012(jet):
    """ Finalizes the jet selection
    https://cds.cern.ch/record/1472547/files/ATL-COM-PHYS-2012-1202.pdf
    """
    if not (jet.pt() > 30 * GeV):
        return False

    if not (abs(jet.eta()) < 4.5):
        return False

    # suppress forward jets
    if (abs(jet.eta()) > 2.4) and not (jet.pt() > 35 * GeV):
        return False

    # JVF cut on central jets below 50 GeV
    if (jet.pt() < 50 * GeV) and (abs(jet.eta()) < 2.4):
        jvf_vec = jet.auxdataConst('std::vector<float, std::allocator<float> >')('JVF')
        jvf = 0 if jvf_vec.empty() else jvf_vec[0]
        if not (abs(jvf) > 0.5):
            return False

    return True


class JetSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def __init__(self, year, **kwargs):
        if year == 2011:
            self.filter_func = jet_selection_2011
        elif year == 2012:
            self.filter_func = jet_selection_2012
        else:
            raise ValueError("No jet selection defined for year %d" % year)
        super(JetSelection, self).__init__(**kwargs)

    def passes(self, event):
        event.jets.select(self.filter_func)
        return True


class JetPreselection(EventFilter):

    def passes(self, event):
        event.jets.select(lambda jet: jet.pt() > 20 * GeV)
        return True


class MCWeight(EventFilter):
    # NOT FULLY CONVERTED TO XAOD YET

    def __init__(self, datatype, tree, **kwargs):
        self.datatype = datatype
        self.tree = tree
        super(MCWeight, self).__init__(**kwargs)

    def passes(self, event):
        # set the event weights
        if self.datatype == datasets.MC:
            truth_event = event.TruthEvent[0]
            self.tree.mc_weight = event.EventInfo.mcEventWeight()
            val_i = ROOT.Long(0)
            if truth_event.pdfInfoParameter(val_i, truth_event.PDFID1):
                self.tree.mcevent_pdf_id1_0 = val_i
            if truth_event.pdfInfoParameter(val_i, truth_event.PDFID2):
                self.tree.mcevent_pdf_id2_0 = val_i
            val_f = carray('f', [0.])
            if truth_event.pdfInfoParameter(val_f, truth_event.X1):
                self.tree.mcevent_pdf_x1_0 = val_f[0]
            if truth_event.pdfInfoParameter(val_f, truth_event.X2):
                self.tree.mcevent_pdf_x2_0 = val_f[0] 
            if truth_event.pdfInfoParameter(val_f, truth_event.scalePDF):
                self.tree.mcevent_pdf_scale_0 = val_f[0]

        elif self.datatype == datasets.EMBED:
            self.tree.mc_weight = event.EventInfo.mcEventWeight()
        return True


class ggFReweighting(EventFilter):
    # NOT CONVERTED TO XAOD YET

    def __init__(self, dsname, tree, **kwargs):
        self.dsname = dsname
        self.tree = tree
        super(ggFReweighting, self).__init__(**kwargs)

    def passes(self, event):
        # self.tree.ggf_weight = reweight_ggf(event, self.dsname)
        return True


class JetIsPileup(EventFilter):
    # NOT converted to XAOD yet
    """
    must be applied before any jet selection
    """
    def __init__(self, **kwargs):
        super(JetIsPileup, self).__init__(**kwargs)
        if not self.passthrough:
            # from externaltools import JVFUncertaintyTool as JVFUncertaintyTool2012
            from ROOT import JVFUncertaintyTool
            self.tool = JVFUncertaintyTool("AntiKt4LCTopo")

    def passes(self, event):
        # collect truth jets
        truejets = VectorTLorentzVector()
        for truejet in event.truejets:
            if truejet.pt() > 10e3:
                truejets.push_back(truejet.p4())
        # test each jet
        for jet in event.jets:
            ispileup = self.tool.isPileUpJet(jet.p4(), truejets)
            jet.ispileup = ispileup
        return True


class JetCopy(EventFilter):
    # NOT CONVERTED TO XAOD YET
    # NOT NEEDED ANYMORE

    def __init__(self, tree, **kwargs):
        super(JetCopy, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        tree = self.tree
        tree.jet_E_original.clear()
        tree.jet_m_original.clear()
        tree.jet_pt_original.clear()
        tree.jet_eta_original.clear()
        tree.jet_phi_original.clear()
        for jet in event.jets:
            tree.jet_E_original.push_back(jet.E)
            tree.jet_m_original.push_back(jet.m)
            tree.jet_pt_original.push_back(jet.pt)
            tree.jet_eta_original.push_back(jet.eta)
            tree.jet_phi_original.push_back(jet.phi)
        return True


class HiggsPT(EventFilter):
    # NOT CONVERTED TO XAOD YET

    def __init__(self, year, tree, **kwargs):
        super(HiggsPT, self).__init__(**kwargs)
        self.tree = tree
        if year == 2011:
            self.status = (2, 10902, 62)
        elif year == 2012:
            self.status = (62, 195)
        else:
            raise ValueError("No HiggsPT defined for year {0}".format(year))
        
    def passes(self, event):
        pt = 0
        higgs = None
        status = self.status
        # find the Higgs
        for mc in event.mc:
            if mc.pdgId() == 25 and mc.status() in status:
                pt = mc.pt()
                higgs = mc
                break
        if higgs is None:
            raise RuntimeError("Higgs not found!")
        self.tree.true_resonance_pt = pt
        # Only consider taus here since there are very soft photons radiated
        # off the taus but included as children of the Higgs
        vertex = higgs.decayVtx()
        children = [
            vertex.outgoingParticle(i) for i in
            xrange(vertex.nOutgoingParticles())]
        true_taus = [TauDecay(mc).fourvect_visible
                     for mc in children
                     if mc.pdgId() in (pdg.tau_plus, pdg.tau_minus)
                     and mc.status() in (2, 11, 195)]
        # The number of anti kt R = 0.4 truth jets with pT>25 GeV, not
        # originating from the decay products of the Higgs boson.
        # Start from the AntiKt4Truth collection. Reject any jet with pT<25
        # GeV. Reject any jet withing dR < 0.4 of any electron, tau, photon or
        # parton (directly) produced in the Higgs decay.
        jets = [jet for jet in event.truejets if jet.pt() >= 25 * GeV
                and not any([tau for tau in true_taus if
                             utils.dR(jet.eta(), jet.phi(),
                                      tau.Eta(), tau.Phi()) < 0.4])]
        # Count the number of remaining jets
        self.tree.num_true_jets_no_overlap = len(jets)
        if len(jets) >=2:
            jet1, jet2 = jets[:2]
            self.tree.true_jet1_no_overlap_pt = jet1.pt()
            self.tree.true_jet2_no_overlap_pt = jet2.pt()
            self.tree.true_dEta_jet1_jet2_no_overlap = abs(jet1.eta() - jet2.eta())
            self.tree.true_mass_jet1_jet2_no_overlap = (jet1.p4() + jet2.p4()).M()
            self.tree.true_dphi_jj_higgs_no_overlap = abs(utils.dphi(higgs.phi(), (jet1.p4() + jet2.p4()).Phi()))
        return True


class BCHSampleRunNumber(EventFilter):
    # NOT CONVERTED TO XAOD YET

    """
    d3pd.RunNumber=195848 tells the tool that the sample was made with mc12b
    pileup conditions. Our BCH samples were made with identical pileup
    conditions, but for reasons unclear, they were assigned
    d3pd.RunNumber=212399, and the pileup tool does not know what to do with
    this RunNumber.
    """
    def passes(self, event):
        event.RunNumber = 195848
        return True


class ClassifyInclusiveHiggsSample(EventFilter):
    # NOT CONVERTED TO XAOD YET

    UNKNOWN, TAUTAU, WW, ZZ, BB = range(5)

    def __init__(self, tree, **kwargs):
        super(ClassifyInclusiveHiggsSample, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        higgs = None
        # find the Higgs
        for mc in event.mc:
            if mc.pdgId() == 25 and mc.status() == 62:
                pt = mc.pt()
                higgs = mc
                break
        if higgs is None:
            raise RuntimeError("Higgs not found!")
        decay_type = self.UNKNOWN
        # check pdg id of children
        # for mc in higgs.iter_children():
        #     if mc.pdgId in (pdg.tau_minus, pdg.tau_plus):
        #         decay_type = self.TAUTAU
        #         break
        #     elif mc.pdgId in (pdg.W_minus, pdg.W_plus):
        #         decay_type = self.WW
        #         break
        #     elif mc.pdgId == pdg.Z0:
        #         decay_type = self.ZZ
        #         break
        #     elif mc.pdgId in (pdg.b, pdg.anti_b):
        #         decay_type = self.BB
        #         break
        # self.tree.higgs_decay_channel = decay_type
        return True
