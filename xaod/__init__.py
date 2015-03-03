# Set up ROOT and RootCore:
import os
import sys

BASE_DIR = os.getenv('DIR_HIGGSTAUTAU_SETUP')
if not BASE_DIR:
    sys.exit('You did not source setup.sh!')

CACHE_DIR = os.path.join(BASE_DIR, 'cache')

import ROOT
ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
# Initialize the xAOD infrastructure: 
RC = ROOT.xAOD.Init()
if not RC.isSuccess():
    raise RuntimeError('Error in the xAOD initialization')

ROOT.xAOD.AuxContainerBase()

# Set up the input files:
# ftemp = ROOT.TFile(os.path.join(CACHE_DIR, 'xaod_struct.root'))
# ttemp = ROOT.xAOD.MakeTransientTree(ftemp)
# ftemp.Close()

# HACK HACK HACK
# DataVector __iter__ does not wrok in ROOT 6 :-(
def __iter__(self):
    for idx in xrange(self.size()):
        yield self.at(idx)
ROOT.DataVector(ROOT.xAOD.Electron_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.Jet_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.Vertex_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.TauJet_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.MissingET_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.Muon_v1).__iter__ = __iter__
ROOT.DataVector(ROOT.xAOD.TruthParticle_v1).__iter__ = __iter__


TOOLS = ROOT.asg.ToolStore()
