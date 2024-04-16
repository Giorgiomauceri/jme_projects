
from bamboo.plots import Plot, CutFlowReport
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op
from bamboo.analysisutils import forceDefine

import src.definitions as defs

from modules.baseModule import NanoBaseJME


class TTbarModule(NanoBaseJME):
    """"""

    def __init__(self, args):
        super(TTbarModule, self).__init__(args)

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []
        yields = CutFlowReport("yields", printInLog=True, recursive=True)
        plots.append(yields)
        yields.add(noSel, 'No Selection')


        # # AK8 Jets
        # ak8Jets = op.sort(
        #     op.select(tree.FatJet, lambda jet: defs.ak8jetDef(jet)), lambda jet: -jet.pt)

        # ak8JetsID = op.sort(
        #     ak8Jets, lambda jet: jet.jetId & 2)


        muons, electrons, clElectrons, ak4Jets, clak4Jets, ak4JetsID, ak4Jetspt40, ak4Jetspt100, ak4Jetsetas2p4, ak4Jetsetag2p4, ak8Jets, clak8Jets, clak4genjets = defs.defineObjects(tree)

        # 1 lepton channel selection
        hasOneSFLepton = noSel.refine('hasOneSFLepton', cut=
                op.AND(op.rng_len(muons) == 1, op.rng_len(clElectrons) ==0, 
                       #muons[0].pt > 25.
            ))

        ### Di-leptonic channel ###

        # 2 muon and 2 electron selection
        hasTwoSFLeptons = noSel.refine('hasTwoSFLeptons', cut=(
            op.OR(
                op.AND(op.rng_len(muons) == 2, op.rng_len(clElectrons) ==0, 
                       muons[0].charge != muons[1].charge, 
                       muons[0].pt > 25., muons[1].pt > 15.),
                op.AND(op.rng_len(clElectrons) == 2, op.rng_len(muons) ==0,
                       clElectrons[0].charge != clElectrons[1].charge, 
                       clElectrons[0].pt > 25., clElectrons[1].pt > 15. )
            )
        ))



        # Neutrino Reconstruction
        bjets = op.select(tree.Jet, lambda jet: jet.btagDeepFlavCvB>=0.8)
        #bjets = []
        #for x in range(op.rng_len()+1):
        #    if (tree.Jet[x].btagDeepFlavCvB>=0.8):
        #        bjets.append(tree.Jet[x])
        bjet =  op.rng_min_element_by(bjets, lambda bjet: op.deltaR(bjet.p4,muons[0].p4))

        
        Mass_W = 80.399
        
        #mu = op.sum(op.pow(Mass_W, 2)/2., op.product(muons[0].pt, tree.PuppiMET.pt)) 
        #A = op.product(muons[0].p4.Pt(), muons[0].p4.Pt())
        #B = op.product(mu, muons[0].p4.Pz())
        #C = op.sum(op.product(mu, mu), op.product(-op.pow(muons[0].p4.E(), 2), op.pow(tree.PuppiMET.pt, 2)))
        #Discr = op.sum(op.pow(B, 2),-op.product(A, C))
        
        #neutrinosol1=op.multiSwitch((Discr<=0,-B/A), (Discr>0, -op.sum(B, op.sqrt(Discr))/A))
        #neutrinosol2=op.multiSwitch((Discr<=0,-B/A), (Discr>0, op.sum(-B, op.sqrt(Discr))/A))
                
        A = Mass_W*Mass_W/2. + tree.PuppiMET.p4.Px()*muons[0].p4.Px() + tree.PuppiMET.p4.Py()*muons[0].p4.Py()      
        Discr = muons[0].p4.E() * muons[0].p4.E() * (A * A - tree.PuppiMET.pt * tree.PuppiMET.pt * muons[0].pt*muons[0].pt)
        neutrinosol1=op.multiSwitch((Discr<=0,op.product(A, muons[0].p4.Pz())/op.pow(muons[0].pt, 2)), 
                                    (Discr>0, op.sum(-op.product(A, muons[0].p4.Pz()), op.sqrt(Discr))/op.pow(muons[0].pt, 2)))
        neutrinosol2=op.multiSwitch((Discr<=0,op.product(A, muons[0].p4.Pz())/op.pow(muons[0].pt, 2)), 
                                    (Discr>0, op.sum(-op.product(A, muons[0].p4.Pz()), -op.sqrt(Discr))/op.pow(muons[0].pt, 2)))
        
        neutrino1 = op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", 
                                (tree.PuppiMET.p4.Px(), tree.PuppiMET.p4.Py(), neutrinosol1, op.sqrt(op.pow(tree.PuppiMET.p4.E(), 2)+op.pow(neutrinosol1, 2))))           
        neutrino2 = op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", 
                                (tree.PuppiMET.p4.Px(), tree.PuppiMET.p4.Py(), neutrinosol2, op.sqrt(op.pow(tree.PuppiMET.p4.E(), 2)+op.pow(neutrinosol2, 2))))            
                      
        Wboson1 = op.sum(neutrino1, muons[0].p4)
        Wboson2 = op.sum(neutrino2, muons[0].p4)
        Topcand1 = op.sum(Wboson1, bjet.p4)
        Topcand2 = op.sum(Wboson2, bjet.p4)
        Topq_M = op.switch(op.abs(op.sum(Topcand1.M(), -172.76))<op.abs(op.sum(Topcand2.M(), -172.76)), Topcand1.E(), Topcand2.E())

        Top_mass_hist = Plot.make1D(f"noSel_Top_MASS", Topq_M, noSel, EqBin(40, 0., 500.), xTitle=f"Top Mass")
        Top_mass_hist_twolep = Plot.make1D(f"hasTwoSFLeptons_Top_MASS", Topq_M, hasTwoSFLeptons, EqBin(40, 0., 500.), xTitle=f"Top Mass")
        Top_mass_hist_onelep = Plot.make1D(f"hasOneSFLepton_Top_MASS", Topq_M, hasOneSFLepton, EqBin(40, 0., 500.), xTitle=f"Top Mass")

###############################################

        plots = []
        plots+=[Top_mass_hist, Top_mass_hist_twolep, Top_mass_hist_onelep]
        return plots


