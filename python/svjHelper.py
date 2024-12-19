import os, math, sys, shutil
from string import Template
from glob import glob
import numpy as np
#import pandas as pd
#import ROOT as rt
#from scipy.interpolate import interp1d
#from scipy.special import exp1

class quark(object):
    def __init__(self,id,mass,charge):
        self.id = id
        self.mass = mass
        self.charge = charge
        self.massrun = mass
        self.color = 3
        self.bf = 1           #theory branching ratios
        self.bf_scaled = 1    #scaled branching ratios by rinv value
        self.on = True        #phase space allowed decay
        self.active = True    #for running nf

    def __repr__(self):
        return str(self.id)+": m = "+str(self.mass)+", mr = "+str(self.massrun)+", on = "+str(self.on)+", bf = "+str(self.bf)



# follows Ellis, Stirling, Webber calculations
class massRunner(object):
    def __init__(self):
        # QCD scale in GeV
        self.Lambda = 0.218

    # RG terms, assuming nc = 3 (QCD)
    def c(self): return 1./math.pi
    def cp(self,nf): return (303.-10.*nf)/(72.*math.pi)
    def b(self,nf): return (33.-2.*nf)/(12.*math.pi)
    def bp(self,nf): return (153.-19.*nf)/(2.*math.pi*(33.-2.*nf))
    def alphaS(self,Q,nf): return 1./(self.b(nf)*math.log(Q**2/self.Lambda**2))

    # derived terms
    def cb(self,nf): return 12./(33.-2.*nf)
    def one_c_cp_bp_b(self,nf): return 1.+self.cb(nf)*(self.cp(nf)-self.bp(nf))

    # constant of normalization
    def mhat(self,mq,nfq):
        return mq/math.pow(self.alphaS(mq,nfq),self.cb(nfq))/self.one_c_cp_bp_b(nfq)

    # mass formula
    def m(self,mq,nfq,Q,nf):
        # temporary hack: exclude quarks w/ mq < Lambda
        alphaq = self.alphaS(mq,nfq)
        if alphaq < 0: return 0
        else: return self.mhat(mq,nfq)*math.pow(self.alphaS(Q,nf),self.cb(nf))*self.one_c_cp_bp_b(nf)

    # operation
    def run(self,quark,nfq,scale,nf):
        # run to specified scale and nf
        return self.m(quark.mass,nfq,scale,nf)

class quarklist(object):
    def __init__(self):
        # mass-ordered
        self.qlist = [
            quark(2,0.0023,0.67), # up
            quark(1,0.0048,0.33), # down
            quark(3,0.095,0.33),  # strange
            quark(4,1.275,0.67),  # charm
            quark(5,4.18,0.33),   # bottom
        ]
        self.scale = None
        self.runner = massRunner()

    def set(self,scale):
        self.scale = scale
        # mask quarks above scale
        for q in self.qlist:
            # for decays
            if scale is None or 2*q.mass < scale: q.on = True
            else: q.on = False
            # for nf running
            if scale is None or q.mass < scale: q.active = True
            else: q.active = False
        # compute running masses
        if scale is not None:
            qtmp = self.get(active=True)
            nf = len(qtmp)
            for iq,q in enumerate(qtmp):
                q.massrun = self.runner.run(q,iq,scale,nf)
        # or undo running
        else:
            for q in self.qlist:
                q.massrun = q.mass

    def reset(self):
        self.set(None)

    def get(self,active=False):
        return [q for q in self.qlist if (q.active if active else q.on)]


class svjHelper(object):
    def __init__(self,svjl):
        with open(os.path.join(os.path.expandvars('$CMSSW_BASE'),'src/SVJ/Production/test/dict_xsec_Zprime.txt'),'r') as xfile:
            self.xsecs = {int(xline.split('\t')[0]): float(xline.split('\t')[1]) for xline in xfile}
        self.quarks_pseudo = quarklist()
        self.quarks_vector = quarklist()
        self.quarks = quarklist()
        
        self.alphaName = ""
        self.generate = None
        
        # parameters for lambda/alpha calculations

        self.n_c = 3
        self.n_f = 2
        self.b0 = 11.0/6.0*self.n_c - 2.0/6.0*self.n_f


    # Calculation of rho mass from Lattice QCD fits (arXiv:2203.09503v2)
    def calcLatticePrediction(self,mPiOverLambda,mPseudo):        
        mVectorOvermPseudo = (1.0/mPiOverLambda)*math.pow(5.76 + 1.5*math.pow(mPiOverLambda,2) ,0.5) 
        mVector = mVectorOvermPseudo*mPseudo
        return mVector


    # has to be "lambdaHV" because "lambda" is a keyword
    def setModel(self,channel,svjgamma,mMediator,mDark,mPseudo,mVector,rinv,alpha,mPiOverLambda,lambdaHV,BRGamma=None,yukawa=None,generate=True,boost=0.,boostvar=None,nMediator=None,sepproc=True):
        # check for issues
        if channel!="s" and channel!="t": raise ValueError("Unknown channel: "+channel)
        # store the basic parameters
        print("option svjgamma: ", svjgamma)
        self.channel = channel
        self.mg_name = "DMsimp_SVJ_s_spin1" if channel=="s" else "DMsimp_SVJ_t" if channel=="t" else ""
        self.generate = generate
        self.mMediator = mMediator
        self.svjgamma = svjgamma
        self.mPiOverLambda = mPiOverLambda
        self.lambdaHV = lambdaHV
        self.mPseudo = self.mPiOverLambda * self.lambdaHV
        self.mVector = self.calcLatticePrediction(self.mPiOverLambda,self.mPseudo)
        self.rinv = rinv
        self.brGamma = BRGamma

        self.nMediator = None
        self.yukawa = None
        self.sepproc = sepproc
        # yukawa not used by pythia "t-channel" generation (only includes strong pair prod)
        # but will still be included in name if provided in model setting
        if self.channel=="t":
            if nMediator<0: nMediator = None
            if nMediator is not None:
                if nMediator!=2 and generate:
                    raise ValueError("Pythia-only generation can only be used for pair production")
            self.nMediator = nMediator

            self.yukawa = yukawa
            if self.yukawa is None: raise ValueError("yukawa value must be provided for madgraph t-channel")

        # boosting
        allowed_boostvars = ["pt","madpt"]
        if boost>0 and boostvar is not None:
            if boostvar not in allowed_boostvars:
                raise ValueError("Unknown boost variable {}".format(boostvar))
            # some filters are implemented in madgraph
            if (boostvar=="madpt") and generate:
                raise ValueError("{} boostvar not compatible with Pythia-only generation".format(boostvar))
            self.boostvar = boostvar
            self.boost = boost
        else:
            self.boost = 0
            self.boostvar = ""

        # get more parameters
        self.xsec = self.getPythiaXsec(self.mMediator)
        self.mMin = self.mMediator-1
        self.mMax = self.mMediator+1
        self.mSqua = self.lambdaHV + 0.2 # dark scalar quark mass (also used for pTminFSR)

        # get limited set of quarks for decays (check mDark against quark masses, compute running)
        self.quarks_pseudo.set(self.mPseudo)
        self.quarks_vector.set(self.mVector)
            

    def getOutName(self,events=0,signal=True,outpre="outpre",part=None,sanitize=False,gridpack=False):
            _outname = outpre
            if signal:
                params = [
                    ("channel", "{}-channel".format(self.channel)),
                ]
                if self.nMediator is not None: params.append(("nMediator", "nMed-{:g}".format(self.nMediator)))
                params.extend([
                ("mMediator", "mMed-{:g}".format(self.mMediator)),
                ("mDark", "mDark-{:g}".format(self.mPseudo)),
                ("rinv", "rinv-{:g}".format(self.rinv)),
                ("alpha", "alpha-{}".format(self.brGamma)),
                ])
                if self.yukawa is not None: _outname += "_yukawa-{:g}".format(self.yukawa)
                if self.boost>0: _outname += "_{}{:g}".format(self.boostvar.upper(),self.boost)
                not_for_gridpack = ["rinv","alpha"]
                if not self.sepproc: not_for_gridpack.append("nMediator")
                if self.boostvar=="pt": not_for_gridpack.append("boost")
                for pname, pval in params:
                    if gridpack and pname in not_for_gridpack: continue
                    _outname += "_"+pval

            if self.generate is not None:
                if self.generate:
                    _outname += "_13TeV-pythia8"
                else:
                    _outname += "_13TeV-madgraphMLM-pythia8"
            if events>0: _outname += "_n-{:g}".format(events)
            if part is not None:
                _outname += "_part-{:g}".format(part)
            if sanitize:
                _outname = _outname.replace("-","_").replace(".","p")

            return _outname
      

    # allow access to all xsecs
    def getPythiaXsec(self,mMediator):
        xsec = 1.0
        # todo: get t-channel cross sections
        if self.channel=="t": return xsec
        # a function of mMediator
        if mMediator in self.xsecs: xsec = self.xsecs[mMediator]
        return xsec

     #implementation of rinv
    def invisibleDecay(self,mesonID,dmID):
        lines = ['{:d}:oneChannel = 1 {:g} 0 {:d} -{:d}'.format(mesonID,self.rinv,dmID,dmID)]
        return lines

    #Here mass insertion for Z' mediated decays to quarks, while br_gamma is used to set the branching ratio to photons for ALP mediated decays 
    def pseudo_scalar_visibleDecay(self,type,mesonID):
        
        theQuarks = self.quarks_pseudo.get()

        if type=="ALPplusMassInsertion":
            denom = sum([q.massrun**2 for q in theQuarks])
            for q in theQuarks:
                q.bf = (1-self.brGamma)*(1.0-self.rinv)*(q.massrun**2)/denom
        else:
            raise ValueError("unknown visible decay type: "+type)
    

        if (type=="ALPplusMassInsertion"):
            lines = ['{:d}:addChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,q.bf,q.id,q.id) for q in theQuarks if q.bf>0]
            lines += ['{:d}:addChannel = 1 {:g} 91 22 22'.format(mesonID,self.brGamma*(1.0-self.rinv))]


        return lines

   
    ### Vector visible decays via Z' mediated decays
    def vector_visibleDecay(self,type,mesonID):
            theQuarks = self.quarks_vector.get()
           
            if type=="Democratic":
                bfQuarks = (1.0-self.rinv)/float(len(theQuarks))
                for iq,q in enumerate(theQuarks):
                    theQuarks[iq].bf = bfQuarks
            else:
                raise ValueError("unknown visible decay type: "+type)
            lines = ['{:d}:addChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,q.bf,q.id,q.id) for q in theQuarks if q.bf>0]
            return lines


    ### Internal decays rho to pi pi
    def vector_internalDecay(self,type,mesonID):
        br = 0 
        if type=="Internal_simple":
            if (mesonID == 4900113):
               br = 1.0  
               decay_prod_1 = 4900211
               decay_prod_2 = 4900211
            if (mesonID == 4900213):
               br = 1.0  
               decay_prod_1 = 4900111
               decay_prod_2 = 4900211
        # lines for decays to quarks                                                                                                                              
        lines_rho_to_pipi = [
             '{:d}:mayDecay=on'.format(mesonID),
             '{:d}:oneChannel = 1 {:g} 91 {:d} -{:d}'.format(mesonID,br,decay_prod_1,decay_prod_2),
        ]
        
        # lines for decays to leptons                                                                                                                              
        return lines_rho_to_pipi


    def getPythiaSettings(self):
        # todo: include safety/sanity checks

        lines_schan = [
            # parameters for leptophobic Z'
            '4900023:m0 = {:g}'.format(self.mMediator),
            '4900023:mMin = {:g}'.format(self.mMin),
            '4900023:mMax = {:g}'.format(self.mMax),
            '4900023:mWidth = 0.01',
            '4900023:oneChannel = 1 0.982 102 4900101 -4900101',
            # SM quark couplings needed to produce Zprime from pp initial state
            '4900023:addChannel = 1 0.003 102 1 -1',
            '4900023:addChannel = 1 0.003 102 2 -2',
            '4900023:addChannel = 1 0.003 102 3 -3',
            '4900023:addChannel = 1 0.003 102 4 -4',
            '4900023:addChannel = 1 0.003 102 5 -5',
            '4900023:addChannel = 1 0.003 102 6 -6',
            # decouple
            '4900001:m0 = 50000',
            '4900002:m0 = 50000',
            '4900003:m0 = 50000',
            '4900004:m0 = 50000',
            '4900005:m0 = 50000',
            '4900006:m0 = 50000',
            '4900011:m0 = 50000',
            '4900012:m0 = 50000',
            '4900013:m0 = 50000',
            '4900014:m0 = 50000',
            '4900015:m0 = 50000',
            '4900016:m0 = 50000',
        ]

        # parameters for bifundamental mediators
        # (keep default flavor-diagonal couplings)
        bifunds = [4900001,4900002,4900003,4900004,4900005,4900006]
        lines_tchan = []
        for bifund in bifunds:
            lines_tchan.extend([
                '{:d}:m0 = {:g}'.format(bifund,self.mMediator),
                '{:d}:mMin = {:g}'.format(bifund,self.mMin),
                '{:d}:mMax = {:g}'.format(bifund,self.mMax),
                '{:d}:mWidth = 0.01'.format(bifund),
            ])
        lines_tchan.extend([
            # decouple
            '4900011:m0 = 50000',
            '4900012:m0 = 50000',
            '4900013:m0 = 50000',
            '4900014:m0 = 50000',
            '4900015:m0 = 50000',
            '4900016:m0 = 50000',
            '4900023:m0 = 50000',
        ])

        if self.generate:
            lines_schan.extend([
                'HiddenValley:ffbar2Zv = on',
            ])
            # pythia can only generate pair prod of bifundamental
            lines_tchan.extend([
                'HiddenValley:gg2DvDvbar = on',
                'HiddenValley:gg2UvUvbar = on',
                'HiddenValley:gg2SvSvbar = on',
                'HiddenValley:gg2CvCvbar = on',
                'HiddenValley:gg2BvBvbar = on',
                'HiddenValley:gg2TvTvbar = on',
                'HiddenValley:qqbar2DvDvbar = on',
                'HiddenValley:qqbar2UvUvbar = on',
                'HiddenValley:qqbar2SvSvbar = on',
                'HiddenValley:qqbar2CvCvbar = on',
                'HiddenValley:qqbar2BvBvbar = on',
                'HiddenValley:qqbar2TvTvbar = on',
            ])

        #Setting hidden spectrum    
        # hidden spectrum:                                                                                                                                                                                        
        # fermionic dark quark,                                                                                                                                                                                           
        # diagonal pseudoscalar meson, off-diagonal pseudoscalar meson, DM stand-in particle,                                                                                                                             
        # diagonal vector meson, off-diagonal vector meson, DM stand-in particle
        lines_decay = [
                '4900101:m0 = {:g}'.format(self.mSqua),
                '4900111:m0 = {:g}'.format(self.mPseudo),
                '4900211:m0 = {:g}'.format(self.mPseudo),
                '51:m0 = 0.0',
                '51:isResonance = false',
                '4900113:m0 = {:g}'.format(self.mVector),
                '4900213:m0 = {:g}'.format(self.mVector),
                '53:m0 = 0.0',
                '53:isResonance = false',
        ]    


            # other HV params
        lines_decay.extend([
            'HiddenValley:Ngauge = {:d}'.format(self.n_c),
            # when Fv has spin 0, qv spin fixed at 1/2
            'HiddenValley:spinFv = 0',
            'HiddenValley:FSR = on',
            'HiddenValley:fragment = on',
            'HiddenValley:alphaOrder = 1',
            'HiddenValley:Lambda = {:g}'.format(self.lambdaHV),
            'HiddenValley:nFlav = {:d}'.format(self.n_f),
            'HiddenValley:probVector = 0.75',
        ])



        # branching - effective rinv (applies to all meson species)
        # pseudoscalars have mass insertion decays (from Z') and ALP decays
        # vectors have democratic decays from Z' or internal decays to pi pi
      
        #all pseudoscalars decays
        lines_decay += self.invisibleDecay(4900111,51)
        lines_decay += self.pseudo_scalar_visibleDecay("ALPplusMassInsertion",4900111)
        lines_decay += self.invisibleDecay(4900211,51)
        lines_decay += self.pseudo_scalar_visibleDecay("ALPplusMassInsertion",4900211)
        
        if self.mPiOverLambda <= 1.5:
            #all vector mesons decays
            lines_decay += self.vector_internalDecay("Internal_simple",4900113)
            lines_decay += self.vector_internalDecay("Internal_simple",4900213)
        else:
            lines_decay += self.invisibleDecay(4900113,53)
            lines_decay += self.vector_visibleDecay("Democratic",4900113)
            lines_decay += self.invisibleDecay(4900213,53)
            lines_decay += self.vector_visibleDecay("Democratic",4900213)
        

        lines = []
        if self.channel=="s": lines = lines_schan + lines_decay
        elif self.channel=="t": lines = lines_tchan + lines_decay

        return lines

    def getJetMatchSettings(self):
        lines = [
            'JetMatching:setMad = off', # if 'on', merging parameters are set according to LHE file
            'JetMatching:scheme = 1', # 1 = scheme inspired by Madgraph matching code
            'JetMatching:merge = on', # master switch to activate parton-jet matching. when off, all external events accepted
            'JetMatching:jetAlgorithm = 2', # 2 = SlowJet clustering
            'JetMatching:etaJetMax = 5.', # max eta of any jet
            'JetMatching:coneRadius = 1.0', # gives the jet R parameter
            'JetMatching:slowJetPower = 1', # -1 = anti-kT algo, 1 = kT algo. Only kT w/ SlowJet is supported for MadGraph-style matching
            'JetMatching:qCut = 125.', # this is the actual merging scale. should be roughly equal to xqcut in MadGraph
            'JetMatching:nJetMax = 2', # number of partons in born matrix element for highest multiplicity
            'JetMatching:doShowerKt = off', # off for MLM matching, turn on for shower-kT matching
        ]

        return lines

    def getMadGraphCards(self,base_dir,lhaid,events=1,cores=1):
        if base_dir[-1]!='/': base_dir = base_dir+'/'

        # helper for templates
        def fill_template(inname, outname=None, **kwargs):
            if outname is None: outname = inname
            with open(inname,'r') as temp:
                old_lines = Template(temp.read())
                new_lines = old_lines.substitute(**kwargs)
            with open(inname,'w') as temp:
                temp.write(new_lines)
            if inname!=outname:
                shutil.move(inname,outname)

        mg_model_dir = os.path.expandvars(base_dir+"mg_model_templates")

        # replace parameters in relevant file
        param_args = dict(
            mediator_mass = "{:g}".format(self.mMediator),
            dark_quark_mass = "{:g}".format(self.mSqua),
        )
        if self.yukawa is not None: param_args["dark_yukawa"] = "{:g}".format(self.yukawa)
        fill_template(
            os.path.join(mg_model_dir,"parameters.py"),
            **param_args
        )

        # use parameters to generate card
        sys.path.append(mg_model_dir)
        from write_param_card import ParamCardWriter
        param_card_file = os.path.join(mg_model_dir,"param_card.dat")
        ParamCardWriter(param_card_file, generic=True)

        mg_input_dir = os.path.expandvars(base_dir+"mg_input_templates")
        modname = self.getOutName(events=events,outpre="SVJ",sanitize=True,gridpack=True)
        template_paths = [p for ftype in ["dat","patch"] for p in glob(os.path.join(mg_input_dir, "*."+ftype))]
        for template in template_paths:
            fname_orig = os.path.join(mg_input_dir,template)
            fname_new = os.path.join(mg_input_dir,template.replace("modelname",modname))
            fill_template(
                fname_orig,
                fname_new,
                modelName = modname,
                totalEvents = "{:g}".format(events),
                cores = "{:g}".format(cores),
                lhaid = "{:g}".format(lhaid),
                # for boosted
                madpt = "{:g}".format(self.boost if self.boostvar=="madpt" else 0.),
                # for t-channel
                procInclusive = "" if not self.sepproc else "#",
                procPair = "" if self.sepproc and self.nMediator==2 else "#",
                procSingle = "" if self.sepproc and self.nMediator==1 else "#",
                procNonresonant = "" if self.sepproc and self.nMediator==0 else "#",
            )

        return mg_model_dir, mg_input_dir
