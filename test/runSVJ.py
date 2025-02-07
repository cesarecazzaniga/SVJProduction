import FWCore.ParameterSet.Config as cms
import sys, os
from SVJ.Production.optSVJ import options, _helper

def replace_last(item,old,new):
    return new.join(item.rsplit(old,1))

# output name definition
_outname = _helper.getOutName(
    options.maxEvents,
    part=options.part,
    signal=options.signal and len(options.scan)==0 and len(options.fragment)==0,
    outpre="outpre",
)
_outname += ".root"

_inname = ""
if len(options.inpre)>0:
    _inname = _outname.replace("outpre",options.inpre)
    if options.maxEvents!=options.maxEventsIn: _inname = _inname.replace("_n-{}_".format(options.maxEvents),"_n-{}_".format(options.maxEventsIn),1)

def fix_inname(inname,options,lhe=False):
    if len(options.indir)>0: inname = options.indir+"/"+inname
    if len(options.redir)>0 and inname.startswith("/store"): inname = options.redir+inname
    if not lhe:
        if not inname.startswith("/store") and not inname.startswith("root:"): inname = "file:"+inname
        # folderization only occurs in remote case
        elif options.useFolders:
            inname = replace_last(inname,"_part","/part")
    return inname

# import process
process = getattr(__import__(options.config,fromlist=["process"]),"process")

# input settings
process.maxEvents.input = cms.untracked.int32(options.maxEvents)
if len(_inname)>0:
    if hasattr(process.source,"fileNames"):
        _inname = fix_inname(_inname,options)
        process.source.fileNames = cms.untracked.vstring(_inname)
    elif hasattr(process,"externalLHEProducer"):
        _inname = _helper.getOutName(events=options.maxEventsIn,outpre=options.inpre,gridpack=True)+".tar.xz"
        _inname = fix_inname(_inname,options,True)
        process.externalLHEProducer.args = cms.vstring(_inname)
        process.externalLHEProducer.nEvents = cms.untracked.uint32(options.maxEvents)
        if options.syst: process.externalLHEProducer.scriptName = cms.FileInPath("SVJ/Production/test/run_svj_tarball_syst.sh")
if process.source.type_()=='EmptySource':
    process.source.firstEvent = cms.untracked.uint32((options.part-1)*options.maxEvents+1)
    if len(options.scan)>0: process.source.numberEventsInLuminosityBlock = cms.untracked.uint32(100)

# output settings
oprocess = process if (not hasattr(process,'subProcesses') or len(process.subProcesses)==0) else process.subProcesses[-1].process()
if len(options.output)==0: options.output = sorted(oprocess.outputModules_())
if len(options._outpre)!=len(options.output):
    raise ValueError("Mismatch between # of output prefixes and # of output modules\n\tOutput prefixes are: "+", ".join(options._outpre)+"\n\tOutput modules are: "+", ".join(options.output))
for iout,output in enumerate(options.output):
    if len(output)==0: continue
    if not hasattr(oprocess,output):
        raise ValueError("Unavailable output module: "+output)
    getattr(oprocess,output).fileName = 'file:'+_outname.replace("outpre",options._outpre[iout])

# reset all random numbers to ensure statistically distinct but reproducible jobs
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
# CHANGED from maxEvents+part because concurrent LHE generation adds thread number to default random seed -> degeneracy
randHelper.resetSeeds(int(str(options.part)+str(options.maxEvents)))

if options.signal:
    if len(options.scan)>0:
        if hasattr(process,'generator'):
            process.generator = getattr(__import__("SVJ.Production."+options.scan+"_cff",fromlist=["generator"]),"generator")
    elif len(options.fragment)>0:
        if hasattr(process,'generator'):
            _params = getattr(__import__("SVJ.Production."+options.fragment,fromlist=["processParameters"]),"processParameters")
            if isinstance(_params,tuple): _params = _params[0] # handle weird python behavior change
            process.generator.PythiaParameters.processParameters = _params
    else:
        # generator settings
        if hasattr(process,'generator'):
            process.generator.crossSection = cms.untracked.double(_helper.xsec)
            process.generator.PythiaParameters.processParameters = cms.vstring(_helper.getPythiaSettings())
            if hasattr(process.generator.PythiaParameters,"JetMatchingParameters"):
                process.generator.PythiaParameters.JetMatchingParameters = cms.vstring(_helper.getJetMatchSettings())
            if options.suep:
                process.generator.UserCustomization = cms.VPSet(
                    _helper.getHookSettings()
                )
                if options.filterHT > 0:
                    process.genHTFilter = cms.EDFilter("GenHTFilter",
                        src = cms.InputTag("ak4GenJetsNoNu"), #GenJet collection as input
                        jetPtCut = cms.double(30.0), #GenJet pT cut for HT
                        jetEtaCut = cms.double(2.5), #GenJet eta cut for HT
                        genHTcut = cms.double(options.filterHT) #genHT cut
                    )
                    if hasattr(process,'ProductionFilterSequence'):
                        process.ProductionFilterSequence += process.pgen
                        process.ProductionFilterSequence += process.genHTFilter
    # gen filter settings
    # pythia implementation of model has 4900111/211 -> -51 51 and 4900113/213 -> -53 53
    # this is a stand-in for direct production of a single stable dark meson in the hadronization
    # stable mesons should be produced in pairs (Z2 symmetry),
    # so require total number produced by pythia to be a multiple of 4
    # do *not* require this separately for 111/211 and 113/213 (pseudoscalar vs. vector)
    if options.filterZ2 and hasattr(process,'ProductionFilterSequence'):
        process.darkhadronZ2filter = cms.EDFilter("MCParticleModuloFilter",
            moduleLabel = cms.InputTag('generator','unsmeared'),
            particleIDs = cms.vint32(51,53),
            multipleOf = cms.uint32(4),
            absID = cms.bool(True),
        )
        process.ProductionFilterSequence += process.darkhadronZ2filter

    # also filter out events with Zprime -> SM quarks
    if options.channel=="s" and hasattr(process,'ProductionFilterSequence'):
        process.darkquarkFilter = cms.EDFilter("MCParticleModuloFilter",
            moduleLabel = cms.InputTag('generator','unsmeared'),
            particleIDs = cms.vint32(4900101),
            multipleOf = cms.uint32(2),
            absID = cms.bool(True),
            min = cms.uint32(2),
            status = cms.int32(23),
        )
        process.ProductionFilterSequence += process.darkquarkFilter

    if options.boost>0 and hasattr(process,'ProductionFilterSequence'):
        # apply GenJet pt cut for boosted search
        if options.boostvar=="pt":
            process.genjetptFilter = cms.EDFilter("GenJetPTFilter",
                ptMin = cms.double(options.boost),
            )
            process.ProductionFilterSequence += process.pgen
            process.ProductionFilterSequence += process.genjetptFilter

    if not options.sepproc and options.nMediator>=0 and hasattr(process,'ProductionFilterSequence'):
        # apply LHE-level filter *before* pythia
        process.nmedfilter = cms.EDFilter("LHENParticleFilter",
            min = cms.int32(options.nMediator),
            max = cms.int32(options.nMediator),
            particleIDs = cms.vint32(4900001,4900002,4900003,4900004,4900005,4900006),
        )
        process.ProductionFilterSequence.insert(0,process.nmedfilter)

def set_if_has(obj,attr,val):
    if hasattr(obj,attr): setattr(obj,attr,val)

if hasattr(process,'generator'):
    if options.quiet:
        set_if_has(process.generator,'maxEventsToPrint',0)
        set_if_has(process.generator,'pythiaHepMCVerbosity',False)
        set_if_has(process.generator,'pythiaPylistVerbosity',0)
        if hasattr(process.generator,'PythiaParameters'):
            pythia_quiet_settings = [
                'Print:quiet = on',
                'Init:showProcesses = off',
                'Init:showMultipartonInteractions = off',
                'Init:showChangedSettings = off',
                'init:showChangedParticleData = off',
            ]
            set_if_has(
                process.generator.PythiaParameters,
                'pythia8CommonSettings',
                getattr(process.generator.PythiaParameters,'pythia8CommonSettings',[])+pythia_quiet_settings
            )
    elif hasattr(process.generator,'maxEventsToPrint'):
        process.generator.maxEventsToPrint = options.printEvents

if options.quiet and hasattr(process,'MessageLogger'):
    for dest in process.MessageLogger.destinations:
        dest_attr = getattr(process.MessageLogger,dest)
        for level in ['INFO','WARNING']:
            existing_pset = getattr(dest_attr,level,cms.untracked.PSet())
            existing_pset.limit = cms.untracked.int32(0)
            setattr(dest_attr,level,existing_pset)

# genjet/met settings - treat DM stand-ins as invisible
_particles = ["genParticlesForJetsNoMuNoNu","genParticlesForJetsNoNu","genCandidatesForMET","genParticlesForMETAllVisible"]
for _prod in _particles:
    if hasattr(process,_prod):
        getattr(process,_prod).ignoreParticleIDs.extend([51,52,53])
if hasattr(process,'recoGenJets') and hasattr(process,'recoAllGenJetsNoNu'):
    process.recoGenJets += process.recoAllGenJetsNoNu
    # to get hadronFlavour at gen level
    process.load("PhysicsTools.PatAlgos.slimming.genParticles_cff")
    process.recoGenJets += process.prunedGenParticlesWithStatusOne
    process.recoGenJets += process.prunedGenParticles
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    process.recoGenJets += process.selectedHadronsAndPartonsForGenJetsFlavourInfos
    process.load("PhysicsTools.JetMCAlgos.AK4GenJetFlavourInfos_cfi")
    process.ak8GenJetFlavourInfos = process.ak4GenJetFlavourInfos.clone(
        jets = "ak8GenJetsNoNu",
        rParam = cms.double(0.8),
    )
    process.recoGenJets += process.ak4GenJetFlavourInfos
    process.recoGenJets += process.ak8GenJetFlavourInfos
    for output in options.output:
        if len(output)==0: continue
        output_attr = getattr(oprocess,output)
        if hasattr(output_attr,"outputCommands"):
            output_attr.outputCommands.extend([
                'keep *_ak4GenJetFlavourInfos_*_*',
                'keep *_ak8GenJetFlavourInfos_*_*',
            ])
if hasattr(process,'genJetParticles') and hasattr(process,'genParticlesForJetsNoNu'):
    process.genJetParticles += process.genParticlesForJetsNoNu
    for output in options.output:
        if len(output)==0: continue
        output_attr = getattr(oprocess,output)
        if hasattr(output_attr,"outputCommands"):
            output_attr.outputCommands.extend([
                'keep *_genParticlesForJets_*_*',
                'keep *_genParticlesForJetsNoNu_*_*',
                'keep *_bsmHtFilter_*_*',
            ])

# DIGI settings
if hasattr(process,"mixData"):
    if options.year.startswith("2016"): puname = "Neutrino_E-10_gun_RunIISummer20ULPrePremix-UL16_106X_mcRun2_asymptotic_v13-v1_PREMIX.pkl"
    elif options.year=="2017": puname = "Neutrino_E-10_gun_RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3_PREMIX.pkl"
    elif options.year=="2018": puname = "Neutrino_E-10_gun_RunIISummer20ULPrePremix-UL18_106X_upgrade2018_realistic_v11_L1v1-v2_PREMIX.pkl"
    if not os.path.isfile(puname):
        print "retrieving "+puname
        os.system("xrdcp -f root://cmseos.fnal.gov//store/user/pedrok/SVJ2017/pileup/"+puname+" .")
        if not os.path.isfile(puname):
            raise Exception("Could not retrieve pileup input list.")
    import cPickle as pickle
    process.mixData.input.fileNames = cms.untracked.vstring(*pickle.load(open(puname,"rb")))

# miniAOD settings
_pruned = ["prunedGenParticlesWithStatusOne","prunedGenParticles"]
_keeps = ["keep (4900001 <= abs(pdgId) <= 4900991 )", "keep (51 <= abs(pdgId) <= 53)"]
if options.suep: 
    _keeps = ["keep 999998 <= abs(pdgId) <= 999999", "++keep  abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 1 || abs(pdgId) == 211", "keep++ abs(pdgId) == 1" ]
    # keep dark pions, darkphotons 
    # higgs already kept
    # keep SM decay products, electrons, muons, pions, uubar
    # keep decays of uubar
for _prod in _pruned:
    if hasattr(process,_prod):
        # keep HV & DM particles
        getattr(process,_prod).select.extend(_keeps)

def add_outputs(output_list):
    if not isinstance(output_list,list): output_list = [output_list]
    for output in options.output:
        if len(output)==0: continue
        output_attr = getattr(oprocess,output)
        if hasattr(output_attr,"outputCommands"):
            output_attr.outputCommands.extend(output_list)

if options.scout and "MINIAOD" in options.config:
    add_outputs([
        'keep *_hltScouting*_*_*',
    ])

if options.hepmc and any(cfg in options.config for cfg in ["LHE","GEN","SIM","DIGI","HLT","RECO","MINI"]):
    add_outputs([
        'keep *_generator_unsmeared_*',
    ])

# multithreading options
if options.threads>0:
    if not hasattr(process,"options"):
        process.options = cms.untracked.PSet()
    process.options.numberOfThreads = cms.untracked.uint32(options.threads)
    process.options.numberOfStreams = cms.untracked.uint32(options.streams if options.streams>0 else 0)

if options.tmi:
    from Validation.Performance.TimeMemoryInfo import customise
    process = customise(process)
    
if options.dump:
    print process.dumpPython()
    sys.exit(0)
