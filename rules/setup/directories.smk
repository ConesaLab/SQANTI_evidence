"""
Ensures the correct output paths in the final direcotry
"""
import attrmap as ap

dir = ap.AttrMap()

### BUSCO database

try:
    assert (ap.utils.to_dict(config.args)["busco_db"]) is not None
    dir.busco_db = str(config.args.busco_db)
except (KeyError, AssertionError):
    dir.busco_db = os.path.join(workflow.basedir, "..", "databases", "busco")

### Output_location
try:
    assert config.required.outdir is not None
    dir.out.base = str(config.required.outdir)
except (KeyError, AssertionError):
    dir.out.base = "evidence_annot"

### Workflow_dirs
dir.envs = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

### Tools directories
dir.tools = str(config.required.toolsdir)
dir.tools_conda = os.path.join(dir.tools, "conda_envs")
dir.tools_tama = os.path.join(dir.tools, "tama")
dir.tools_db = os.path.join(dir.tools, "databases")
dir.tools_omark = os.path.join(dir.tools_db, "omark")
dir.tools_busco = os.path.join(dir.tools_db, "busco")
dir.tools_sqanti = os.path.join(dir.tools, "sqanti")
dir.tools_gaqet2 = os.path.join(dir.tools, "gaqet2")
dir.tools_index = os.path.join(dir.tools, "minimap2")
dir.tools_reference = os.path.join(dir.tools, "split_reference")


### Output_dirs
dir.logs = os.path.join(dir.out.base, "logs")
## Ab_initio
dir.out.ab_initio = os.path.join(dir.out.base, "ab_initio")
dir.out.ab_busco = os.path.join(dir.out.ab_initio, "busco")
dir.out.ab_augustus = os.path.join(dir.out.ab_initio, "augustus")
dir.out.ab_augustus_model = os.path.join(dir.out.ab_augustus, "model")
dir.out.ab_augustus_training = os.path.join(dir.out.ab_augustus, "training")

## IsoSeq3
dir.out.isoseq = os.path.join(dir.out.base, "isoseq")
dir.out.isoseq_lima = os.path.join(dir.out.isoseq, "lima")
dir.out.isoseq_refine = os.path.join(dir.out.isoseq, "refine")
dir.out.isoseq_cluster = os.path.join(dir.out.isoseq, "cluster")
dir.out.isoseq_mapping = os.path.join(dir.out.isoseq, "mapping")
dir.out.isoseq_collapsed = os.path.join(dir.out.isoseq, "collapsed")

## Hints
dir.out.evidence_driven = os.path.join(dir.out.base, "evidence_driven")
dir.out.ed_sqanti = os.path.join(dir.out.evidence_driven, "sqanti")
dir.out.ed_hints = os.path.join(dir.out.evidence_driven, "hints")
dir.out.ed_augustus = os.path.join(dir.out.evidence_driven, "augustus")

# Quality control
dir.out.qc = os.path.join(dir.out.base, "quality_control")
dir.out.qc_omark = os.path.join(dir.out.qc, "omark")
dir.out.qc_busco = os.path.join(dir.out.qc, "busco")
dir.out.qc_agat = os.path.join(dir.out.qc, "agat")
dir.out.qc_gaqet2 = os.path.join(dir.out.qc, "gaqet2")