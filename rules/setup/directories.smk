"""
Ensures the correct output paths in the final direcotry
"""
import attrmap as ap

dir = ap.AttrMap()

### BUSCO database

try:
    assert (ap.utils.to_dict(config.args)["busco_db"]) is not None
    dir.busco_db = str(config.args.busco_db)
except (KeyError, AssertionError, AttributeError):
    dir.busco_db = os.path.join(workflow.basedir, "..", "databases", "busco")

### Output_location
dir.out.base = str(config.project.outdir)

### Workflow_dirs
dir.envs = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

### Tools directories
dir.tools = str(config.project.toolsdir)
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

## IsoQuant
dir.out.isoquant = os.path.join(dir.out.base, "isoquant")

## Hints
dir.out.evidence_driven = os.path.join(dir.out.base, "evidence_driven")
dir.out.ed_sqanti = os.path.join(dir.out.evidence_driven, "sqanti")
dir.out.ed_hints = os.path.join(dir.out.evidence_driven, "hints")
dir.out.ed_augustus = os.path.join(dir.out.evidence_driven, "augustus")

# Quality control
dir.out.qc = os.path.join(dir.out.base, "quality_control")
dir.out.qc_gffcompare = os.path.join(dir.out.qc, "gffcompare")
dir.out.qc_gaqet2 = os.path.join(dir.out.qc, "gaqet2")
