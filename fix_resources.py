import os
import glob
import re

pattern = re.compile(
    r'[ \t]*resources:\s*'
    r'slurm_extra\s*=\s*f"\'--qos=\{config\.resources\.([a-zA-Z0-9_]+)\.qos\}\'",\s*'
    r'cpus_per_task\s*=\s*config\.resources\.([a-zA-Z0-9_]+)\.cpus,\s*'
    r'mem\s*=\s*config\.resources\.([a-zA-Z0-9_]+)\.mem(?:,)?\s*'
    r'(?:runtime\s*=\s*config\.resources\.([a-zA-Z0-9_]+)\.time)?',
    re.MULTILINE
)

def repl(match):
    # Keep the initial resources: tag, but we need to reconstruct the indentation
    # Let's just do a simpler line-by-line replacement.
    pass

