
localrules: download_omark_db

rule download_omark_db:
    output:
        os.path.join(dir.tools_omark,f"{config.evaluation.omark_db}.h5")
    params:
        db=config.evaluation.omark_db
    log:
        os.path.join(dir.logs,"download_omark_db.log")
    shell:
        """
        wget https://omabrowser.org/All/{params.db}.h5 -O {output} &> {log}
        """