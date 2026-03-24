
localrules: install_sqanti, install_tama, download_omark_db

rule install_sqanti:
    output:
        touch(os.path.join(dir.tools_sqanti,"sqanti_installed.done"))
    conda:
        f"{dir.envs}/git.yaml"    
    shell:
        """
        git clone https://github.com/ConesaLab/SQANTI3.git {dir.tools_sqanti}
        """

rule install_gaqet2:
    output:
        touch(os.path.join(dir.tools_gaqet2,"gaqet2_installed.done"))
    conda:
        f"{dir.envs}/git.yaml"
    shell:
        """
        if [ -d {dir.tools_gaqet2} ]; then rm -rf {dir.tools_gaqet2}; fi
        git clone https://github.com/victorgcb1987/GAQET2.git {dir.tools_gaqet2}
        """

rule download_omark_db:
    output:
        os.path.join(dir.tools_omark,f"{config.qc.omark_db}.h5")
    params:
        db=config.qc.omark_db
    log:
        os.path.join(dir.logs,"download_omark_db.log")
    shell:
        """
        wget https://omabrowser.org/All/{params.db}.h5 -O {output} &> {log}
        """