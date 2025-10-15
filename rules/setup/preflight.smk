
rule input_check:
    input:
        config = "config.yaml"
    output:
        touch("input_checked.txt")
    script:
        "../scripts/input_check.py {input.config}"
