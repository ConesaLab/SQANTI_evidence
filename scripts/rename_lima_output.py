import os,sys
import logging

# Use the pipeline logger if available, otherwise create a basic one
logger = logging.getLogger('pipeline')
if not logger.handlers:
    # Fallback for standalone execution
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

def main():
    dir=os.path.dirname(snakemake.input[0])
    samples_file = snakemake.params[0]
    brk_2_samples = {}

    with open(samples_file) as f:
        samples = f.readlines()
        for line in samples:
            line = line.strip("\n")
            brk_2_samples[line.split(",")[0]] = line.split(",")[1]
    
    for brk in os.listdir(dir):
        directory = os.path.join(dir, brk)
        if os.path.isdir(directory):
            try:
                sample_name = brk_2_samples[brk]
                path=os.path.join(dir, sample_name)
                try:
                    os.mkdir(path)
                except:
                    pass
                for file in os.listdir(directory):
                    new_name = file.replace(f"fl.{brk}", f"{sample_name}.fl")
                    logger.info(f"Renaming {file} to {new_name}")
                    logger.debug(f"FULL PATH: {os.path.join(directory, file)}, to {os.path.join(path,new_name)}")
                    os.rename(os.path.join(directory, file), 
                            os.path.join(path,new_name))
                os.rmdir(directory)
            except KeyError:
                if brk not in list(brk_2_samples.values()):
                    logger.error(f"Sample {brk} not found in samples file")
                    sys.exit(1)
            
                


if __name__=="__main__":
    main()
    # with open(snakemake.output[0], "w") as f:
    #     f.write("Lima has been correctly renamed\n")
    # f.close()