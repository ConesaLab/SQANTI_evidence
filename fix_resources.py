import re
import glob

for f in glob.glob("rules/*.smk"):
    with open(f, "r") as file:
        content = file.read()
    
    content = content.replace("}}\'", "}\'")
    
    with open(f, "w") as file:
        file.write(content)
