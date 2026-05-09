import statistics
with open ("file.fq") as f: 
    lines = f.read().splitlines()
quality_line = [3] 
qualities = [ord(x) - 33 for x in quality_line]
print(statistics.median(qualities))
#ord -33 zmienia znak jakości na wartość phred

