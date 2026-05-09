#odbliczanie mediany
import statistics
with open ("file.fq") as f: 
    lines = f.read().splitlines()
quality_line = [3] 
qualities = [ord(x) - 33 for x in quality_line]
print(statistics.median(qualities))
#ord -33 zmienia znak jakości na wartość phred
#spodziewana liczba błędów
first_20 = qualities[:20]
last_20 = qualities[-20:]
errors_first = sum([10**(-q/10) for q in first_20])
errors_last = sum([10**(-q/10) for q in last_20])
print(f"Spodziewana liczba błędów w pierwszych 20 znakach: {errors_first}")

#wyciąganie barkodów
barcode = header.split("#")[1].split("/")[0] #1 bierze część po #, split("/")[0] bierze część przed / usuwa /1 i zostawia sam indeks

