#odbliczanie mediany
import statistics
with open ("file.fq") as f: 
    lines = f.read().splitlines()
quality_line = [3] 
qualities = [ord(x) - 33 for x in quality_line]
print(statistics.median(qualities))
#ord -33 zmienia znak jakości na wartość phred
#spodziewana liczba błędów to suma prawdopodobieństw błędu dla każdej zasady, gdzie prawdopodobieństwo błędu to 10^(-Q/10), Q to wartość jakościowa dla danej zasady. Im wyższe Q, tym mniejsze prawdopodobieństwo błędu. Dlatego sumujemy te wartości dla pierwszych 20 znaków i ostatnich 20 znaków, aby oszacować liczbę błędów w tych częściach sekwencji.
first_20 = qualities[:20]
last_20 = qualities[-20:]
errors_first = sum([10**(-q/10) for q in first_20])
errors_last = sum([10**(-q/10) for q in last_20])
print(f"Spodziewana liczba błędów w pierwszych 20 znakach: {errors_first}")

#wyciąganie barkodów
barcode = header.split("#")[1].split("/")[0] #1 bierze część po #, split("/")[0] bierze część przed / usuwa /1 i zostawia sam indeks

#PHRED - każda wartość ma przypisaną wartość jakościową, phred, mówi ona o tym, jak bardzo sekwenator ufa, że dana zasada została odczytana poprawnie. Im wyższe Q tym mniejsze prawdopodoieństwo błędu
#QC wykonuje się przed mapowaniem, analizą ekspresji, CHIP-seq i variant callingiem, żeby ocenić jakość danych i zdecydować czy można je dalej analizować, czy trzeba je odrzucić lub poprawić.
#Typowe błędy w FastQC (pozwala wykrywać błędy przed dalszą analizą) to: niska jakość na końcach sekwencji, obecność adapterów, zanieczyszczenia, nierównomierna dystrybucja jakości, duplikaty, itp. FastQC generuje raporty i wykresy, które pomagają zidentyfikować te problemy i podjąć odpowiednie działania naprawcze.
#adapter to techniczna sekwencja dodawana podczas przygotowywania biblioteki. Może utrudniać mapowanie i zaburzać analizę wariantów. Narzędzia: Cutadapt, trimmomatic, Fastp - cutadapt: cutadapt \  -a AGATCGGAAGAGC        adapter do usunięcia z końca 3' / -q 20                   przytnij zasady o jakości poniżej Q20 / -m 20                   odrzuć odczyty krótsze niż 20 nt / -o trimmed.fastq        plik wynikowy / input.fastq             plik wejściowy //  trimmomatic / SE                         tryb single-end / -phred33                   kodowanie jakości Phred+33 / ILLUMINACLIP               usuwanie adapterów / LEADING:20                 przytnij początek, jeśli Q < 20 / TRAILING:20                przytnij koniec, jeśli Q < 20 / SLIDINGWINDOW:4:20         okno 4 nt, średnia jakość min. 20 / MINLEN:20                  odrzuć odczyty krótsze niż 20 nt
#Fastp tworzy raport HTML - fastp \ -i input.fastq \ -o trimmed.fastq \ -q 20 \ -l 20 \ -h fastp_report.html \ -j fastp_report.json \\ -i       plik wejściowy / -o       plik wynikowy / -q 20    minimalna jakość zasad / -l 20    minimalna długość odczytu po przycięciu / -h       raport HTML / -j       raport JSON

#Najpierw robie fastqc potem np. fastp, potem znowu fastqc żeby sprawdzić czy jakość się poprawiła. FastQC generuje raporty i wykresy, które pomagają zidentyfikować te problemy i podjąć odpowiednie działania naprawcze. Fastp tworzy raport HTML - fastp \ -i input.fastq \ -o trimmed.fastq \ -q 20 \ -l 20 \ -h fastp_report.html \ -j fastp_report.json \\ -i       plik wejściowy / -o       plik wynikowy / -q 20    minimalna jakość zasad / -l 20    minimalna długość odczytu po przycięciu / -h       raport HTML / -j       raport JSON
#Demultipleksing to rozdzielanie odczytów na próbki na podstawie indeksów/barkodów. 

#Mapowanie to dopasowanie odczytów do znanej sekwencji, jeśli nie ma sekwencji referencyjnej to robię składanie de novo - łączysz odczyty w kontigi. Mapery używają heurystyk opartych o ziarna, czyli krótkie fragmenty odczytu. Najpierw wyszukuje się dokładne lub prawie dokładne dopasowanie a potem rozszerza do pełnego. Jeśli seed jest krótki można przypasować do tysięcy miejsc. 

#STRUKTURY INDEKSOWANIA SEKWENCJI: 
#   BWT + FM-index - BWT (Burrows-Wheeler Transform) to algorytm transformacji tekstu, który przekształca sekwencję w sposób, który ułatwia kompresję i wyszukiwanie. FM-index to struktura danych oparta na BWT, która umożliwia szybkie wyszukiwanie wzorców w dużych sekwencjach. Używane w mapperach takich jak Bowtie2, BWA, STAR, HISAT2. SZYBKO WYSZUKUJE KRÓTKIE SEKWENCJE W SKOMPRESOWANEJ REFERENCJI!!
#   Hierarchiczny FM-index oparty na grafach:  Używany przez HISAT2, pozwala na efektywne mapowanie do referencji z dużą zmiennością, np. genomów z wieloma haplotypami. JEST ZOPTYMALIZOWANY PRZEZ RNA-SEQ I POZWALA MAPOWAĆ ODCZYTY PRZEZ DŁUGIE INTRONY. 
#   Tablica sufiksów: Używana przez mappery takie jak STAR, pozwala na szybkie wyszukiwanie wzorców w dużych sekwencjach, ale jest bardziej pamięciochłonna niż BWT/FM-index. SZYBKO WYSZUKUJE DŁUGIE SEKWENCJE, ALE ZAJMUJE DUŻO PAMIĘCI.

#BWA-MEM do wgs/wes, czyli DNA-seq. Mapuje odczyty do genomu, ale nie jest zoptymalizowany pod splicing. 
#HISAT2 - najlepszy do RNA-seq, bo jest zoptymalizowany pod splicing. Pozwala mapować odczyty przez długie introny i połączenia egzon-egzon.

#MAPQ - jakość mapowania. MAPQ = -10 * log10(p), gdzie p to prawdopodobieństwo, że odczyt został źle zmapowany. Im wyższy MAPQ, tym większa pewność mapowania. MAPQ=0 oznacza, że odczyt może być zmapowany do wielu miejsc lub że mapper nie jest pewny mapowania. MAPQ=255 oznacza, że mapper nie obliczył jakości mapowania. W praktyce często używa się progu MAPQ > 30, aby wybrać odczyty o wysokiej jakości mapowania.

#ETAPY ANALIZY RNA-SEQ:
# 1. QC surowych odczytów
# 2. Usunięcie adapterów i filtrowanie jakościowe
# 3. Mapowanie do genomu/transkryptomu, np. HISAT2/STAR
# 4. Ocena mapowania
# 5. Zliczanie odczytów na geny, np. featureCounts
# 6. Normalizacja
# 7. Filtrowanie genów o niskiej ekspresji
# 8. Analiza różnicowej ekspresji, np. edgeR, DESeq2, limma-voom
# 9. Wizualizacja: MDS/PCA, heatmapy, MA/MD plot
# 10. Interpretacja biologiczna, np. analiza wzbogacenia genów

#dUTP pozwala zachować informację o kierunkowości transkryptu, czyli wiadomo z której nici pochodził RNA. 
#Cel analizy różnicowej ekspresji to znalezienie genów, których ekspresja różni się między warunkami. Metoda identyfikacji genów, których poziom ekspresji koreluje z określonymi zmiennymi biologicznymi. 
#CPM = liczba odczytów genu / liczba wszystkich zmapowanych odczytów × 1 000 000
#CPM koryguje głębokość sekwencjonowania, ale nie koryguje długości genu.
#Do porównania tych samych genów między warunkami wystarczy CPM, ponieważ długość genu nie zmienia się między próbkami. Do porównywania różnych genów w tej samej próbce potrzebna jest normalizacja względem długości, np. RPKM, ponieważ dłuższe geny zbierają więcej odczytów.
#geny z bardzo małą liczbą odczytów:
#    mają dużą zmienność losową
#   często są niewiarygodne statystycznie
#   zwiększają liczbę testów wielokrotnych
#   mogą pogarszać korekcję FDR

#WYKRES MDS - próbki z tej samej grupy biologicznej leżą blisko siebie, różne warunki są oddzielone - próbki grupują się według lane/batch, a nie według warunku biologicznego
#MD/MA - oś X : średnia ekspresja, oś Y : log fold change. dla genów o niskiej ekspresji wariancja jest większa, więc potrzeba większej zmiany ekspresji, aby wynik był istotny.

#ChIP-seq: identyfikacja miejsc w genomie, do których wiąże się konkretne białko albo gdzie występuje konkretna modyfikacja histonowa. Np. FOXA1 w linii HepG2
#Pipeline:
# 1. Fragmentacja chromatyny
# 2. Immunoprecypitacja przeciwciałem przeciwko badanemu białku
# 3. Oczyszczenie DNA
# 4. Sekwencjonowanie
# 5. Mapowanie odczytów
# 6. Peak calling
# 7. Analiza motywów i lokalizacji peaków

#Próba kontrolna, input - pomaga odróżnić prawdziwe wzbogacenie od tła technicznego. Uwzględnia bias fragmentacji, regiony łatwo sekwencjonujące się, regiony powtarzalne, niespecyficzne wiązanie, lokalne różnice dostępności chromatyny
#PEAK - region genomu, w którym jest lokalne wzbogacenie odczytów z eksperymentu ChIP względem tla. 

#Tag Directory to katalog tworzony przez HOMER, w którym program zapisuje przetworzone informacje o odczytach z pliku BAM. Zamiast za każdym razem analizować cały BAM od nowa, HOMER używa tego katalogu jako swojej wewnętrznej reprezentacji danych. Pierwszy krok anlizy ChIP-seq 
#np. makeTagDirectory treatment_tags chr1_treatment.bam
#użyć findPeaks, podając próbę kontrolną jako tło, a potem przekonwertować wynik do formatu BED.
#   findPeaks treatment_tags -style factor -i control_tags -o foxa1_peaks.txt
#   findPeaks treatment_tags -style factor -i control_tags -o foxa1_peaks.txt = Znajdź regiony, gdzie odczytów z FOXA1 jest istotnie więcej niż w kontroli.
#przekonwertować wynik do .bed, np. skryptem pos2bed.pl = 
#   pos2bed.pl foxa1_peaks.txt > foxa1_peaks.bed
#W materiałach jest polecenie, żeby użyć makeUCSCfile dla obu katalogów tagów i wygenerować pliki .bedGraph
#   makeUCSCfile treatment_tags -o treatment.bedGraph 
#INTERPRETACJA W IGV: Jeśli peak jest prawdziwy, w treatment powinien być wyraźny sygnał, a w control dużo słabszy albo brak. Jeśli podobny sygnał jest w treatment i control, może to być artefakt/tło.
#Do sprawdzenia, gdzie znajdują się peaki, używa się:
#    annotatePeaks.pl foxa1_peaks.bed hg38 > foxa1_peaks_annotated.txt
#użyć findMotifsGenome.pl na peakach BED, z genomem hg38, wykorzystując 10 wątków. Trzeba potem znaleźć dwa motywy de novo z najniższym p-value.
#   findMotifsGenome.pl foxa1_peaks.bed hg38 motifs_output -p 10
#W materiałach jest polecenie, żeby użyć annotatePeaks.pl z -m i -size, aby sprawdzić występowanie motywów w promieniu 100 bp od środka peaku.
#    annotatePeaks.pl foxa1_peaks.bed hg38 \
#    -m motifs_output/homerResults/motif1.motif \
#    -size 200 \
#    > peaks_with_motif1.txt
#W materiałach jest polecenie, żeby z -m i -hist zrobić histogram rozmieszczenia motywu wokół centrum peaku ±500 bp.
#   annotatePeaks.pl foxa1_peaks.bed hg38 \
#   -m motifs_output/homerResults/motif1.motif \
#   -size 1000 \
#   -hist 10 \
#   > motif1_histogram.txt
#Centralne położenie motywu względem centrum peaków sugeruje, że motyw odpowiada rzeczywistemu miejscu wiązania badanego czynnika transkrypcyjnego. Jeśli motyw jest rozproszony albo przesunięty, może oznaczać motyw kofaktora, szerszy region regulatorowy albo mniej precyzyjne peaki.

#CAŁY WORKFLOW HOMER:
# # 1. Aktywuj środowisko
# conda activate chip

# # 2. Zaindeksuj BAM-y
# samtools index chr1_treatment.bam
# samtools index chr1_control.bam

# # 3. Utwórz Tag Directories
# makeTagDirectory treatment_tags chr1_treatment.bam
# makeTagDirectory control_tags chr1_control.bam

# # 4. Sprawdź podstawowe informacje
# cat treatment_tags/tagInfo.txt
# cat control_tags/tagInfo.txt

# # 5. Wykryj peaki z kontrolą jako tłem
# findPeaks treatment_tags -style factor -i control_tags -o foxa1_peaks.txt

# # 6. Policz liczbę peaków
# grep -v "^#" foxa1_peaks.txt | wc -l

# # 7. Zamień peaki na BED
# pos2bed.pl foxa1_peaks.txt > foxa1_peaks.bed

# # 8. Wygeneruj bedGraph do IGV
# makeUCSCfile treatment_tags -o treatment.bedGraph
# makeUCSCfile control_tags -o control.bedGraph

# # 9. Zannotuj peaki
# annotatePeaks.pl foxa1_peaks.bed hg38 > foxa1_peaks_annotated.txt

# # 10. Szukaj motywów
# findMotifsGenome.pl foxa1_peaks.bed hg38 motifs_output -p 10

# # 11. Sprawdź występowanie konkretnego motywu w peakach ±100 bp
# annotatePeaks.pl foxa1_peaks.bed hg38 \
#   -m motifs_output/homerResults/motif1.motif \
#   -size 200 \
#   > peaks_with_motif1.txt

# # 12. Histogram motywu wokół centrum peaków ±500 bp
# annotatePeaks.pl foxa1_peaks.bed hg38 \
#   -m motifs_output/homerResults/motif1.motif \
#   -size 1000 \
#   -hist 10 \
#   > motif1_histogram.txt

#Motywy sekwencyjne nadreprezentowane w peakach mogą oznaczać:
# miejsce wiązania badanego czynnika transkrypcyjnego
# miejsce wiązania kofaktora
# współdziałanie kilku czynników transkrypcyjnych
# element regulatorowy, np. enhancer/promotor
# potwierdzenie specyficzności eksperymentu

#Składanie genomu to łączenie krótkich odczytów w dłuższe sekwencje:
#odczyty → contigi → scaffoldy → chromosomy
#Contig to ciągła sekwencja złożona z nakładających się odczytów. Scaffold to zbiór contigów uporządkowanych i zorientowanych na podstawie dodatkowych informacji, np. mate-pair, Hi-C, map optycznych.

#Overlap: szukanie nakładania się odczytów
# Layout: budowa grafu połączeń między odczytami
# Consensus: wyznaczenie sekwencji konsensusowej


#LICZENIE FLAG W SAM
with open("reads.sam") as f:
    for line in f:
        if line.startswith("@"):
            continue

        fields = line.strip().split("\t")
        flag = int(fields[1])

        if flag & 4:
            unmapped += 1
#WARUNKI: 
# niezmapowany
flag & 4

# poprawnie sparowany
flag & 2

# nić minus
flag & 16

# primary, czyli NIE secondary i NIE supplementary
not (flag & 256) and not (flag & 2048)

# supplementary
flag & 2048

# niepoprawnie sparowany supplementary
(flag & 2048) and not (flag & 2)
