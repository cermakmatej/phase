# phase
Bakalářská práce na téma fázování ze sekvenačních dat

### Instalace závislostí
Program využívá knihovny numpy a tqdm, závislosti se nacházejí v souboru `requirements.txt`.

Pro správné fungování programu je potřeba mít nainstalovaný nástroj **bcftools** (verze 1.22) (https://samtools.github.io/bcftools/).

Kromě toho je nutné do `bcftools` přidat plugin **`+process`**.  
Zdrojový soubor pluginu `process.c` se nachází v adresáři `utils`.

### Instalace pluginu `+process`

Nejjednodušší postup je následující:

1. Přesuňte soubor `process.c` do adresáře s pluginy vaší lokální instalace `bcftools`.
2. V tomto adresáři spusťte příkaz:

   ```bash
   make

### Program
Program se nachází v adresáři `program`. Lze spustit v návaznosti na zmíněný plugin pomocí příkazu:
   ```bash
bcftools +process -a alignment.bam -f reference.fasta -s variants.vcf \
| python program/phase.py --variant_file variants.vcf > phased.txt
```

### Simulace
Simulace kromě **bcftools** využívá nástroje **art** (verze 2.5.8) (https://www.niehs.nih.gov/research/resources/software/biostatistics/art), **bwa** (verze 0.7.17-r1188) (https://bio-bwa.sourceforge.net/) a **samtools** (verze verze 1.15.1) (https://www.htslib.org/). Nástroj **art** je přiložen v adresáři `utils` .Závislosti jsou uvedeny v souboru `sim_requirements.txt`.

Celá simulace lze spustit pomocí skriptu `simulation_complete.sh`. První krok simulace, náhodný výběr variant z katalogu **GnomAD** (verze v2.1.1) vyžaduje velký soubor `gnomad.genomes.r2.1.1.sites.22.vcf.bgz` (https://gnomad.broadinstitute.org/downloads#v2). Proto je přiložen soubor `simulated_indiv_chr22.vcf` s již vybranými variantami.  

### Výsledky
Výsledky se nacházejí v adresáři `results`. V adresáři `real_data` se nacházejí tři soubory: `variants22.vcf` jsou vyfiltrované a určené varianty, `hets22.vcf` jsou pouze heterozygotní varianty a `real_output22.txt` je výstup algoritmu pro tyto varianty. Práce pracuje se sekvenačními daty v souboru `RMNISTHS\_30xdownsample.bam`, který lze stáhnout z repozitáře (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/).

V adresáři `simulated_data` se nachází soubor `simulated_hets22.vcf`, který obsahuje heterozygotní simulované varianty. V adresáři `true_phased` se nacházejí VCF soubory s uměle zfázovanými variantami pro experimenty s rostoucím zastoupením vložených repetic. Pro každý z těchto experimentů najdeme výstup fázování v adresáři `outputs` a vyhodnocení fázování v adresáři `evaluations`.
