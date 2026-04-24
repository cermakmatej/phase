# phase
Bakalářská práce na téma fázování ze sekvenačních dat

### Instalace závislostí

Program využívá knihovny numpy a tqdm, závislosti se nacházejí v souboru requirements.txt.

Pro správné fungování programu je potřeba mít nainstalovaný nástroj **bcftools** (https://samtools.github.io/bcftools/).

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
| python program/phase.py --variant_file variants.fasta > phased.txt
