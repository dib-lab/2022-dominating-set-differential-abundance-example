rule download_host:
    output:
        bison='inputs/host/bison.fna.gz',
        pig='inputs/host/pig.fna.gz',
        cow='inputs/host/cow.fna.gz',
        gp='inputs/host/gp.fna.gz',
        rabbit='inputs/host/rabbit.fna.gz',
        rat='inputs/host/rat.fna.gz',
        sheep='inputs/host/sheep.fna.gz',
        mouse='inputs/host/mouse.fna.gz',
        dog='inputs/host/dog.fna.gz',
        wolf='inputs/host/wolf.fna.gz'
    shell:'''
    wget -O {output.rat} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/895/GCF_000001895.5_Rnor_6.0/
GCF_000001895.5_Rnor_6.0_genomic.fna.gz
    wget -O {output.pig} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11
.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
    wget -O {output.sheep} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/742/125/GCF_002742125.1_Oar_ram
bouillet_v1.0/GCF_002742125.1_Oar_rambouillet_v1.0_genomic.fna.gz
    wget -O {output.rabbit} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/625/GCF_000003625.3_OryCun
2.0/GCF_000003625.3_OryCun2.0_genomic.fna.gz
    wget -O {output.mouse} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38
.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
    wget -O {output.gp} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/735/GCF_000151735.1_Cavpor3.0/
GCF_000151735.1_Cavpor3.0_genomic.fna.gz
    wget -O {output.bison} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/754/665/GCF_000754665.1_Bison_U
MD1.0/GCF_000754665.1_Bison_UMD1.0_genomic.fna.gz
    wget -O {output.cow} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.
2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
    wget -O {output.dog} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.3_CanFam3.1
/GCF_000002285.3_CanFam3.1_genomic.fna.gz
    wget -O {output.wolf} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/922/845/GCA_007922845.1_UniMelb_
Wolf_Refassem_1/GCA_007922845.1_UniMelb_Wolf_Refassem_1_genomic.fna.gz
    '''
