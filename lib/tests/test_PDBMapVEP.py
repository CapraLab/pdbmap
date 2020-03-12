import pytest
import gzip
import logging
from lib import PDBMapGlobals
from lib import PDBMapVEP

logging.basicConfig(level=logging.DEBUG)

def test_PDBMapVEP():
    # First create a little VCF file
    small_vcf_filename = '/tmp/clinvar_excerpt.vcf'
    small_vcf = open(small_vcf_filename,'wb')
    info_linecount = 0
    with gzip.open('/dors/capra_lab/data/clinvar/2019-10-31/clinvar_20191021.vcf.gz','r') as f:
        for line in f:
             small_vcf.write(line)
             if not line.startswith(b'#'):
                 info_linecount += 1
                 if info_linecount > 10:
                     break
    small_vcf.close()

    vep_interface = PDBMapVEP()
    assert vep_interface, "Fundamental troubles constructing PDBMapVEP"

    for x in vep_interface.launch_vep(small_vcf_filename,input_format='vcf'):
        print("You got a line back from vep! -> %s",x)

  

test_PDBMapVEP()
            
