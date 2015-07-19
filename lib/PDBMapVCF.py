#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapVCF.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-12-10
# Description    : Extends PyVCF to customize how the parser
#                : handles 1000 Genomes and VEP fields.
#=============================================================================#

import vcf

class PDBMapVCF(vcf.Reader):

  def __init__(self,*args,**kwargs):
    # Initialize the standard Reader
    parser = vcf.Reader(*args,**kwargs)
    
    # Remove variants on haplotype and non-standard chromosomes
    for record in parser:
      print record.CHROM
    parser = [record for record in parser if record.CHROM[3:] not in range(1,23)+['X','Y','MT']]

    # Determine Info headers
    info_headers = parser.infos.keys()
    # Determine Consequence headers
    csq_headers  = parser.infos['CSQ'].desc.split(': ')[-1].split('|')

    # Clean up the underlying data structure
    self.nsnps = 0
    for record in parser:
      self.nsnps += 1

      # Ensure that an end is specified, default to start+1
      if "END" not in record.INFO:
        record.INFO["END"] = int(record.POS) + 1

      # Ensure that each record contains all fields
      for header in info_headers:
        if header not in record.INFO:
          record.INFO[header] = None

      # Process fields which may or may not be included in VEP output
      if 'SNPSOURCE' in record.INFO and record.INFO['SNPSOURCE']:
        record.INFO['SNPSOURCE'] = ';'.join(record.INFO['SNPSOURCE'])
      else: record.INFO['SNPSOURCE'] = None

      # Make any necessary population allele frequency corrections
      freqs = ['AMR_AF','AFR_AF','EUR_AF','EAS_AF','SAS_AF']
      if 'ASN_AF' in record.INFO:
        record.INFO['EAS_AF'] = record.INFO['ASN_AF']
        record.INFO['SAS_AF'] = record.INFO['ASN_AF']
      for freq in freqs:
        if freq not in record.INFO: record.INFO[freq] = None
        recordINFO[freq] = 0.0 if not record.INFO[freq] else record.INFO[freq]

      # Ensure the ancestral allele is encoded properly
      if 'AA' in record.INFO and record.INFO['AA']:
        record.INFO['AA'] = record.INFO['AA'].upper()
      else: 
        record.INFO['AA'] = None

      # Enforce biallelic assumption
      # Record only the first alternate allele count
      if 'AC' in record.INFO:
        record.INFO['AC'] = record.INFO['AC'][0]
      else: 
        record.INFO['AC'] = None

      # Ensure necessary fields are present in record.INFO
      if 'AN'      not in record.INFO: record.INFO['AN']      = None
      if 'AF'      not in record.INFO: record.INFO['AF']      = None
      if 'VT'      not in record.INFO: record.INFO['VT']      = None
      if 'SVTYPE'  not in record.INFO: record.INFO['SVTYPE']  = None
      if 'SVLEN'   not in record.INFO: record.INFO['SVLEN']   = None
      if 'AVGPOST' not in record.INFO: record.INFO['AVGPOST'] = None
      if 'RSQ'     not in record.INFO: record.INFO['RSQ']     = None
      if 'ERATE'   not in record.INFO: record.INFO['ERATE']   = None
      if 'THETA'   not in record.INFO: record.INFO['THETA']   = None
      if 'LDAF'    not in record.INFO: record.INFO['LDAF']    = None

      # Allele frequency is sometimes reecorded as a tuple or list
      # Extract the first (only) element and cast to float
      if type(record.INFO['AF']) in [type((None,)),type([])]:
        record.INFO['AF'] = float(record.INFO['AF'][0])

      # Add attribute fields to INFO
      record.INFO["ID"]     = record.ID
      record.INFO["CHROM"]  = record.CHROM
      record.INFO["START"]  = int(record.POS)
      record.INFO["REF"]    = record.REF
      record.INFO["ALT"]    = record.ALT[0] # Enforce biallelic assumption
      record.INFO["QUAL"]   = record.QUAL
      record.INFO["FILTER"] = str(record.FILTER)

      # Determine and record the derived allele
      if record.INFO["REF"] == record.INFO["AA"]:
        record.INFO["DA"] = record.INFO["ALT"]
      else:
        record.INFO["DA"] = record.INFO["REF"]

      # If VEP consequences was not included, ensure that some fields are defined
      if CSQ not in record.__dict__:
        record.INFO["EXISTING"] = None
        record.INFO["GENE"]     = None
        record.INFO["HGNC"]     = None
        record.INFO["ID"] = "%(CHROM)s:%(START)d-%(END)d"%(record.INFO)
      # Parse VEP consequences if included
      else:
        record.CSQ = self._parse_csq(csq_headers,record.INFO['CSQ'])

        # Move some of the information to INFO
        record.INFO["EXISTING"] = record.CSQ[0]['Existing_variation']
        record.INFO["GENE"]     = record.CSQ[0]['Gene']

        if record.CSQ[0]['SYMBOL_SOURCE'] == 'HGNC':
          record.INFO["HGNC"] = record.CSQ[0]["SYMBOL"]
        else: record.INFO["HGNC"] = None

        # Generate unique idenfitiers for new variants
        if not record.INFO["ID"]:
          if record.INFO["EXISTING"]:
            record.INFO["ID"] = record.INFO["EXISTING"]
          else:
            record.INFO["ID"] = "%(CHROM)s:%(START)d-%(END)d"%(record.INFO)

        # Add some variant info to the consequences
        for csq in record.CSQ:
          csq["ID"]    = record.INFO["ID"]
          csq["CHROM"] = record.INFO["CHROM"]
          csq["START"] = record.INFO["START"]
          csq["END"]   = record.INFO["END"]

  def _parse_csq(self,csq_header,csq_values):
    """ Creates a dictionary from VEP CSQ desc and each row of values """
    # grab column headers
    res = []
    for row in csq_values:
      csq = {}
      ## Edits to make while parsing fields (inner loop)
      for i,field in enumerate(row.split('|')):
        csq[csq_header[i]] = field
        # Check the consequence and reformat
        if csq_header[i] == "Consequence":
          cons = field.split('&')
          ## We are now allowing Synonymous SNPs to be mapped ##
          # if not any([con in self._parse_csq.nonsyn for con in cons]):
          #   csq = None
          #   break # Ignore this row. Invalid consequence.
          csq['Consequence'] = ';'.join(cons)
        # Set any empty strings to None and continue
        elif csq[csq_header[i]] == '':
          csq[csq_header[i]] = None
          continue
        # Reformat the reference and alternate amino acids
        elif csq_header[i] == "Amino_acids":
          if not field: # No Value
            ref,alt = None,None
          elif len(field) < 2: # Only one value
            ref,alt = field,field
          else: # Both specified
            ref,alt = field.split('/')[0:2]
          csq["Ref_AminoAcid"] = ref
          csq["Alt_AminoAcid"] = alt
        # Reformat the reference and alternate codons
        elif csq_header[i] == "Codons":
          if not field: # No value
            ref,alt = None,None
          elif len(field) < 2: # Only one value
            ref,alt = field,field
          else: # Both specified
            ref,alt = field.split('/')[0:2]
          csq["Ref_Codon"] = ref
          csq["Alt_Codon"] = alt
        # Reformat the existing variation
        elif csq_header[i] == "Existing_variation":
          csq["Existing_variation"] = csq["Existing_variation"].replace('&',',')
        # Reformat the domains
        elif csq_header[i] == "DOMAINS":
          csq['DOMAINS'] = csq['DOMAINS'].replace('&',',')
        # Convert canonical flag to boolean
        elif csq_header[i] == "CANONICAL":
          csq["CANONICAL"] = 1 if csq["CANONICAL"] == "YES" else 0
        ##Transform valid 0-indexed positions to 1-indexed positions
        #UPDATE: Despite initial appearances, positions seem to be correct
        elif csq_header[i] == "Protein_position":
          if csq['Protein_position'] and '?' not in csq['Protein_position']:
            csq["Protein_position"] =  int(csq['Protein_position'].split('-')[0])# + 1
          else: csq["Protein_position"] = None
        elif csq_header[i] == "cDNA_position":
          if csq['cDNA_position'] and '?' not in csq['cDNA_position']:
            csq["cDNA_position"] =  int(csq['cDNA_position'].split('-')[0])# + 1
          else: csq["cDNA_position"] = None
        elif csq_header[i] == "CDS_position":
          if csq['CDS_position'] and '?' not in csq['CDS_position']:
            csq["CDS_position"] =  int(csq['CDS_position'].split('-')[0])# + 1
          else: csq["CDS_position"] = None

      if not csq: continue # Invalid conseqeunce. Skip the entry.

      ## Edits to make after parsing fields completely (outer loop)
      # Transcript is not canonical if no value was given
      if "CANONICAL" not in csq: csq["CANONICAL"] = 0
      # Define ref/alt amino acids if not specified
      if "Amino_acids" not in csq or not csq["Amino_acids"]:
        csq["Ref_AminoAcid"] = None
        csq["Alt_AminoAcid"] = None
      # Define ref/alt codons if not specified
      if "Codons" not in csq or not csq["Codons"]:
        csq["Ref_Codon"] = None
        csq["Alt_Codon"] = None
      # For any field not specified, add with value None
      for header in csq_header:
        if header not in csq:
          csq[header] = None
        elif csq[header] == '':
          csq[header] = None
      res.append(csq)
    return res
    # The Variant Effect Predictor will return all consequences for a
    # variant if any one of those consequences passes the filter. This
    # list is used to filter the erroneous consequences of otherwise
    # relevant variants.
  _parse_csq.nonsyn = ['initiator_codon_variant','inframe_insertion',
                      'inframe_deletion','missense_variant','stop_gained',
                      'stop_lost','frameshift_variant','splice_region_variant',
                      'transcript_ablation']
