#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapVep.py
# Author         : Chris Moth, refactored from R. Michael Sivley's PDBMapVep
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-18 / 2019-11-05
#=============================================================================#

""" Interface to the ENSEMBL Variant Effect Predictor, which returns
    the transcript consequences of genomic level (chr and position) changes"""

# See main check for cmd line parsing
import sys
import os
import gzip
import logging
import string
import vcf
from typing import Iterator,Dict
import subprocess as sp
from lib.PDBMapGlobals import PDBMapGlobals
from lib.PDBMapSQLdb import PDBMapSQLdb
from lib import bed # PyVCF emulator for BED files
LOGGER = logging.getLogger(__name__)
SQL_EXTENDED_LOCK_WAIT = "SET innodb_lock_wait_timeout=200;\n"

class PDBMapVEP():
    """The ENSEMBL Variant Effect Predictor is a command line utility with
    many powerful options  https://www.ensembl.org/vep

    class PDBMapVEP manages the VEP command line parameters for PDBMap applications
    which typically employ VEP to convert genomic .vcf data ("clinvar", "exac", etc)
    to ENST-referenced transcript-level genetic impacts.  These impacts are frequently
    loaded into SQL databases

    VEP output is filtered to remove VEP-returned data that is spurious for the PDBMap
    application suite
    """

    def __init__(self,config_dict: Dict[str,str]  = None):
        """Construct a PDBMapVEP() interface by providing a dictionary containing:
              vep: The full path to the ENSEMBL vep binary program

            and also one of these 3:

                vep_cache_dir: The ENSEMBL-recommended datasource for vep

            or: if "$ENSEMBL_REGISTRY" present in the environment, vep is run with --registry flag

            or: dbhost/dbuser/dbpass access parameters to the ENSEMBL SQL database

         """

        self._config_dict = config_dict
        if not self._config_dict:
           LOGGER.info("Initializing PDBMapVEP with PDBMapGlobals.config dictionary")
           self._config_dict = PDBMapGlobals.config
        if not self._config_dict or 'vep' not in self._config_dict:
            raise Exception("The vep executable is missing from the config dictionary, or invalid")
        self.vep_executable = self._config_dict['vep']
        self.vep_cache_dir = None

        self.ensembl_registry_filename = None

        self.dbhost = None
        self.dbuser = None
        # self.dbname = None
        self.dbpass = None
 
        if not self.vep_executable or not os.path.exists(self.vep_executable):
            msg = "VEP executable program path is invalid: %s"%self.vep_executable
            LOGGER.critical(msg)
            sys.exit(msg)

        # import pdb; pdb.set_trace()
        # Vep preferred access method #1
        if 'vep_cache_dir' in self._config_dict:
            self.vep_cache_dir = self._config_dict['vep_cache_dir']
            if not os.path.exists(self.vep_cache_dir):
                msg = "VEP cache was not found at provided directory: %s"%self.vep_cache_dir
                LOGGER.critical(msg)
                sys.exit(msg)
            LOGGER.info("PDBMapVEP() constructed with\nvep_executable=%s\nvep_cache_dir= %s", self.vep_executable, self.vep_cache_dir)
        else: # MEthod #2 Access the Ensembl SQL files like the rest of the ENSEMBL PERL API
            self.ensembl_registry_filename = os.getenv('ENSEMBL_REGISTRY')
            if self.ensembl_registry_filename:
                LOGGER.info("vep communications settings to be taken from: %s",self.ensembl_registry_filename)
            else: # #3 Access the Ensembl SQL via specific mysql connection parameters
                for db_param in ['dbhost','dbuser','dbpass']:
                    if db_param not in self._config_dict:
                        msg = "DB parameter %s missing from dictionary, and neither vep_cache_dir nor $ENSEMBL_REGISTRY are defined. Halting."%db_param
                        LOGGER.critical(msg)
                        sys.exit(msg)
                LOGGER.warning("No VEP (--cache) --dir was provided to PDBMapVEP()\n" +
                           "per ENSEMBL recommendations. VEP run will use SQL queries.")
                sys.exit("Terminating for now until VEP cahcing is better understood")
                



    def launch_vep(self,
                input_filename: str,
                input_format: str = 'default',
                generator_echo_filename: str = None) -> Iterator[str]:

        """1) Build a VEP ready command line,
        2) launch vep
        3) yield strings for each vcf output line that passes filter for
           being consequential to protein sequence (have CSQ)"""

        if not input_filename or not input_filename.strip():
            input_filename = 'stdin'
            LOGGER.critical(\
    "PDBMapVEP.run_VEP() requires an input filename.  stdin is not an option at present")

        vep_cmd = [self.vep_executable, '-i', input_filename]
        # Specify the input type.  By default, VEP will auto-detect the input format
        if input_format and input_format != 'default':
            vep_cmd.extend(['--format', input_format])
        # Use the local VEP cache if it exists
        if self.vep_cache_dir:
            vep_cmd.extend([
                '--cache',   # VEP will emply the local cache files
                '--dir',     # Tell VEP to use the cache that follows
                # VEP will look in this directory for downloaded cache elements:w
                self.vep_cache_dir,
                '--offline', # VEP will not reach out to SQL database(s)
                '--pubmed',  # output publications, a feature only available with a downloaded cache
            ])
        elif self.ensembl_registry_filename:
            vep_cmd.extend(['--registry',self.ensembl_registry_filename])
        else: # VEP will use an (ideally local) ENSEMBL database
            vep_cmd.extend([
                '--database',
                '--host',self.vep_dbhost,
                '--pass',self.vep_dbpass])

        # Disable the progress bars
        vep_cmd.extend(['--no_progress'])
        # Increase buffer size 100x to improve runtime (default 5000)
        vep_cmd.extend(['--buffer_size', '100000']) # ~5GB ?RAM? per 100,000

        # Annotate with functional info/prediction.  From vep documentation:
        # Predict whether an AA substitution affects protein function based on
        # sequence homology and the physical properties of amino acids.
        # 's' outputs a score only ('-p' would output a 'prediction term')
        vep_cmd.extend(['--sift', 's'])

        # Human only predicts possible impact on AA substution on
        # structure and function using straightforward physical and comparative
        # considerations.  's' option outputs score only
        vep_cmd.extend(['--polyphen', 's'])

        # Annotate with variant, gene, protein, and domain identifiers
        vep_cmd.extend(['--check_existing', '--symbol', '--protein', '--uniprot', '--domains'])

        # Annotate transcript with canonical boolean 
        vep_cmd.extend(['--canonical'])
        # Annotate transcript with canonical boolean and biotype (a gene or transcript classification)
        vep_cmd.extend(['--biotype'])

        # Retain only coding-region variants
        vep_cmd.extend(['--coding_only'])
        #  output format
        vep_cmd.extend(['--vcf'])

        #  SOMETHING IS WRONG WITH db_version
        #  DO NOT CHECK THIS IN 
        LOGGER.critical("*** REMOVE --db_version SOON ***")
        vep_cmd.extend(['--db_version','96'])


        # Send to stdout and don't generate a summary stats file
        vep_cmd.extend(['--no_stats'])
        vep_cmd.extend(['-o', 'stdout'])
        LOGGER.info("Invoking VEP with commands: %s", ' '.join(vep_cmd))

        # Call VEP and capture stdout in realtime
        VEP_process = sp.Popen(vep_cmd, stdout=sp.PIPE, bufsize=1)
        echo_f = None
        if generator_echo_filename:
            echo_f = gzip.open(generator_echo_filename, 'wb') # Open cache for writing

        def CSQ_in_INFO(vep_output_line):
            vcf_split_tabs = vep_output_line.split(b'\t')
            if len(vcf_split_tabs) < 8:
                LOGGER.info("Mal-formed (too few tabs) VEP output line",vep_output_line)
                return False
            # The INFO column is the final (right-most) of the vep output columns 
            return vcf_split_tabs[7].find(b';CSQ=') > 0

        vep_lines_yielded = 0
        for vep_output_line in iter(VEP_process.stdout.readline, b''):
            # print('VEP OUTPUT LINE IS %s'%vep_output_line.decode("utf-8"))
            # Filter variants without consequence (CSQ) annotations
            if vep_output_line.startswith(b'#') or CSQ_in_INFO(vep_output_line):
                # Echo to file, if caller requested, before yielding
                if echo_f:
                    echo_f.write(vep_output_line)
                vep_lines_yielded += 1
                yield vep_output_line.decode('utf-8')

        LOGGER.info("%d lines (strings) yielded from VEP"%vep_lines_yielded)

        if echo_f:
            echo_f.close() # Close the cache
            del echo_f
        LOGGER.critical("*** REMOVE --db_version SOON ***")

        """
        except KeyboardInterrupt:
            msg = "ERROR (PDBMapVep) Keyboard interrupt. Canceling VEP...\n"
            sys.stderr.write(msg)
            raise
        except:
            msg = "ERROR (PDBMapVep) Unknown VEP error. Examine stacktrace...\n"
            sys.stderr.write(msg)
        raise

        finally
        """

        VEP_process.kill() # kill any running subprocess, Does not seem to throw exceptions...



    def _parse_csq(self,csq_header,vcf_record):
        """ Creates a dictionary from VEP CSQ desc and each row of values """
        # The Variant Effect Predictor will return all consequences for a
        # variant if any one of those consequences passes the filter. This
        # list is used to filter the erroneous consequences of otherwise
        # relevant variants.
        _parse_csq_nonsyn = ['initiator_codon_variant','inframe_insertion',
                            'inframe_deletion','missense_variant','stop_gained',
                            'stop_lost','frameshift_variant','splice_region_variant',
                            'transcript_ablation']

        # grab column headers
        res = []
        if 'CSQ' in vcf_record.INFO and vcf_record.INFO['CSQ']:
            csq_values = vcf_record.INFO['CSQ']
        elif 'vep' in vcf_record.INFO and vcf_record.INFO['vep']:
            csq_values = vcf_record.INFO['vep']
        else:
            LOGGER.warning("No INFO['CSQ'] or INFO['VEP'] in vcf_record %s"%str(vcf_Record))
            csq_values = []
        for row in csq_values:
            csq = {}
            ## Edits to make while parsing fields (inner loop)
            for i,field in enumerate(row.split('|')):
                csq[csq_header[i]] = field
                # Check the consequence and reformat
                if csq_header[i] == "Consequence":
                    cons = field.split('&')
                    ## We are now allowing Synonymous SNPs to be mapped ##
                    # if not any([con in _parse_csq_nonsyn for con in cons]):
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
                    elif len(field.split('/')) < 2: # Only one amino acid recorded
                        ref,alt = field,field
                    else: # Both specified
                        ref,alt = field.split('/')[0:2]
                    csq["Ref_AminoAcid"] = ref
                    csq["Alt_AminoAcid"] = alt
                # Reformat the reference and alternate codons
                elif csq_header[i] == "Codons":
                    if not field: # No value
                        ref,alt = None,None
                    elif len(field.split('/')) < 2: # Only one codon recorded
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
    def _complete_vep_record_parse_and_clean(self,record,info_headers,csq_headers):
        """Perform numerous changes to VCF records initially parsed from VEP output, so that
           the record is ready for addition to SQL"""

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
        else: 
            record.INFO['SNPSOURCE'] = None

        # Make any necessary population allele frequency corrections
        # Enforce biallelic assumption
        population_allele_freqs = ['AMR_AF','AFR_AF','EUR_AF','EAS_AF','SAS_AF','ASN_AF']
        for allele_freq in population_allele_freqs:
            if (allele_freq not in record.INFO) or (not record.INFO[allele_freq]): 
                record.INFO[allele_freq] = 0.0 # Ensure all fields present
            
            # If the type of this .info is a tuple of None or list, grab the FIRST 
            # component of the tuple or list as the allele frequency tuple or list
            if type(record.INFO[allele_freq]) in [type((None,)),type([])]:
                record.INFO[allele_freq] = float(record.INFO[allele_freq][0])

        # Ensure the ancestral allele is encoded properly
        if 'AA' in record.INFO and record.INFO['AA']:
            # Format: AA|REF|ALT|INDELTYPE
            record.INFO['AA'] = record.INFO['AA'].split('|')[0].upper()
        else: 
            record.INFO['AA'] = None

        # Extract the format and flatten genotypes (comma delimit) at this site
        record.INFO['FORMAT'] = record.FORMAT
        record.INFO['GT'] = ','.join([str(sample['GT']) for sample in record.samples])

        # Enforce biallelic assumption
        # Record only the first alternate allele count
        record.INFO['AC'] = record.INFO['AC'][0] if 'AC' in record.INFO and record.INFO['AC'] else None

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

        # Reformat the variant-type value
        if type(record.INFO["VT"]) == type(tuple()):
            record.INFO["VT"] = ','.join(vt for vt in list(record.INFO["VT"]))

        # Allele frequency is sometimes reecorded as a tuple or list
        # Enforce biallelic assumption
        # Extract the first (only) element and cast to float
        if type(record.INFO['AF']) in [type((None,)),type([])]:
            record.INFO['AF'] = float(record.INFO['AF'][0])

        # Add attribute fields to INFO
        record.INFO["ID"]     = record.ID
        record.INFO["CHROM"]  = record.CHROM
        record.INFO["POS"]  = int(record.POS)
        record.INFO["REF"]    = record.REF
        record.INFO["ALT"]    = record.ALT[0] # Enforce biallelic assumption
        record.INFO["QUAL"]   = record.QUAL
        record.INFO["FILTER"] = str(record.FILTER)

        # Infer the derived allele from the ancestral
        # If ancestral not specified, assume minor allele is derived
        if record.INFO["ALT"] == record.INFO["AA"]:
            record.INFO["DA"] = record.INFO["REF"]
        else:
            record.INFO["DA"] = record.INFO["ALT"]      
        
        # If CSQ missing for this record or all records, dummy code
        if not csq_headers or not (
              ('CSQ' in record.INFO and record.INFO["CSQ"]) or 
              ('vep' in record.INFO and record.INFO["vep"])):
            record.CSQ = [{"ID":   record.INFO["ID"],
                            "CHROM": record.INFO["CHROM"],
                            "POS": record.INFO["POS"],
                            "END":   record.INFO["END"],
                            "Feature":          None,
                            "ENSP":             None,
                            "CANONICAL":        None,
                            "Allele":           None,
                            "Consequence":      None,
                            "cDNA_position":    None,
                            "CDS_position":     None,
                            "Protein_position": None,
                            "Ref_AminoAcid":    None,
                            "Alt_AminoAcid":    None,
                            "Ref_Codon":        None,
                            "Alt_Codon":        None,
                            "PolyPhen":         None,
                            "SIFT":             None,
                            "BIOTYPE":          None,
                            "DOMAINS":          None,
                            "SWISSPROT":        None}]
            record.INFO['EXISTING']  = None
            record.INFO["GENE"]      = None
            record.INFO["HGNC"]      = None
            record.INFO["SWISSPROT"] = None
        else:
            # Extract the consequences
            record.CSQ = self._parse_csq(csq_headers,record)
            # Add some "consequence" info to the variant info
            record.INFO["EXISTING"] = record.CSQ[0]['Existing_variation']
            record.INFO["GENE"] = record.CSQ[0]['Gene']
            if record.CSQ[0]['SYMBOL_SOURCE'] == 'HGNC':
                record.INFO["HGNC"] = record.CSQ[0]["SYMBOL"]
            else: record.INFO["HGNC"] = None
            # Correct missing variant IDs
            if not record.INFO["ID"]:
                if record.INFO["EXISTING"]:
                  record.INFO["ID"] = record.INFO["EXISTING"]
                else:
                  record.INFO["ID"] = "%(CHROM)s:%(POS)d-%(END)d"%(record.INFO)
            # Add some variant info to the consequences
            for csq in record.CSQ:
                csq["ID"] = record.INFO["ID"]
                csq["CHROM"] = record.INFO["CHROM"]
                csq["POS"] = record.INFO["POS"]
                csq["END"] = record.INFO["END"]
        return record

    @staticmethod
    def _my_vcf_file_Reader(vcf_filename):
        vcf_file_compressed = None # <- This works great with vcf.Reader() calls involving .gz (compressed automatically) and all else (not compressed)
        if vcf_filename.endswith(".bgz"): # Why does Gnomad use .bgz for .gz files?  So irritating!  But, we override and all works out
            vcf_file_compressed = True
        return vcf.Reader(filename=vcf_filename,
            prepend_chr=True,
            compressed=vcf_file_compressed,
            encoding='UTF-8')
        
    def vcf_file_to_sql_table(self,table_name: str,vcf_filename: str,buffer_size:int =1,local_db_config: Dict[str,str] = None):
        """ Creates a supplementary SQL table as echo of  VCF datafile
            The first 7 VCF tab-delimited data columns are extracted.
            The 8th column, INFO, is parsed apart and SQL columns created
            for each element

            For more details, be sure to google the VCF4.x specification .pdf
            and PyVCF library documentation """

        if not local_db_config:
            local_db_config = {}
            for dbparam in ['dbhost','dbname','dbuser','dbpass']:
                assert dbparam in self._config_dict,"%s missing from configuration dictionary"%dbparam
                local_db_config[dbparam] = self._config_dict[dbparam]

        for dbparam in ['dbhost','dbname','dbuser','dbpass']:
            assert dbparam in local_db_config,"%s missing from local_db_config dictionary"%dbparam

        # Open a PyVCF parser  and initialize vcf_reader.infos from the ## header
        # The parser will later iterate over vcf records in the input vcf_filename
        #
        #  rameter prepend_chr=True, prepends the string 'chr' to all the chromosome
        #   values found (as in chr1/chrX/chrEtc)
        vcf_reader = self._my_vcf_file_Reader(vcf_filename)

        # The first 7 mandatory vcf columns are consistent from vcf file to vcf file
        # var_header  = ["chr","start","name","ref","alt","qual","filter"]
        # 2019-Nov.  Chris Moth strongly prefers to maintain source data field names
        first_7_columns  = ["chrom","pos","id","ref","alt","qual","filter"]
        first_7_comments = ["Chromosome","position","identifier",
                            "reference base(s)","alternate base(s)","quality",
                            "filter status PASS/MISSING/etc"]
                            
        # SQL types for these first 7 standard mandatory fields
        first_7_types  = ["VARCHAR(100)","BIGINT","VARCHAR(100)","VARCHAR(100)"]
        first_7_types += ["VARCHAR(100)","DOUBLE","VARCHAR(100)"]

        primary_key_components = ['chrom','pos','ref','vcf_record_number']

        # Replace INFO keys that overlap with the first_7 column headers
        for info_key in list(vcf_reader.infos.keys()):
            if info_key.lower() in first_7_columns:
              # If the INFO has a key same as first 7
              # then rename that as key2 to avoid confusion
              vcf_reader.infos["%s2"%info_key] = vcf_reader.infos[info_key]
              del vcf_reader.infos[info_key] # Remove the duplicating info_key from the dictionary

        # Extract info fields
        info_header = list(vcf_reader.infos.keys())
        # replace single quotes in the comment strings with two single adjacent quotes, for SQL
        info_comments = [info.desc.replace("'","''") for info in list(vcf_reader.infos.values())]
        # Extract and convert info types from VCF INFO ## Meta information to SQL types
        type_conv = {"Integer":"BIGINT","Float":"DOUBLE","Flag":"TINYINT(1)","String":"TEXT"}
        info_types  = [type_conv[info.type] for info in list(vcf_reader.infos.values())]

        # replace all punctuation with underscores in the info keys that could foul up SQL
        punctuation_to_underscore = str.maketrans(string.punctuation,'_'*len(string.punctuation))
        info_header = [f.translate(punctuation_to_underscore)  for f in info_header]

        header = first_7_columns + info_header #+ csq_header
        comments = first_7_comments + info_comments #+ csq_header
        # Use the first row to infer data types for the INFO fields
        record = next(vcf_reader)
        types    = first_7_types + info_types #+ csq_types

        # Set default values for each type
        # In our Mariadb 5.5, TEXT has no DEFAULT
        sql_defaults = {"BIGINT":0,"DOUBLE":0.0,"TINYINT(1)":0,"VARCHAR(100)":"''"}
        sql_notnull_types  = set() # Not sure why Mike had notnull on these: was-> set["TINYINT","TEXT","VARCHAR(100)","DOUBLE"])
        # don't worry about type conversion for now, Nones are causing an issue
        formatter = {"BIGINT":"%s","DOUBLE":"%s","TEXT":"%s","TINYINT(1)":"%s","VARCHAR(100)":"%s"}
        # Generate a create statement for this table
        # First, generate the column definitions
        # Each line defines:
        # a table column name, backquoted to avoid conflicts with SQL reserved words
        # the sql type of the column
        # A default value
        # A non null specifier
        table_columns_def = (
            ["vcf_record_number BIGINT NOT NULL AUTO_INCREMENT"] + 
            ["`%s` %s %s %s COMMENT '%s'"%(header[i].lower(), \
               types[i], \
               ("DEFAULT %s"%sql_defaults[types[i]]) if types[i] in sql_defaults else "", \
               "NOT NULL" if (header[i].lower in primary_key_components) or (types[i] in sql_notnull_types) else "", \
               comments[i]) \
                        for i in range(len(header))]
            )

        # Include as many non-VARCHAR(65000)/TEXT columns in primary key as allowed (16)
        # Additional INFO (like END) may justify duplicates within the standard VCF fields
        dbname_dot_table_name = local_db_config['dbname'] + "." + table_name

        # Dropping the table is a bad idea if 25 Slurm tasks attempting to load it
        # with PDBMapSQLdb() as db:
        #     db.execute("DROP TABLE IF EXISTS %s;"%dbname_dot_table_name)
        create_statement = "CREATE TABLE IF NOT EXISTS %s.%s\n(%s, PRIMARY KEY(%s), UNIQUE KEY(vcf_record_number))\n%s\n%s"
        query = create_statement%(local_db_config['dbname'],
            table_name,
            ',\n'.join(table_columns_def), # The columns with their types and defaults
            ',\n'.join(primary_key_components),
            "CHARACTER SET utf8mb4 COLLATE utf8mb4_general_ci",
            "COMMENT 'created from %s'"%vcf_filename)
        LOGGER.info("Creating table %s with 'query'=\n%s"%(table_name,query))

        # Create the table
        with PDBMapSQLdb() as db:
            db.execute(query)


        # Populate the table with contents of VCF file reformatted for clean INSERT
        # Strategy is one INSERT with many many value tuples
        def record2row(record,infos):#,csq_header=None):
            row  = [record.CHROM,record.POS,record.ID]
            row += [record.REF,record.ALT,record.QUAL,record.FILTER]
            # Only retain the most frequent alternate allele (force biallelic)
            row += [record.INFO[f] if f in record.INFO else None for f in list(infos.keys())]
            # Replace any empty lists with None
            row =  [r if type(r)!=list or len(r)<1 else r[0] for r in row]
            row =  [r if r!=[] else None for r in row]
            # Change None to DEFAULT in 
            return row
        # Add the first record to insert rows
        unwritten_vcf_rows   = [record2row(record,vcf_reader.infos)]
        rows_inserted = 0
        query  = SQL_EXTENDED_LOCK_WAIT 
        query += "INSERT INTO %s "%dbname_dot_table_name
        query += "(%s) VALUES "%','.join(['`%s`'%h for h in header])
        query += "(%s)"%','.join([formatter[f] for f in types])

        def write_unwritten_vcf_rows():
            nonlocal rows_inserted, query, unwritten_vcf_rows
            # try:
            # import pdb; pdb.set_trace()

            with PDBMapSQLdb() as db:
                db.set_session_transaction_isolation_read_committed()
                rows_inserted += db.executemany(query,unwritten_vcf_rows)
            # We flushed the rows so now we re-initialize for more rows
            unwritten_vcf_rows  = []
            return rows_inserted


        for i,record in enumerate(vcf_reader):
            # Buffered upload - Flush as specified in config
            if not (i+1) % buffer_size: # Flush this one!
                write_unwritten_vcf_rows()
                unwritten_vcf_rows  = []

            # More often, we keep building up our big INSERT statement
            unwritten_vcf_rows.extend([record2row(record,vcf_reader.infos)])

        if unwritten_vcf_rows:
            write_unwritten_vcf_rows()


        # Return the number of rows uploaded
        return rows_inserted

    def vcf_reader_from_file_supplemented_with_vep_outputs(self,vcf_filename,vep_echo_filename=None):
        """ Pipe a VCF file into VEP and yield records suitable for additional to SQL """
        # vep_echo_filename  = '/tmp/VEP_stdout_%s'%os.path.basename(vcf_filename)
        # if vep_echo_filename.split('.')[-1] != 'gz':
        #     vep_echo_filename += '.gz'
        LOGGER.info("Piping %s to VEP."%vcf_filename)
        if vep_echo_filename:
            LOGGER.info("Echoing VEP output to %s"%vep_echo_filename)
        # Launch the vep, and return a vcf reader that parses vep outputs
        vcf_reader = vcf.Reader(self.launch_vep(vcf_filename,'vcf',vep_echo_filename),prepend_chr=True,encoding='UTF-8')
        vcf_reader.CSQorVEP = None
        for key in ['CSQ','vep']:
            if key in vcf_reader.infos:
                vcf_reader.CSQorVEP = key
                break
        if not vcf_reader.CSQorVEP:
            LOGGER.critical("Neither CSQ nor VEP found in the vcf_reader.infos field.  Please check integratity of input vcf file")
            # Mike Sivley note:
            # This can happen if the VCF already includes INFO columns and
            # the first row does not pass the consequence filtering criteria.
            # The actual data rows will have CSQ  / VEP information.
            # Manually add the VEP CSQ field to the VCF parser.
            vcf_reader.infos["CSQ"] = vcf.parser._Info("CSQ",'.',"String","Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|CLIN_SIG|SOMATIC|PHENO|PUBMED")
            vcf_reader.CSQorVEP = 'CSQ'
        return vcf_reader

    


    def vcf_reader_from_file_without_vep(self,vcf_filename):
        vcf_reader = self._my_vcf_file_Reader(vcf_filename)

        vcf_reader.CSQorVEP = None
        for key in ['CSQ','vep']:
            if key in vcf_reader.infos:
                vcf_reader.CSQorVEP = key
                break
        if not vcf_reader.CSQorVEP:
            LOGGER.critical("You are skipping the VEP step - but there is no CSQ in your .vcf file %s INFO definition!!!"%vcf_filename)
            sys.exit(1)
            # Later perhaps:   This is not wired in yet...  Not sure...  If no CSQ present, insert dummy structure
            vcf_reader.CSQorVEP = 'CSQ'
            vcf_reader.infos['CSQ'] = bed.Consequence() # dummy object
        return vcf_reader

    def yield_completed_vcf_records(self,vcf_reader):
        # Determine Info headers
        info_headers = list(vcf_reader.infos.keys())
        # Determine Consequence headers
        csq_headers  = vcf_reader.infos[vcf_reader.CSQorVEP].desc.split(': ')[-1].split('|')

        snpcount = 0
        parse_count = 0
        for record in vcf_reader: 
            parse_count += 1
            LOGGER.debug("vcf parsed record #%d = %s"%(parse_count,str(record)))
            # Do not load non-PASS variants
            if not record.FILTER:
                # If position refers to a haplotype chromosome, ignore
                if record.CHROM in ['chr1','chr2','chr3',
                                'chr4','chr5','chr6','chr7','chr8',
                                'chr9','chr10','chr11','chr12','chr13',
                                'chr14','chr15','chr16','chr17','chr18',
                                'chr19','chr20','chr21','chr22',
                                'chrX','chrY','chrMT']:
                    snpcount += 1
                    yield self._complete_vep_record_parse_and_clean(record,info_headers,csq_headers)
                else:
                    LOGGER.debug("Skipping record due to CHROM=%s pos=%d filter=%s"%(record.CHROM,record.POS,record.FILTER))
            else:
                LOGGER.debug("Skipping record that does not pass filter CHROM=%s pos=%d"%(record.CHROM,record.POS))
        ## We are now allowing Synonymous SNPs to be mapped ##
        LOGGER.info("Exiting PDBMapVEP:yield_completed_vcf_records() after yielding %d Total SNPs (syn+nonsyn) records"%snpcount)
        # print "Nonsynonymous SNPs in %s: %d"%(fname,nscount)


    def vep_records_to_genomic_tables(self,records_from_vep_generator,data_label:str,buffer_size:int = 100):
        """ Uploads genomic data to GenomicData and GenomicConsequence tables
            from VEP record outputs """

        uploaded_records_from_vep_generator = 0
        unwritten_GenomicDataINFOs = []
        unwritten_GenomicConsequenceCSQs = []

        def write_unwritten_GenomicDataINFOs():
            with PDBMapSQLdb() as db:
                nonlocal unwritten_GenomicDataINFOs
                db.activate_dict_cursor()
                db.set_session_transaction_isolation_read_committed()

                insert_GenomicData   = (SQL_EXTENDED_LOCK_WAIT+"INSERT IGNORE INTO GenomicData "
                     "(label,chrom,pos,end,id,variation,vtype,svtype,ref_allele,alt_allele,"
                     "svlen,quality,avgpost,rsq,erate,theta,ldaf,ac,an,aa,da,maf,amr_af,asn_af,"
                     "eas_af,sas_af,afr_af,eur_af,ens_gene,hgnc_gene,"
                     "snpsource,format,gt) VALUES "
                     "(%(LABEL)s,%(CHROM)s,%(POS)s,%(END)s,%(ID)s,%(EXISTING)s,%(VT)s,"
                     "%(SVTYPE)s,%(REF)s,%(ALT)s,%(SVLEN)s,%(QUAL)s,%(AVGPOST)s,%(RSQ)s,"
                     "%(ERATE)s,%(THETA)s,%(LDAF)s,%(AC)s,%(AN)s,%(AA)s,%(DA)s,"
                     "%(AF)s,%(AMR_AF)s,%(ASN_AF)s,%(EAS_AF)s,%(SAS_AF)s,%(AFR_AF)s,%(EUR_AF)s,"
                     "%(GENE)s,%(HGNC)s,%(SNPSOURCE)s,%(FORMAT)s,%(GT)s)")

                rows_inserted = db.executemany(insert_GenomicData,unwritten_GenomicDataINFOs)
                # We flushed the rows so now we re-initialize for more rows
                unwritten_GenomicDataINFOs = []
                return rows_inserted

        def write_unwritten_GenomicConsequenceCSQs():
            with PDBMapSQLdb() as db:
                nonlocal unwritten_GenomicConsequenceCSQs
                db.activate_dict_cursor()
                db.set_session_transaction_isolation_read_committed()

                # Upload each consequence to GenomicConsequence
                insert_GenomicConsequence  = (SQL_EXTENDED_LOCK_WAIT+"INSERT IGNORE INTO GenomicConsequence "
                     "(label,chrom,pos,end,id,transcript,protein,uniprot,canonical,allele,"
                     "consequence,cdna_pos,cds_pos,protein_pos,ref_amino_acid,"
                     "alt_amino_acid,ref_codon,alt_codon,polyphen,sift,biotype,"
                     "domains) VALUES "
                     "(%(LABEL)s,%(CHROM)s,%(POS)s,%(END)s,%(ID)s,"
                     "%(Feature)s,%(ENSP)s,%(SWISSPROT)s,%(CANONICAL)s,%(Allele)s,"
                     "%(Consequence)s,%(cDNA_position)s,%(CDS_position)s,"
                     "%(Protein_position)s,%(Ref_AminoAcid)s,%(Alt_AminoAcid)s,"
                     "%(Ref_Codon)s,%(Alt_Codon)s,%(PolyPhen)s,%(SIFT)s,"
                     "%(BIOTYPE)s,%(DOMAINS)s)")

                rows_inserted = db.executemany(insert_GenomicConsequence,unwritten_GenomicConsequenceCSQs)
                # We flushed the rows so now we re-initialize for more rows
                unwritten_GenomicConsequenceCSQs = []
                return rows_inserted


        for vcf_record in records_from_vep_generator:
            vcf_record.INFO['LABEL'] = data_label
            # Upload each consequence to GenomicData
            unwritten_GenomicDataINFOs += [vcf_record.INFO]

            if (len(unwritten_GenomicDataINFOs) % buffer_size) == 0:
                rows_affected = write_unwritten_GenomicDataINFOs()
            # assert rows_affected <= 1,"Faiure with %s"%(insert_GenomicData%vcf_record.INFO,)
    
            # Iterate over our cleaned up, parsed, supplemented CSQs                
            for csq in vcf_record.CSQ:
                csq["LABEL"] = data_label

                # if csq['Consequence'].find("mis") >= 0:
                #     import pdb; pdb.set_trace()
                unwritten_GenomicConsequenceCSQs += [csq]
                if (len(unwritten_GenomicConsequenceCSQs) % buffer_size) == 0:
                   rows_affected = write_unwritten_GenomicConsequenceCSQs()
                # assert rows_affected <= 1,"Faiure with %s"%(insert_GenomicConsequence%csq)

            uploaded_records_from_vep_generator += 1

        if len(unwritten_GenomicDataINFOs) > 0:
            write_unwritten_GenomicDataINFOs()
        if len(unwritten_GenomicConsequenceCSQs) > 0:
            write_unwritten_GenomicConsequenceCSQs()

        return uploaded_records_from_vep_generator

# Main check
if __name__ == "__main__":
    sys.stderr.write("Class definition. Should not be called from command line.")
    sys.exit(1)
