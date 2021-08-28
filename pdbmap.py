#!/usr/bin/env python3
#
# Project        : PDBMap
# Filename       : PDBMap.py
# Author         : R. Michael Sivley - refactored by Chris Moth 2019 November
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu  chris.moth@vanderbilt.edu
# Date           : 2014-02-17                  2019-11-19
#
# =============================================================================#

"""PDBMap is a command-line interface for loading data
and providing access to the PDBMap* library modules.
PDBMap library modules map codons in the human exome to
their corresponding amino acids in known and predicted
protein structures. Using this two-way mapping, PDBMap
allows for the mapping of genomic annotations onto protein
structures and the mapping of protein structural annotation 
onto genomic elements. The methods to build and utilize this
tool may be called from this master class."""

# See main check for cmd line parsing
import argparse, configparser
import traceback
import sys, os, csv, time, pdb, glob, gzip, shutil, getpass
import subprocess as sp
from multiprocessing import cpu_count
import numpy as np
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from warnings import filterwarnings, resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from lib import PDBMapVEP
from lib.PDBMapSQLdb import PDBMapSQLdb
from lib import PDBMapProtein
from lib import PDBMapAlignment, PDBMapData
from lib import PDBMapModel, PDBMapSwiss
from lib.PDBMapVisualize import PDBMapVisualize
from lib import amino_acids
import logging
from logging.handlers import RotatingFileHandler
from logging import handlers

sh = logging.StreamHandler()
LOGGER = logging.getLogger()
LOGGER.addHandler(sh)

log_format_string = '%(asctime)s %(levelname)-4s [%(filename)16s:%(lineno)d] %(message)s'
date_format_string = '%H:%M:%S'
log_formatter = logging.Formatter(log_format_string, date_format_string)

LOGGER.setLevel(logging.DEBUG)
sh.setLevel(logging.INFO)
sh.setFormatter(log_formatter)


class PDBMap():
    def __init__(self, config_dict):
        self._config_dict = config_dict

        """
      # Initialize
      if idmapping:
        PDBMapProtein.load_idmapping(idmapping)
      if sec2prim:
        PDBMapProtein.load_sec2prim(sec2prim)
      if sprot:
        PDBMapProtein.load_sprot(sprot)
      if pdb_dir:
        self.pdb = True
        self.pdb_dir = pdb_dir
      if swiss_dir and swiss_summary:
        self.swiss = True
        PDBMapSwiss.load_swiss_INDEX_JSON(swiss_dir,swiss_summary)
      if modbase2016_dir and modbase2016_summary:
        self.modbase = True
        PDBMapModel.load_modbase(modbase2016_dir,modbase2016_summary)
      if modbase2013_dir and modbase2013_summary:
        self.modbase = True
        PDBMapModel.load_modbase(modbase2013_dir,modbase2013_summary)
      if vep:
        self.vep       = vep
        self.vep_cache = vep_cache
      if reduce:
        self.reduce = reduce
      if probe:
        self.probe = probe
      """

    def load_unp(self, unp, label=None, use_pdb=True, use_modbase=True, update=False):
        """ Loads all known structures associated with UniProt ID """
        if self.pdb and use_pdb:
            pdb_label = label if label else 'pdb'
            pdb_label = "%s_update" % pdb_label if update else pdb_label
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, slabel=pdb_label)
            pdbids = list(set(PDBMapProtein.unp2pdb(unp)))
            for pdbid in pdbids:
                print(" # Processing (%s) PDB %s # " % (pdb_label, pdbid))
                self.load_pdb(pdbid, label=pdb_label, io=io)
                sys.stdout.flush()  # Force stdout flush after each PDB
        if self.modbase and use_modbase:
            mod_label = label if label else 'modbase'
            mod_label = "%s_update" % mod_label if update else mod_label
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, slabel=mod_label)
            modelids = PDBMapModel.unp2modbase(unp)
            models = [PDBMapModel.get_info(modelid) for modelid in modelids]
            for model in models:
                print(" # (%s) Processing ModBase %s #" % (mod_label, model['modelid']))
                self.load_model(model, label=mod_label, io=io, update=update)
                sys.stdout.flush()  # Force stdout flush after each model
        if not pdbids and not models:
            msg = "  WARNING (PDBMap) No PDB structures or Modbase models found for %s\n" % unp
            LOGGER.warning(msg)

    def load_pdb(self, pdbid, pdb_fname=None, label="", io=None, update=False):
        """ Loads a given PDB into the PDBMap database """
        if not io:
            # Create a PDBMapIO object
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, slabel=label)
        # Check if PDB is already in the database
        if io.structure_in_db(pdbid, label):
            if not update:  # silence if updating
                print("  VALID (PDBMap) %s already in database." % pdbid)
                return 0
        # Load the PDB structure
        if not pdb_fname:
            pdb_fname = "%s/structures/pdb%s.ent.gz" % (self.pdb_dir, pdbid.lower())
            print("  # Fetching %s" % pdbid)
            if not os.path.exists(pdb_fname):
                msg = "  ERROR (PDBMap) Cannot fetch %s. %s Not in PDB mirror.\n" % (pdbid, pdb_fname)
                LOGGER.error(msg)
                return 1
        # Locate all biological assemblies
        biounit_fnames = glob.glob("%s/biounits/%s.pdb*.gz" % (self.pdb_dir, pdbid.lower()))
        try:  # Load the structure
            p = PDBMapParser()
            s = p.get_structure(pdbid, pdb_fname, biounit_fnames=biounit_fnames, io=io)
            io.set_structure(s)
            io.upload_structure()
        except Exception as e:
            msg = "  ERROR (PDBMap) %s: %s\n\n" % (pdbid, str(e))
            LOGGER.exception(msg)
            return 1
        msg = "  VALID (PDBMap) %s complete.\n" % pdbid
        LOGGER.info(msg)
        return 0

    def load_swiss_to_MySQL(self, modelid, label="", io=None):
        if not io:
            # Create a PDBMapIO object
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, slabel=label)

        # Load the dictionary of information about the modelid
        model_summary = PDBMapSwiss.get_info(modelid)

        # Check if model is already in the database
        model_fname = PDBMapSwiss.get_coord_file(modelid);
        if io.swiss_in_db(modelid, label):
            msg = "  VALID (SWISS) %s (%s) already in database.\n" % (modelid, label)
            print(msg)
            return 0
        print("Will attempt to add new swiss %s.  Fetching: %s" % (modelid, model_fname))
        if not os.path.exists(model_fname):
            model_fname += '.gz'  # check for compressed copy
        if not os.path.exists(model_fname):
            msg = "  ERROR (load_swiss_to_MySQL) %s not in local Swiss mirror.\nExpected file %s\n" % (
            modelid, model_fname)
            LOGGER.error(msg)
            return 1

        remark3_metrics = PDBMapSwiss.load_REMARK3_metrics(modelid)

        try:
            s = PDBMapParser.getBiopythonStructureOrFail(modelid, model_fname)
            m = PDBMapSwiss(s, model_summary)
            s = PDBMapParser.process_structure_dssp_unp2hgnc(m, model_summary, model_fname, m.unp)

            io.set_structure(m)
            io.upload_swiss(remark3_metrics)
        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            emsg = str(e)
            print(emsg)
            msg = "  ERROR (pdbmap.py: load_swiss_to_MySQL(%s):\n%s" % (modelid, emsg)
            LOGGER.error(msg)
            traceback.print_tb(exc_traceback, limit=2, file=sys.stderr)
            return 1
        msg = "  VALID (pdbmap.py: load_swiss_to_MySQL(%s)\n" % modelid
        LOGGER.info(msg)
        return 0

    def load_model(self, model_summary, label="", io=None, update=False):
        """ Loads a given ModBase model into the PDBMap database """

        if not io:
            # Create a PDBMapIO object
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, slabel=label)

        # Check if model is already in the database
        modelid = model_summary['modelid']  # extract ModBase model ID
        model_fname = model_summary['filename']
        if io.model_in_db(modelid, label):
            if not update:  # silence if updating
                print("  VALID (PDBMap) %s (%s) already in database.\n" % (modelid, label))
                return 0

        # Query UniProt ID if not provided
        if 'unp' not in model_summary:
            unp = PDBMapProtein.ensp2unp(modelid.split('.')[0].split('_')[0])[0]
        else:
            unp = model_summary['unp']

        # Load the ModBase model
        if not model_fname:
            model_fname = PDBMapModel.get_coord_file(modelid.upper())
            print("  # Fetching %s" % modelid)
            if not os.path.exists(model_fname):
                model_fname += '.gz'  # check for compressed copy
            if not os.path.exists(model_fname):
                msg = "  ERROR (load_model) %s not in ModBase mirror.\n" % modelid
                LOGGER.error(msg)
                return 1
        try:
            p = PDBMapParser()
            print("   # Loading %s (%s) from %s..." % (modelid, unp, model_fname.split('/')[-1]))
            m = p.get_model(model_summary, model_fname, unp=unp)
            io.set_structure(m)
            io.upload_model()
        except Exception as e:
            msg = "  ERROR (load_model) %s: %s\n" % (modelid, str(e))
            LOGGER.error(msg)
            return 1
        msg = "  VALID (load_model) %s complete.\n" % modelid
        LOGGER.info(msg)
        return 0

    def load_vcf_data(self, dlabel, vcf_filename):
        """ Load a single vcf format file into a sql table
            Parameters:
                dlabel:       Label for sql table name (ex. clinvar)
                vcf_filename: VCF format input file
                
            Dry-run testing parameters:
                usevep:       Set to false for dry-run testing and VEP will not be run
                upload_to_sql Set to false if the SQL tables should not be created """

        if not os.path.exists(vcf_filename):
            msg = 'vcf input file %s not found' % vcf_filename
            LOGGER.critical(msg)
            sys.exit(msg)

        pdbmap_vep = PDBMapVEP(self._config_dict)

        PDBMapSQLdb.set_access_dictionary(self._config_dict)
        nrows_inserted = 0
        sql_table_name = dlabel  # <-- for now - do more later
        nrows_inserted = pdbmap_vep.vcf_file_to_sql_table(sql_table_name, vcf_filename,
                                                          args.buffer_size if args.buffer_size else 1)
        return nrows_inserted

    def vcf_data_to_vep_to_genomic_tables(self, dlabel, vcf_filename):
        """Send a single vcf format file to 
            Parameters:
                dlabel:       Label for label column in GenomicData and GenomicConsequence tables (ex. clinvar)
                vcf_filename: VCF format input file that will be read, and optionally additionally annotated by ENSEMBL VEP
                usevep:       Set to False when the CSQ annotations in the source file are sufficient, and vep run unnecesssary (example EXAC)"""

        if not os.path.exists(vcf_filename):
            msg = 'vcf input file %s not found' % vcf_filename
            LOGGER.critical(msg)
            sys.exit(msg)

        pdbmap_vep = PDBMapVEP(self._config_dict)  ## Need to specify vep_cache
        if args.novep:
            # The caller feels that the CSQs are adequence without adding vep annotations
            LOGGER.info("Calling pdbmap_vep.vcf_reader_from_file_without_vep(%s)", vcf_filename)
            vcf_reader = pdbmap_vep.vcf_reader_from_file_without_vep(vcf_filename)
        else:
            # Feed the caller's file to VEP, and parse the file with added vep annotations
            vep_echo_filename = os.path.join(os.path.dirname(args.logfile), os.path.basename(vcf_filename) + ".vep.out")
            LOGGER.info("Calling pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(%s,%s)", vcf_filename,
                        vep_echo_filename)
            vcf_reader = pdbmap_vep.vcf_reader_from_file_supplemented_with_vep_outputs(vcf_filename, vep_echo_filename)

        completed_vep_records_generator = pdbmap_vep.yield_completed_vcf_records(vcf_reader)

        LOGGER.info("Calling PDBMapVEP.vep_records_to_genomic_tables")
        uploaded_records_count = pdbmap_vep.vep_records_to_genomic_tables(completed_vep_records_generator, dlabel)
        LOGGER.info(
            "PDBMapVEP.pipe_vcf_to_vep_and_yield_complete_records returned %d (records uploaded)" % uploaded_records_count)

    def load_data(self, dname, dfile, indexing=None, usevep=True, upload=True):
        """ Loads a data file into the PDBMap database """
        if usevep:
            d = PDBMapData(vep=self.vep, vep_cache=self.vep_cache, dname=dname)
        else:
            d = PDBMapData(dname=dname)
        if not os.path.exists(dfile):
            dfile = "%s.ped" % dfile  # Test if PEDMAP basename
            if not os.path.exists(dfile):
                msg = "  ERROR (PDBMap) File does not exist: %s" % dfile
                raise Exception(msg)
        io = PDBMapIO(args.dbhost, args.dbuser, args.dbpass, args.dbname, dlabel=dname)
        # Determine file type
        ext = dfile.split('.')[-1].lower()
        if ext == 'gz':
            ext = dfile.split('.')[-2].lower()
        # Process and accordingly
        if ext == 'vcf':
            if upload:
                print("\nUploading VCF to supplemental database...")
                nrows = d.load_vcffile(dfile, io, args.buffer_size)
                print("%d VCF records uploaded to supplemental database before processing" % nrows)
            generator = d.load_vcf(dfile, usevep)
        elif ext in ["bed", "txt", "csv"]:
            if usevep:
                print("\nNote: You have provided a %s file and requested VEP analysis." % ext)
                print("      There are three acceptable values for the 'name' column.")
                print("       1. REF/ALT - SNP alleles will be used as input to VEP")
                print("       2. rsID    - SNP names will be used as input to VEP.")
                print("       3. HGVS    - SNPs HGVS will be used as input to VEP")
                print("      We highly recommend option 1 when possible. Option 2 may")
                print("       exclude rare or otherwise unlabeled SNPs.\n")
            # Determine the column delimiter by file type
            delim = '\t'
            if ext != "bed":
                delim = ' ' if ext == 'txt' else ','
            if not indexing:
                indexing = 'ucsc' if ext == 'bed' else 'pdbmap'
            print("Using %s indexing for %s." % (indexing, dfile))
            dfile, id_type = d.load_bedfile(dfile, io, delim, indexing, usevep)
            print("Creating BED generator...")
            generator = d.load_bed(dfile, id_type, usevep, indexing)
        elif ext in ["ped", "map"]:
            generator = d.load_pedmap(dfile)
        else:
            msg = "  ERROR (PDBMap) Unsupported file type: %s" % ext
            raise Exception(msg)
        # Pass the relevant generator to be uploaded
        nrows = io.upload_genomic_data(generator, dname)
        return (nrows)

    def intersect_data(self, dname, slabel=None, dtype="Genomic", quick=False):
        """ Intersects a loaded dataset with the PDBMap structural domain """
        io = PDBMapIO(args.dbhost, args.dbuser, args.dbpass, args.dbname, dlabel=dname, slabel=slabel)
        i = PDBMapIntersect(io)
        # Note: Only all-structures <-> genomic data intersections supported
        if quick:
            nrows = i.quick_intersect(dname, slabel, dtype)
        else:
            nrows = i.intersect(dname, slabel, dtype, args.buffer_size)
        return (nrows)  # Return the number of intersections

    def visualize(self, entity, biounits=[], struct_label='pdb',
                  data_label='1kg', anno_list=['maf'], spectrum_range=[], colors=[]):
        """ Visualizes a PDBMap structure, model, or protein """
        io = PDBMapIO(args.dbhost, args.dbuser, args.dbpass, args.dbname,
                      slabel=struct_label, dlabel=data_label)
        v = PDBMapVisualize(io, args.pdb_dir)
        entity_type = io.detect_entity_type(entity) if not entity == 'all' else 'all'
        if entity_type == 'structure' and not biounits:
            if io.is_nmr(entity):
                biounits = [-1]
            else:
                # Query all biological assemblies, exclude the asymmetric unit
                query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
                res = io.secure_query(query, (struct_label, entity,), cursorclass='Cursor')
                biounits = [r[0] for r in res]
        elif entity_type == 'model' and not biounits:
            biounits = [-1]
        eps, mins = False, False
        synonymous_flag = False
        if any(['.synonymous' in a for a in anno_list]):
            # Replace synonymous with DAF and set the flag
            synonymous_flag = True
            idx = ['.synonymous' in a for i, a in enumerate(anno_list)].index(True)
            anno_list[idx] = anno_list[idx].replace('.synonymous', '')
            print("\n%s will be plotted for synonymous variants." % anno_list[idx])
        if 'popdaf' in anno_list:
            idx = anno_list.index('popdaf')
            anno_list = anno_list[0:idx] + anno_list[idx + 1:]
            anno_list += ['daf', 'amr_daf', 'eas_daf', 'sas_daf', 'afr_daf', 'eur_daf']
            sr = spectrum_range[idx]
            spectrum_range = spectrum_range[0:idx] + spectrum_range[idx + 1:]
            spectrum_range += [sr for i in range(6)]
        if 'popmaf' in anno_list:
            idx = anno_list.index('popmaf')
            anno_list = anno_list[0:idx] + anno_list[idx + 1:]
            anno_list += ['maf', 'amr_af', 'eas_af', 'sas_af', 'afr_af', 'eur_af']
            sr = spectrum_range[idx]
            spectrum_range = spectrum_range[0:idx] + spectrum_range[idx + 1:]
            spectrum_range += [sr for i in range(6)]
        if 'dbscan' in anno_list:
            idx = anno_list.index('dbscan')
            anno_list = anno_list[0:idx] + anno_list[idx + 1:]
            eps, mins = spectrum_range[idx]
            spectrum_range = spectrum_range[0:idx] + spectrum_range[idx + 1:]
            if len(anno_list):  # more than one DBSCAN specification
                msg = "ERROR (PDBMap) Cannot run other annotations with DBSCAN"
                raise Exception(msg)
        try:
            if entity_type in ['structure', 'model']:
                for biounit in biounits:
                    v.visualize_structure(entity, biounit, anno_list, eps, mins, spectrum_range, colors=colors,
                                          syn=synonymous_flag)
            elif entity_type == 'unp':
                v.visualize_unp(entity, anno_list, eps, mins, spectrum_range, colors=colors, syn=synonymous_flag)
            elif entity_type == 'all':
                v.visualize_all(anno_list, eps, mins, spectrum_range, colors=colors, syn=synonymous_flag)
            elif entity_type:
                print("%s matched with UniProt ID: %s" % (entity.upper(), entity_type))
                entity = entity_type  # An HGNC ID was detected and converted to UNP ID
                v.visualize_unp(entity, anno_list, eps, mins, spectrum_range, colors=colors, syn=synonymous_flag)
            else:
                msg = "Sorry, but the specified entity is not in the PDBMap database.\n"
                LOGGER.error(msg)
                return 1
        except Exception as e:
            msg = "ERROR (PDBMap) Visualization failed: %s" % str(e)
            raise

    def summarize(self):
        """ Returns summary statistics for the PDBMap database """
        print("Basic summary statistics for PDBMap. Not implemented.")

    def refresh_cache(self, args, io, force=False):
        """ Refreshes all mirrored data """
        if args.sifts:
            print("Refreshing local SIFTS cache...", end=' ')
            if os.path.exists(args.sifts):
                # Record the most recent modification time for the SIFTS directory
                mtime = os.stat(args.sifts)[-2]
            else:
                mtime = None
            script_path = os.path.dirname(os.path.realpath(args.sifts))
            get_sifts = "cd %s; ./get_sifts.sh" % (script_path)
            os.system(get_sifts)
            # Upload if any modifications were made to the SIFTS directory
            if not mtime or mtime != os.stat(args.sifts)[-2]:
                print("  Updating SIFTS in PDBMap...", )
                sys.stdout.flush()  # flush buffer before long query
                rc = io.load_sifts(args.sifts, args.conf_file)
                print(rc)
                print("%s rows added" % "{:,}".format(int(rc)))
        if args.sprot:
            print("Refreshing local SwissProt cache...")
            script_path = os.path.dirname(os.path.realpath(args.sprot))
            get_sprot = "cd %s; ./get_swissprot.sh" % (script_path)
            os.system(get_sprot)
        if args.idmapping:
            print("Refreshing local UniProt ID Mapping cache...")
            script_path = os.path.dirname(os.path.realpath(args.idmapping))
            get_idmapping = "cd %s; ./get_idmapping.sh" % (script_path)
            os.system(get_idmapping)
        if args.pdb_dir:
            print("Refreshing local PDB cache...")
            script_path = os.path.realpath(args.pdb_dir)
            get_pdb = "cd %s; ./get_pdb.sh" % (script_path)
            os.system(get_pdb)
        if args.modbase2016_dir:
            print("Refreshing local ModBase 2016 cache...")
            script_path = os.path.realpath(args.modbase2016_dir)
            get_modbase = "cd %s; ./get_modbase_2016.sh" % (script_path)
            os.system(get_modbase)
        if args.modbase2013_dir:
            print("Refreshing local ModBase 2013 cache...")
            script_path = os.path.realpath(args.modbase2013_dir)
            get_modbase = "cd %s; ./get_modbase_2013.sh" % (script_path)
            os.system(get_modbase)


class Command_runner:
    """Perform the command-line requested command."""

    @staticmethod
    def load_vcf_cmd(pdbmap, args, remaining_args):
        """pdbmap: A pre-initialized PDBMap instance
        args: The original command-line arguments
        remaining_args.  The still unparsed-arguments"""

        def usage_failure():
            msg = "usage: pdbmap.py load_vcf <data_file> <data_label> [data_file data_label] ...\n"
            msg += "alt:   pdbmap.py load_vcf --dlabel=<data_name> <data_file> [data_file] ...\n"
            LOGGER.critical(msg)
            sys.exit(msg)

        if len(remaining_args) < 1:
            usage_failure()

        if args.dlabel:  # In the "alt:" usage, we pair the --dlabel with each data_file
            dfiles = list(zip(remaining_args, [args.dlabel] * len(remaining_args)))
        else:
            # extract all the filename/label pairs
            if len(remaining_args) % 2:
                LOGGER.critical("The odd count of remaining arguments is not right for data_file data_label pairs")
                usage_failure()

            data_file_names = remaining_args[0::2]
            data_labels = remaining_args[1::2]
            dfiles = list(zip(data_file_names, data_labels))

        LOGGER.info("Processing the following filename/data_label pairs:\n%s" % str(dfiles))

        if args.noupload:
            LOGGER.info("Skipping creation of sql mirror table %s because args.noupload is set")
        else:
            nrows = 0
            for vcf_filename, dlabel in dfiles:
                LOGGER.info("Calling pdbmap.load_vcf_data(%s,%s)" % (dlabel, vcf_filename))
                nrows += pdbmap.load_vcf_data(dlabel, vcf_filename)
            LOGGER.info("%d data rows were uploaded from %s" % (nrows, str(dfiles)))

        LOGGER.info(
            "Processing vcf data in ENSEMBL VEP.  Adding completed VEP outputs to GenomicConsequence and GenomicData tables")
        nrows = 0
        for vcf_filename, dlabel in dfiles:
            pdbmap.vcf_data_to_vep_to_genomic_tables(dlabel, vcf_filename)


# Command line usage
if __name__ == "__main__":

    # Print the ASCII header
    header = """ \
__  __  __            
|__)|  \|__)|\/| _  _  
|   |__/|__)|  |(_||_) 
                   |   """
    print(header)

    first_cmdline_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # first_cmdline_parser.add_argument("cmd",
    #    type=str,
    #    action = 'store',
    #    help="Command choice")
    #    choices=['refresh', 'load_pdb', 'load_unp', 'load_data', 'load_vcf',
    #             'load_swiss', 'delete_swiss', 'intersect', 'visualize','confighelp'])
    default_config_filename = os.path.join("..", "config", "global.config")
    first_cmdline_parser.add_argument("-c", "--config", default=default_config_filename,
                                      help="Global config file", metavar="FILE", type=str)
    first_cmdline_parser.add_argument("-u", "--userconfig", default=getpass.getuser() + ".config",
                                      help="User specific config options to selectively override global config",
                                      metavar="FILE")
    first_cmdline_parser.add_argument(
        "-v", "--verbose",
        help="Include routine info log entries on stderr", default=True, action="store_true")
    first_cmdline_parser.add_argument(
        "-d", "--debug",
        help="Include routine info AND 'debug' log entries on stderr",
        default=False, action="store_true")

    first_cmdline_parser.add_argument(
        "--logfile",
        help="Logfile name override.  Useful for slurm runs where logfiles of each session would otherwise overwrite each other.",
        default=None)

    subparsers = first_cmdline_parser.add_subparsers(help='add --help after cmd for sub-command specific help')
    subparser_load_vcf = subparsers.add_parser('load_vcf', help='vcf help')  # ,action='store_true')
    subparser_load_vcf.add_argument("--dlabel", help='optional dataset label, ex. clinvar')
    subparser_load_vcf.add_argument("--novep", action='store_true',
                                    help="Disables VEP consequence prediction. If no CSQ provided, all SNPs uploaded.")
    subparser_load_vcf.add_argument("--buffer_size", type=int, default=5000,
                                    help="Size of MariaDB buffer (in rows/records) when applicable")
    subparser_load_vcf.add_argument("--noupload", action='store_true',
                                    help="Disables upload of the original data file to a supplementary database. Potential information loss.")
    subparser_load_vcf.set_defaults(func=Command_runner.load_vcf_cmd)

    subparser_confighelp = subparsers.add_parser('confighelp',
                                                 help='List all the optional config file override parameters')
    subparser_confighelp.set_defaults(confighelponly=True)
    # subparser_confighelp.add_argument("--xlabel",help='xoptional dataset label, ex. clinvar')

    ###########################################################################
    # Parse the typical command line parameters
    # Save the command-specific additional parameters in remaining_args
    ###########################################################################
    args, remaining_args = first_cmdline_parser.parse_known_args()

    if not args.logfile:
        args.logfile = "pdbmap." + args.func.__name__ + ".log"

    LOGGER.info("Log will be echoed to file %s", args.logfile)
    needRoll = os.path.isfile(args.logfile)

    fh = RotatingFileHandler(args.logfile, maxBytes=(1048576 * 10), backupCount=5)

    fh.setFormatter(log_formatter)

    if args.debug:
        fh.setLevel(logging.DEBUG)
        sh.setLevel(logging.DEBUG)
    elif args.verbose:
        fh.setLevel(logging.INFO)
        sh.setLevel(logging.INFO)
    else:
        fh.setLevel(logging.INFO)

    LOGGER.addHandler(fh)

    if needRoll:
        fh.doRollover()

    # We allow config file parameterization to be overridden at command line
    # That happens down below.  For now, declare the dict of available
    # comamnd line overrides, and dump to user if they just want confighelp
    config_file_duplicate_parameters = {
        'dbhost': "MariaDB server hostname/ip address",
        'dbuser': "MariaDB username",
        'dbpass': "MariaDB password",
        'dbname': "Database name to 'USE'",

        "pdb_dir": "Directory containing PDB structures",
        "swiss_dir": "Directory containing Swissmodel models",
        "swiss_summary": "Swiss Model summary JSON file",
        "modbase2016_dir": "Directory containing ModBase 2016 models",
        "modbase2016_summary": "ModBase 2016 summary file",
        "modbase2013_dir": "Directory containing ModBase 2013 models",
        "modbase2013_summary": "ModBase 2013 summary file",
        "idmapping": "UniProt ID -> EnsEMBL transcript map file location",
        "sec2prim": "UniProt secondary -> primary AC map file location",
        "sprot": "Swiss-Prot file location",
        "vep": "Variant Effect Predictor location",
        "vep_cache": "Directory containing the VEP cache"
    }

    if 'confighelponly' in vars(args):
        print("Configutation settings inform PDBMap where to locate data files, and")
        print("how to access the database.  You provide options in 3 ways.  Either")
        print("via global options in the -c global.config file, overrides in")
        print("the -u user.config file, or lastly by available command line overrides:")
        for parameter_key_and_help in config_file_duplicate_parameters.items():
            print("  --%-19s %s" % parameter_key_and_help)
        sys.exit(0)

    ###########################################################################
    # Create a Config File Parser to read the -c and -u command line 
    # requested config files
    # In most cases, the config files provide the complete dictionary
    # needed to initialize class PDBMap().  However, command line overrides
    # are also (code below) possible
    ###########################################################################
    config_file_parser = configparser.ConfigParser(allow_no_value=True)

    config_file_list = []
    if args.config:
        config_file_list += [args.config]
    if args.userconfig:
        config_file_list += [args.userconfig]

    LOGGER.info("Attempting to read config from " + str(config_file_list))
    parsed_file_list = config_file_parser.read(config_file_list)
    LOGGER.info("Successfully read config from " + str(parsed_file_list))
    if len(parsed_file_list) < 1:
        LOGGER.warning("""No -c or -u config files were processed.  Resources must be specified on command line""")

    # config.items() returns a list of (name, value) pairs - convert to dict
    config_dict = dict(config_file_parser.items("Genome_PDB_Mapper"))

    # It is an odd thing that we poulate from a UserSpecific section - fix later
    user_specific_overrides = None
    try:
        user_specific_overrides = config_file_parser.items("UserSpecific")
    except configparser.NoSectionError:
        LOGGER.info("No UserSpecific section was found in the configuration files")

    if user_specific_overrides:
        config_dict.update(dict(user_specific_overrides))

    # In addition to the global config, and user config files,
    # config dictionary settings may be overridden with command line options
    # Gather those now.  End-user help for these is given with the "confighelp" command
    # and not otherwise - to avoid clutter

    LOGGER.info("Parsing any config file overrides given on the command line")

    overrides_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    for dict_key, description in config_file_duplicate_parameters.items():
        overrides_parser.add_argument("--" + dict_key, help=description, type=str)

    override_args, override_remaining_args = overrides_parser.parse_known_args(remaining_args)
    # Make a dictionary of command line overrides, filter out non-overridden arguments
    override_args_dict = {dict_key: value \
                          for dict_key, value in vars(override_args).items() if value is not None}

    LOGGER.info("%d command line configuration overrides" % len(override_args_dict))

    if len(override_args_dict) > 0:
        config_dict.update(override_args_dict)

    LOGGER.info("Calling command func=%s with remaining args: %s" % (str(args.func), str(remaining_args)))

    pdbmap = PDBMap(config_dict)

    args.func(pdbmap, args, remaining_args)

    sys.exit(1)

    if args.load_vcf:
        print("loading vcf")
        sys.exit("loading vcf")

    """
    if required_config_items:
      missing_keys = [name for name in required_config_iteDoDwu/
drwx------ 3 root   root             17 Oct 17 16:00 systemd-private-1461403187af4efe9331fee92595b522-chronyd.servicems if name not in config_dict]

      if len(missing_keys):
        LOGGER.error("Can't proceed without configuration file options set for: %s", str(missing_keys))
        sys.exit(1)
    """
    """
    if args.conf_filecmdline:
      config = configparser.ConfigParser()
      config.read([args.conf_file])
      defaults.update(dict(config.items("Genome_PDB_Mapper")))
    conf_file = args.conf_file
    """

    # Setup the Command Line Parser
    parser = argparse.ArgumentParser(parents=[first_cmdline_parser], description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    defaults['create_new_db'] = ('True' == defaults['create_new_db'])
    parser.add_argument("-v", "--version", action="version",
                        version="PDBMap version 1.8")

    for duplicate_parameter in config_file_duplicate_parameters:
        parser.add_argument("--" + duplicate_parameter,
                            help=argparse.SUPPRESS
                            )

    parser.add_argument("--create_new_db", action='store_true',
                        help="Create a new database prior to uploading PDB data")
    parser.add_argument("-f", "--force", action='store_true',
                        help="Force configuration. No safety checks.")
    parser.add_argument("--pdbid",
                        help="Force pdbid")
    parser.add_argument("--unp",
                        help="Force unp")
    parser.add_argument("--slabel",
                        help="Structural label for this session")
    parser.add_argument("--dlabel",
                        help="Data label for this session")
    parser.add_argument("--indexing",
                        help="Indexing used by data file(s) [pdbmap,ensembl,ucsc]")
    parser.add_argument("--buffer_size", type=int,
                        help="Size of mysql buffer (in rows/records) when applicable")
    parser.add_argument("--ppart", type=int,
                        help="Used to manage parallel subprocesses. Do not call directly.")
    parser.add_argument("--ppidx", type=int,
                        help="Used to manage parallel subprocesses. Do not call directly.")
    parser.add_argument("-j", "--cores", type=int,
                        help="Number of available processors")
    parser.set_defaults(**defaults)

    args = parser.parse_args(remaining_argv)

    args.conf_file = conf_file
    parser.get_default("vep")
    args.create_new_db = bool(args.create_new_db)
    args.force = bool(args.force)
    args.cores = int(args.cores)

    if args.create_new_db and not args.force:
        print("You have opted to create a new database: %s." % args.dbname)
        if input("Are you sure you want to do this? (y/n):") != 'y':
            print("Aborting...")
        else:
            print("Creating database tables...")
            io = PDBMapIO(args.dbhost, args.dbuser,
                          args.dbpass, args.dbname, createdb=True)
            print("\nDatabase created. Please set create_new_db to False.")
            print("\nIt is strongly recommended that you now refresh the local resource cache.")
            if input("Would you like to refresh the cache now? (y/n):") == 'y':
                print("Refreshing local cache...")
                args.cmd = "refresh"  # continue on to cache refresh
            else:
                sys.exit(0)

    # Initialize PDBMap, refresh mirrored data if specified
    if args.cmd == "refresh":
        io = PDBMapIO(args.dbhost, args.dbuser,
                      args.dbpass, args.dbname)
        try:
            PDBMap().refresh_cache(args, io)
            print("\nRefresh complete. Exiting...")
            sys.exit(0)
        except:
            print("\nAn unknown error occurred. Terminating...")
            sys.exit(1)

    ## load_pdb ##
    elif args.cmd == "load_pdb":
        pdbmap = PDBMap(idmapping=args.idmapping, sec2prim=args.sec2prim,
                        pdb_dir=args.pdb_dir, sprot=args.sprot)
        if len(args.args) < 1:
            msg = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_pdb pdb_file [pdb_file,...]\n"
            msg += "   or: pdbmap.py -c conf_file load_pdb all\n"
            msg += "   or: pdbmap.py -c conf_file load_pdb update"
            print(msg);
            sys.exit(0)
        elif args.args[0] in ('all', 'pdb'):
            update = True if args.args[0] == "update" else False
            if args.slabel:
                msg = "WARNING (PDBMap) PDB all/update does not support custom structure labels.\n"
                LOGGER.warning(msg)
                args.slabel = None
            args.slabel = args.slabel if args.slabel else "pdb"
            # All human Swiss-Prot-containing PDB structures in the local mirror
            all_pdb_ids = list(PDBMapProtein._pdb2unp.keys())
            fname = "%s/structures/pdb%s.ent.gz"
            all_pdb_files = [fname % (args.pdb_dir, pdbid.lower()) for pdbid in all_pdb_ids]
            # Remove any PDB files not contained in the local PDB mirror
            all_pdb_files = [f for f in all_pdb_files if os.path.exists(f)]
            msg = "WARNING (PDBMap) Uploading %d Human Swiss-Prot PDB structures.\n" % len(all_pdb_files)
            LOGGER.warning(msg)
            n = len(all_pdb_files)
            # If this is a parallel command with partition parameters
            if args.ppart != None and args.ppidx != None:
                psize = n / args.ppart  # floor
                if (args.ppart - 1) == args.ppidx:
                    all_pdb_files = all_pdb_files[args.ppidx * psize:]
                else:
                    all_pdb_files = all_pdb_files[args.ppidx * psize:(args.ppidx + 1) * psize]
                msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of PDB\n" % (args.ppidx + 1, args.ppart)
                LOGGER.warning(msg)
                for i, pdb_file in enumerate(all_pdb_files):
                    pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
                    print("## Processing (%s) %s (%d/%d) ##" % (args.slabel, pdbid, i + (args.ppidx * psize) + 1, n))
                    pdbmap.load_pdb(pdbid, pdb_file, label=args.slabel, update=update)
            else:
                msg = "WARNING(PDBMap) Uploading all %d PDB IDs.\n" % n
                LOGGER.warning(msg)
                for i, pdb_file in enumerate(all_pdb_files):
                    pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()

                    print("## Processing (pdb) %s (%d/%d) ##" % (pdbid, i, n))
                    pdbmap.load_pdb(pdbid, pdb_file, label=args.slabel, update=update)
        elif len(args.args) == 1:
            # Process one PDB
            pdb_file = args.args[0].strip()
            if not os.path.exists(pdb_file):
                # Not a file, its a PDB ID
                args.pdbid = pdb_file
                # ppProtein._pdb2unp.keys()
                pdb_file = None
            elif not args.pdbid:
                args.pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
            if not args.slabel:
                args.slabel = "manual"
            print("## Processing (%s) %s ##" % (args.slabel, args.pdbid))
            pdbmap.load_pdb(args.pdbid, pdb_file, label=args.slabel)
        else:
            # Process many PDB IDs
            pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(), pdb_file) for pdb_file in args.args]
            n = len(pdbs)
            for i, (pdbid, pdb_file) in enumerate(pdbs):
                if not os.path.exists(pdb_file):
                    pdb_file = None  # Not a file, its a PDB ID
                if not args.slabel:
                    args.slabel = "manual"
                print("## Processing (%s) %s (%d/%d) ##" % (args.slabel, pdbid, i, n))
                pdbmap.load_pdb(pdbid, pdb_file, label=args.slabel)

    ## load_unp ##
    elif args.cmd == "load_unp":
        pdbmap = PDBMap(idmapping=args.idmapping, sec2prim=args.sec2prim,
                        sprot=args.sprot, pdb_dir=args.pdb_dir,
                        modbase2016_dir=args.modbase2016_dir,
                        modbase2016_summary=args.modbase2016_summary,
                        modbase2013_dir=args.modbase2013_dir,
                        modbase2013_summary=args.modbase2013_summary)
        if len(args.args) < 1:
            msg = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_unp unpid [unpid,...]\n"
            msg += "   or: pdbmap.py -c conf_file load_unp all\n"
            msg += "   or: pdbmap.py -c conf_file load_unp update"
            print(msg);
            sys.exit(0)
        elif args.args[0] in ('all', 'update'):
            update = True if args.args[0] == "update" else False
            if args.slabel:
                msg = "WARNING (PDBMap) UniProt all/update does not support custom structure labels.\n"
                LOGGER.warning(msg)
                args.slabel = None
            # All Human, Swiss-Prot proteins with EnsEMBL transcript cross-references
            print("\nIdentifying all Swiss-Prot IDs with mapped Ensembl transcripts...")
            all_unp = [unp for unp in PDBMapProtein._unp2enst \
                       if unp in PDBMapProtein.sprot and \
                       PDBMapProtein.unp2species[unp] == "HUMAN"]
            n = len(all_unp)
            print("Total Swiss-Prot IDs to process: %d. Beginning..." % n)
            # If this is a parallel command with partition parameters
            if args.ppart != None and args.ppidx != None:
                psize = n / args.ppart  # floor
                if (args.ppart - 1) == args.ppidx:
                    all_unp = all_unp[args.ppidx * psize:]
                else:
                    all_unp = all_unp[args.ppidx * psize:(args.ppidx + 1) * psize]
                msg = "WARNING (PDBMap) Subprocess uploading partition %d/%d of Swiss-Prot\n" % (
                    args.ppidx + 1, args.ppart)
                LOGGER.warning(msg)
                for i, unp in enumerate(all_unp):
                    print("\n## Processing (%s) %s (%d/%d) ##" % (args.slabel, unp, i + (args.ppidx * psize) + 1, n))
                    pdbmap.load_unp(unp, label=args.slabel, update=update)
            # This is a standard, full-set load_unp command
            else:
                msg = "WARNING (PDBMap) Uploading all %d Swiss-Prot UniProt IDs.\n" % n
                LOGGER.warning(msg)
                for i, unp in enumerate(all_unp):
                    print("\n## Processing (%s) %s (%d/%d) ##" % (args.slabel, unp, i, n))
                    pdbmap.load_unp(unp, label=args.slabel, update=update)
        elif len(args.args) == 1:
            # Process one UniProt ID
            unp = args.args[0]
            print("\n## Processing (%s) %s ##" % (args.slabel, unp))
            pdbmap.load_unp(unp, label=args.slabel)
        else:
            # Process many UniProt IDs
            n = len(args.args)
            for i, unp in enumerate(args.args):
                print("\n## Processing (%s) %s (%d/%d) ##" % (args.slabel, unp, i, n))
                pdbmap.load_unp(unp, label=args.slabel)

    ## load_model ##
    elif args.cmd == "delete_swiss":
        io = PDBMapIO(args.dbhost, args.dbuser,
                      args.dbpass, args.dbname, slabel='swiss')
        io.delete_all_swiss()
    elif args.cmd == "load_swiss":
        pdbmap = PDBMap(idmapping=args.idmapping, sec2prim=args.sec2prim,
                        sprot=args.sprot,
                        swiss_dir=args.swiss_dir,
                        swiss_summary=args.swiss_summary)
        if len(args.args) < 1:
            msg = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_swiss swiss_summary summary model[,model,...]\n"
            msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_swiss swiss_summary modeldir/*\n"
            msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_swiss all"
            print(msg);
            sys.exit(1)
        if args.args[0] in ['all', '*', '.']:
            args.slabel = args.slabel if args.slabel else "swiss"
            models = PDBMapSwiss.get_swiss_modelids()  # get all  the swiss model unique IDs
            n = len(models)
            # If this is a parallel command with partition parameters
            if args.ppart != None and args.ppidx != None:
                psize = n / args.ppart  # floor
                if (args.ppart - 1) == args.ppidx:
                    models = models[args.ppidx * psize:]
                else:
                    models = models[args.ppidx * psize:(args.ppidx + 1) * psize]
                msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of Swiss Model\n" % (
                    args.ppidx + 1, args.ppart)
                LOGGER.warning(msg)
                for i, modelid in enumerate(models):
                    print("## Processing (%s) %s (%d/%d) ##" % (
                        args.slabel, swissUniqueIdentifier, i + (args.ppidx * psize) + 1, n))
                    pdbmap.load_swiss_to_MySQL(modelid, label=args.slabel, io=None)
            else:
                msg = "WARNING (PDBMap) Uploading all %d Swiss Model  models.\n" % n
                LOGGER.warning(msg)
                for i, modelid in enumerate(models):
                    print("\n## Processing (%s) %s (%d/%d) ##" % (args.slabel, modelid, i, n))
                    pdbmap.load_swiss_to_MySQL(modelid, label=args.slabel, io=None)


        else:
            print("Code for individual swiss models must be done later")
            sys.exit()

    ## load_model ##
    elif args.cmd == "load_model":
        pdbmap = PDBMap(idmapping=args.idmapping, sec2prim=args.sec2prim,
                        sprot=args.sprot,
                        modbase2016_dir=args.modbase2016_dir,
                        modbase2016_summary=args.modbase2016_summary,
                        modbase2013_dir=args.modbase2013_dir,
                        modbase2013_summary=args.modbase2013_summary)
        if len(args.args) < 1:
            msg = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary model[,model,...]\n"
            msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary modeldir/*\n"
            msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model all"
            print(msg);
            sys.exit(1)
        if args.args[0] in ['all', '*', '.']:
            args.slabel = args.slabel if args.slabel else "modbase"
            models = PDBMapModel.get_models()  # get all 2013 and 2016 ModBase models
            n = len(models)
            # If this is a parallel command with partition parameters
            if args.ppart != None and args.ppidx != None:
                psize = n / args.ppart  # floor
                if (args.ppart - 1) == args.ppidx:
                    models = models[args.ppidx * psize:]
                else:
                    models = models[args.ppidx * psize:(args.ppidx + 1) * psize]
                msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of ModBase\n" % (args.ppidx + 1, args.ppart)
                LOGGER.warning(msg)
                for i, row in enumerate(models):
                    print("## Processing (%s) %s (%d/%d) ##" % (
                        args.slabel, row['modelid'], i + (args.ppidx * psize) + 1, n))
                    pdbmap.load_model(row, label=args.slabel, io=None)
            else:
                msg = "WARNING (PDBMap) Uploading all %d ModBase 2013 and 2016 models.\n" % n
                LOGGER.warning(msg)
                for i, row in enumerate(models):
                    print("\n## Processing (%s) %s (%d/%d) ##" % (args.slabel, row['modelid'], i, n))
                    pdbmap.load_model(row, label=args.slabel, io=None)
        # Parse user-supplied homology models
        else:
            if not args.slabel:
                args.slabel = 'manual'
            if not len(args.args) > 1:
                msg = "ERROR (PDBMap): Must include ModBase-formatted summary file with load_model"
                raise Exception(msg)
            if "*" in args.args[1]:
                model_summary = [args.args[0]]
                models = glob.glob(args.args[1])
            else:
                model_summary = [args.args[0]]
                models = args.args[1:]
            for ms in model_summary:
                with open(ms, 'r') as fin:
                    fin.readline()  # burn the header
                    reader = csv.reader(fin, delimiter='\t')
                    n = len(models)
                    for i, row in enumerate(reader):
                        if not row[3].startswith("ENSP"):
                            i -= 1  # Not Ensembl. Decrement.
                            continue
                        row.append(models[i])
                        row.append(args.unp)
                        # If the summary file conforms to 2013 standard, reconcile with 2016
                        if len(row) != len(PDBMapModel._info_fields) - 2:
                            row.insert(1, None)
                            row.insert(1, None)
                        row = dict(list(zip(PDBMapModel._info_fields, row)))
                        print("## Processing (%s) %s (%d/%d) ##" % (args.slabel, row['modelid'], i, n))
                        # Convert the model summary row to a dictionary
                        # Load the ModBase model
                        pdbmap.load_model(row, label=args.slabel, io=None)

    ## load_data ##
    elif args.load_vcf:
        if args.dlabel:
            print('dlabel=%s' % args.dlabel)

    elif args.cmd == "load_data":
        if len(args.args) < 1:
            msg = "usage: pdbmap.py -c conf_file [--novep] load_data <data_file> <data_name> [data_file data_name] ...\n"
            msg += "alt:   pdbmap.py -c conf_file [--novep] --dlabel=<data_name> load_data <data_file> [data_file] ...\n"
            print(msg);
            sys.exit(1)
        pdbmap = PDBMap(vep=args.vep, vep_cache=args.vep_cache)
        # Process many data file(s) (set(s))
        if not args.dlabel:  # Assign individual labels
            dfiles = list(zip(args.args[0::2], args.args[1::2]))
        else:  # Assign shared label
            dfiles = list(zip(args.args, [args.dlabel for i in range(len(args.args))]))
        nrows = 0
        for dfile, dname in dfiles:
            print("## Processing (%s) %s ##" % (dname, dfile))
            nrows += pdbmap.load_data(dname, dfile, args.indexing, not args.novep, not args.noupload)
            print(" # %d data rows uploaded." % nrows)

    ## visualize ##
    elif args.cmd == "visualize":
        if len(args.args) < 3:
            msg = "usage: pdbmap.py -c conf_file visualize <entity> <data_name> <feature[,...]> <biounit[,...]> [minval:maxval,...] [color1,color2,...;...]\n"
            print(msg);
            sys.exit(1)
        pdbmap = PDBMap(idmapping=args.idmapping)
        entity = args.args[0]
        struct_label = 'pdb' if not args.slabel else args.slabel
        data_label = args.args[1]
        anno_list = args.args[2].split(',')
        if len(args.args) > 3 and args.args[3].lower() not in ['all', '.', ' ']:
            biounits = args.args[3].split(',')
        else:
            biounits = []
        spectrum_range = []
        if len(args.args) > 4:
            spectrum_range = [tuple([float(x) for x in p.split(':')]) for p in args.args[4].split(',')]
        colors = []
        if len(args.args) > 5:
            colors = [tuple([x for x in p.split(':')]) for p in args.args[5].split(',')]
        print("## Visualizing (%s) %s[%s]+%s.%s" % (struct_label,
                                                    entity, ','.join(biounits), data_label, ','.join(anno_list)),
              end=' ')
        if not colors:
            print(','.join(["(%0.2f..%0.2f)" % r for r in spectrum_range]))
        else:
            for i, r in enumerate(spectrum_range):
                for j, c in enumerate(range(int(r[0]), int(r[1]) + 1)):
                    print("%0.2f: %s;" % (c, colors[i][j]), end=' ')
                print('|')
            print('')
        pdbmap.visualize(entity, biounits, struct_label, data_label, anno_list, spectrum_range, colors)

    ## intersect ##
    elif args.cmd == "intersect":
        if not (args.slabel and args.dlabel):
            msg = "usage: pdbmap.py -c conf_file --slabel=<slabel> --dlabel=<data_name> intersect [quick]\n"
            print(msg);
            sys.exit(1)
        pdbmap = PDBMap()
        dname = args.dlabel
        if dname == 'all':
            dname = None
        slabel = args.slabel
        if slabel == 'all':
            slabel = None
        quick = True if len(args.args) > 0 and args.args[0].lower() in ['1', 'true', 'yes', 'quick', 'fast'] else False
        # nrows = QUICK_THRESH+1 if len(args.args) < 3 else int(args.args[2])
        if dname and slabel:
            print("## Intersecting %s with %s ##" % (dname, slabel))
        elif dname:
            print("## Intersecting %s with all structures/models ##" % dname)
        elif slabel:
            print("## Intersecting all genetic datasets with %s ##" % slabel)
        else:
            print("## Intersecting all genetic datasets with all structures/models ##")
        # quick = True if nrows < QUICK_THRESH else False
        print([" # (This may take a while) #", " # Using quick-intersect #"][int(quick)])
        nrows = pdbmap.intersect_data(dname, slabel, quick=quick)
        print(" # %d intersection rows uploaded." % nrows)

    ## no command specified ##
    else:
        msg = "PDBMap must be called with a command. Use -h for help.\n"
        LOGGER.error(msg)

    print('')
