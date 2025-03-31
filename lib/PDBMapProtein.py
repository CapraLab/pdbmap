#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapProtein.py
# Author         : Chris Moth (2017-2022) and R. Michael Sivley before that
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2022-08-25
# Description    : Provides protein-related utilities, primarily external
#                : database ID matching to UniProt IDs. At this time, the
#                : class only provides services, and should not be
#                : instantiated.
# =============================================================================#

# See main check for cmd line parsing
import sys
import csv
import gzip
from typing import Dict, List, Union
import logging

LOGGER = logging.getLogger(__name__)


class PDBMapProtein:
    """
    Loads the tab-delimited HUMAN_9606_idmapping_sprot.dat.gz file from the filesystem and provides a set of
    class member functions (not object instance functions!) to cross-reference between different types of
    identifiers associated with proteins and transcripts.

    HUMAN_9606_idmapping_sprot.dat.gz is prepared by the DOWNLOAD_uniprot.bash script found in the
    pipeline_download_scripts/ repo at the GitHub site.
    """

    # _sec2prim is initialized from an entirely different file
    # uniprot_sec2prim_ac.txt which is a simple 2 column format
    # derived from Uniprot's source sec_ac.txt file

    # Each retired/secondary Uniprot AC maps to a single current Uniprot AC
    # These secondary Uniprot ACs never iso -n isoform extensions, and only
    # map to non isoform unps
    _sec2prim = {}  # Each retired/secondary Uniprot AC maps to a single current Uniprot AC

    # These dictionaries of dictionaries map between uniprot IDs in column 1
    # and IDs of type 'UD_type' in column 3.
    # For example, to find all the ENST000 Ensembl transcript IDs mapped to uniprot ID 'ABC123', you would look
    # at the list returned by:
    # _uniprot_to_ID['ABC123']['Ensembl_TRS']  To find the uniprot ID associated with an ENST000 Ensembl transcript
    # you would retrieve a matching list of uniprot IDs with:   _ID_to_uniprot['Ensembl_TRS']['ENST00000390362']
    # There are many type of cross-references in the ID mapping file
    # On initial load, we skip anywhere column 2 is not in this set here
    ID_TYPES_OF_INTEREST = {
        "UniProtKB-ID",
        "UniParc",
        "Gene_Name",  # One uniprot canonical ID may have multiple Gene_Name entries, but usually only 1
        "GeneID",
        "Ensembl_TRS",
        "Ensembl_PRO",
        "PDB",
        "RefSeq",
        "RefSeq_NT"
    }

    _uniprot_to_ID: Dict[str, Dict[str, list[str]]] = {}

    # Initialize the reverse lookups with empty entries for the ID_TYPES_OF_INTEREST
    _ID_to_uniprot: Dict[str, Dict[str, list[str]]] = {
        id_type: {} for id_type in ID_TYPES_OF_INTEREST
    }

    # I also need a dictionary to deal with user-supplied ENST transcript sequences
    # which lack version numbers .N appended.  In only a very  few cases
    # does the uniprot IDMapping file map an ENST123456 id to two IDs like
    # ENST123456,2  and ENST124356.3
    _unversioned_ensembl_ID_to_versioned_ensembl_ID: Dict[str, str] = {}

    # Refseq IDs can have .N suffixes, with the same challenges
    # as ensembl IDs

    _unversioned_refseq_ID_to_versioned_refseq_ID: Dict[str, str] = {}

    def __init__(self):
        msg = "ERROR (PDBMapProtein) This class should not be instantiated."
        raise Exception(msg)

    @classmethod
    def unp2uniparc(cls, uniprot_id: str):
        # Return UniParc associated with Uniprot AC
        # The pre-August-2022 versions would throw KeyError if uniprot_id was missing.
        # Here I just assert.
        assert uniprot_id in PDBMapProtein._uniprot_to_ID, \
            "Uniprot ID %s was not seen in the load of the idmapping file" % uniprot_id

        assert 'UniParc' in PDBMapProtein._uniprot_to_ID[uniprot_id], \
            "No UniParc cross reference was seen for uniprot id %s in the idmapping file " % uniprot_id

        return PDBMapProtein._uniprot_to_ID[uniprot_id]['UniParc']

    @classmethod
    def isCanonicalByUniparc(cls, uniprot_id: str):
        split_unp = uniprot_id.split('-')
        # unps without dashes are always canonical
        if len(split_unp) == 1:
            return True
        base_unp = split_unp[0]
        # If we lack a uniparc entry, then we can't claim canonical
        # if not PDBMapProtein.unp2uniparc(unp):
        #    return False

        return PDBMapProtein.unp2uniparc(uniprot_id) == PDBMapProtein.unp2uniparc(base_unp)

    @classmethod
    def all_unps(cls):
        for unp in PDBMapProtein._uniprot_to_ID.keys():
            yield unp

    @classmethod
    def best_unp(cls, uniprot_id: str) -> str:
        """
        Return an isoform-specific uniprot ID if one exists for a given uniprot ID

        Typical usage: Given uniprot ID 124356, return 124356-7 if that isoform specific, suffixed, ID has the
        same aminio acid sequence.  'samesness' is determined if both uniprot IDs have the same UniParc sequence IDs

        @param uniprot_id:  A uniprot ID.  Typically lacking an isoform specific -N suffix.
        @return: The isoform-specific uniprot_id with same Amino Acid sequence as the provided uniprot_id
        """

        # Caller has goofed if they sent us a None
        assert uniprot_id is not None

        if '-' in uniprot_id:  # If we already have a full '-' isoform specific uniprot ID, we cannot improve that
            return uniprot_id

        # If we can identify a canonical full unp isoform with dash return that
        uniparc_id_for_base_uniprot_id = PDBMapProtein.unp2uniparc(uniprot_id)[0]
        assert uniparc_id_for_base_uniprot_id

        # Given a uniparc ID, we now grab all the renated uniprot IDs
        # If there is an isoform-specific uniprot ID that matches the Uniparc ID of the uniprot_id
        # then _that_ one (with the dash) is the "best" one
        all_uniprot_ids_for_uniparc = PDBMapProtein._ID_to_uniprot['UniParc'].get(uniparc_id_for_base_uniprot_id, [])

        assert len(all_uniprot_ids_for_uniparc) > 0  # Total fail if the reverse cross-reference is missing somehow

        for possible_isoform_specific_uniprot_id in all_uniprot_ids_for_uniparc:
            if '-' in possible_isoform_specific_uniprot_id:  # Success: full isoform with dash for the base isoform
                return possible_isoform_specific_uniprot_id

        return uniprot_id

    @classmethod
    def best_unps(cls, uniprot_ids: List[str]) -> List[str]:
        return [PDBMapProtein.best_unp(uniprot_id) for uniprot_id in uniprot_ids]

    @classmethod
    def _refseq_id_2unp(cls, refseq_id: str, idmapping_refseq_key: str) -> List[str]:
        """

        @param refseq_id:
        @param idmapping_refseq_key:  Must be 'RefSeq_NT' or 'RefSeq_NP'
        @return:
        """

        # Sorry about nonsense cases - but dismiss them right away.
        if not refseq_id or refseq_id.startswith("NA"):
            return []

        # Make sure that the caller knows that transcript and protein keys are all that is allowed from
        # source column 2 for this query.
        assert idmapping_refseq_key in ['RefSeq_NT', 'RefSeq']
        # Return UniProt IDs associated with a refseq ID transcript

        # If this refseq_id_ has no version number, or a version unknown to the idmaping file, then try
        # to convert it to a current versioned refseq ID
        if refseq_id not in PDBMapProtein._ID_to_uniprot[idmapping_refseq_key]:
            unversioned_refseq_id = refseq_id.split('.')[0]
            if unversioned_refseq_id in PDBMapProtein._unversioned_refseq_ID_to_versioned_refseq_ID:
                refseq_id = PDBMapProtein._unversioned_refseq_ID_to_versioned_refseq_ID[unversioned_refseq_id]

        if refseq_id in PDBMapProtein._ID_to_uniprot[idmapping_refseq_key]:
            return PDBMapProtein.best_unps(
                PDBMapProtein._ID_to_uniprot[idmapping_refseq_key].get(refseq_id, []))

        # If the refseq_id_ has a .N version number at the end, then search for a versionless
        # RefSeq ID to hopefully match to uniprot ID(s)
        version_location = refseq_id.find('.')
        if version_location > -1:
            refseq_versionless = refseq_id[0:version_location]
            return PDBMapProtein.best_unps(
                PDBMapProtein._ID_to_uniprot[idmapping_refseq_key].get(refseq_versionless, []))

        return []  # No luck - no uniprot IDs match to this refseq ID

    @classmethod
    def refseqNT2unp(cls, refseq_transcript_id: str):
        """
        Return a list of uniprot IDs that cross-reference to a RefSeq ID in the IDMapping file.
        If the refseq ID has .NN version number and that is found in idmapping, then exact uniprot ID Xrefs will be 
        returned

        If no matches are found, then the version suffix is tripped off, and the cross references is searched fr again.  

        @param refseq_transcript_id: NM_1234... or XM_1234
        @return: 
        """
        return PDBMapProtein._refseq_id_2unp(refseq_transcript_id, idmapping_refseq_key='RefSeq_NT')

    @classmethod
    def refseqNP2unp(cls, refseq_protein_id: str):
        # Return UniProt IDs associated with RefSeqNP protein
        return PDBMapProtein._refseq_id_2unp(refseq_protein_id, idmapping_refseq_key='RefSeq')

    @classmethod
    def unp2refseqNT(cls, uniprot_id: str) -> List[str]:
        """
        Given a uniprot ID, return the entire list of refseq transcripts cross-referenced in the uniprot file
        We _try_ to use an isoform specific id first...because those refseq mapppings tend to be hgher quality
        sequence matches

        @param uniprot_id:
        @return: list of refseq IDs
        """
        best_refseq_ids = PDBMapProtein._uniprot_to_ID.get(
            PDBMapProtein.best_unp(uniprot_id), {}).get('RefSeq_NT', [])

        # If we can't get any refseqID(s) with the passeed-in uniprot ID.
        # But the isoform specific uniprot ID passed in is  cananonical
        # Then, grab a refseq ID from the canonical uniprot ID (no '-')
        # This is likely wrong - but at least you get something...
        # The uniprot refseq mappings are very strange...
        if not best_refseq_ids and PDBMapProtein.isCanonicalByUniparc(uniprot_id):
            best_refseq_ids = PDBMapProtein._uniprot_to_ID.get(
                uniprot_id[0:6], {}).get('RefSeq_NT', [])

        # We want an XM_ to sort to the end if there are multiples
        # This could easily be an empty list we are returning.
        best_refseq_ids.sort()
        return best_refseq_ids

    @classmethod
    def unp2refseqNP(cls, uniprot_id: str) -> List[str]:
        """
        Given a uniprot ID, return the entire list of refseq proteins cross-referenced in the uniprot file
        We _try_ to use an isoform specific id first...because those refseq mapppings tend to be hgher quality
        sequence matches

        @param uniprot_id:
        @return: list of refseq IDs
        """
        best_refseq_ids = PDBMapProtein._uniprot_to_ID.get(
            PDBMapProtein.best_unp(uniprot_id), {}).get('RefSeq', [])

        # We want an XM_ to sort to the end if there are multiples
        # This could easily be an empty list we are returning.
        best_refseq_ids.sort()
        return best_refseq_ids

    @classmethod
    def unp2uniprotKB(cls, uniprot_id: str):
        # Return  the UNIPROTKB entry for this protein RefSeq protein associatied with UniProt ID
        uniprot_knowledgebase_id_list = PDBMapProtein._uniprot_to_ID.get(uniprot_id, {}).get('UniProtKB-ID', [])
        return uniprot_knowledgebase_id_list[0] if uniprot_knowledgebase_id_list else ''

    @classmethod
    def unp2enst(cls, uniprot_id: str):
        # Return Ensembl Transcript IDs associated with UniProt ID
        ensembl_transcripts = PDBMapProtein._uniprot_to_ID.get(uniprot_id, {}).get('Ensembl_TRS', [])

        # If - in unp name it is isform specific, If also canonical, append  canonical ENST Transcript ID from base
        # uniprot ID
        additional_uniprot_id_to_mine = None

        if '-' in uniprot_id:
            if PDBMapProtein.isCanonicalByUniparc(uniprot_id):
                additional_uniprot_id_to_mine = uniprot_id.split('-')[0]
        else:  # If a base unp is given, then see add transcripts associated wtih canonical unp
            unp_canonical_isoform = PDBMapProtein.best_unp(uniprot_id)
            if '-' in unp_canonical_isoform:
                additional_uniprot_id_to_mine = unp_canonical_isoform

        if additional_uniprot_id_to_mine:
            additional_transcripts = PDBMapProtein._uniprot_to_ID.get(
                additional_uniprot_id_to_mine, {}).get('Ensembl_TRS', [])
            for additional_transcript in additional_transcripts:
                if additional_transcript not in ensembl_transcripts:
                    ensembl_transcripts.append(additional_transcript)

        # if ensembl_transcripts:
        return ensembl_transcripts

        # else Could this unp be secondary somehow?
        # 2022-August-29  Chris Moth removed all secondary uniprot untangling
        # unpbase = unp.split('-')[0]
        # if unpbase in PDBMapProtein._sec2prim:
        #    return PDBMapProtein._unp2enst.get(PDBMapProtein._sec2prim[unpbase], [])
        #
        # return []

    @classmethod
    def unp2hgnc(cls, uniprot_id: str) -> Union[str, None]:
        # Return HGNC gene name associated with UniProt ID
        # For this lookup, always truncate off isoform information
        uniprot_id_base = uniprot_id.split('-')[0]
        gene_name_list = PDBMapProtein._uniprot_to_ID.get(uniprot_id_base, {}).get('Gene_Name', [])
        return gene_name_list[0] if gene_name_list else None
        # Removing all secondary references
        # if unpbase in PDBMapProtein._sec2prim:
        #     if PDBMapProtein._sec2prim[unpbase] in PDBMapProtein._unp2hgnc:
        #        return PDBMapProtein._unp2hgnc[PDBMapProtein._sec2prim[unpbase]]
        # LOGGER.warning("unp of %s not found in PDBMapProtein._unp2hgnc or _sec2prim\n" % unpbase)
        # return None

    @classmethod
    def hgnc2unp(cls, gene_name: str) -> Union[List[str], None]:
        """
        @param gene_name: Example KCNC1
        @return: a corresponding List[] of uniprot identiferss, usually with only one element,
                 with isoform appended if possible.  Example for PCNC1: ['P48547-1']
        """
        # Uniprot identifier for a HGNC gene name
        #
        uniprot_id_list = PDBMapProtein._ID_to_uniprot['Gene_Name'].get(gene_name, [])

        return PDBMapProtein.best_unp(uniprot_id_list[0]) if uniprot_id_list else None

    @classmethod
    def unp2gene_id(cls, uniprot_id: str) -> int:
        # Return the gene ID number from the IDMapping File
        # This is useful for forming direct-access gene
        # information URLs
        # For this lookup, always truncate off isoform information
        uniprot_id_base = uniprot_id.split('-')[0]
        gene_id_list = PDBMapProtein._uniprot_to_ID.get(uniprot_id_base, {}).get('GeneID', [])
        if gene_id_list:
            return int(gene_id_list[0])
        LOGGER.warning("unp of %s not found in PDBMapProtein._unp2gene_id\n" % uniprot_id_base)
        return 0  # No gene ID found

    @classmethod
    def gene2gene_id(cls, gene_name: str) -> Union[int, None]:
        # Return the gene ID number from the IDMapping File
        # This is useful for forming direct-access gene
        # information URLs
        # For this lookup, always truncate off isoform information
        uniprot_id = PDBMapProtein._ID_to_uniprot.get('Gene_Name', {}).get(gene_name, [])
        if uniprot_id:
            uniprot_id_base = uniprot_id[0].split('-')[0]
            return PDBMapProtein.unp2gene_id(uniprot_id_base)
        else:
            LOGGER.warning("gene_name %s not found in IDMapping file", gene_name)
        return None

    @classmethod
    def unp2ensp(cls, uniprot_id: str):
        # Return Ensembl Transcript IDs associated with UniProt ID
        ensembl_proteins = PDBMapProtein._uniprot_to_ID.get(uniprot_id, {}).get('Ensembl_PRO', [])

        # If - in unp name it is isform specific, If also canonical, append  canonical ENST Transcript ID from base
        # uniprot ID
        additional_uniprot_id_to_mine = None
        if '-' in uniprot_id:
            if PDBMapProtein.isCanonicalByUniparc(uniprot_id):
                additional_uniprot_id_to_mine = uniprot_id.split('-')[0]
        else:  # If a base unp is given, then see add proteins associated wtih canonical unp
            unp_canonical_isoform = PDBMapProtein.best_unp(uniprot_id)
            if '-' in unp_canonical_isoform:
                additional_uniprot_id_to_mine = unp_canonical_isoform

        if additional_uniprot_id_to_mine:
            additional_proteins = PDBMapProtein._uniprot_to_ID.get(
                additional_uniprot_id_to_mine, {}).get('Ensembl_PRO', [])
            for additional_protein in additional_proteins:
                if additional_protein not in ensembl_proteins:
                    ensembl_proteins.append(additional_protein)

        # if ensembl_proteins:
        return ensembl_proteins

    @classmethod
    def versioned_ensembl_transcript(cls, ensembl_transcript_or_protein_id: str) -> str:
        """
        Given an ensembl ID without a version suffix, such as:
        ENST00000646695
        _IF_ this ID is mentioned in the uniprot IDmapping file, then return the same ID
        but with the version appended as in ENST00000646695.2

        @param ensembl_transcript_or_protein_id: An ensembl transcript or protein identifier
        @return: the same identifier, but with a suffix added if known to the uniprot idmapping file
        """

        if ensembl_transcript_or_protein_id in PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID:
            return PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID[ensembl_transcript_or_protein_id]
        else:
            return ensembl_transcript_or_protein_id

    @classmethod
    def ensp2unp(cls, ensembl_transcript_id: str) -> List[str]:
        """
        Return the uniprot ID(s) cross-referenced to an ensembl transcript ID
        Versionless ensembl transcript IDs (lacking a .N suffix) will be converted
        to versioned.

        @param ensembl_transcript_id: ENST0000124356 etc or with .N suffix
        @return: The cross-referenced uniprot IDs, isoform-specific when available
                 It should never happen that the list has more than one element
                 But often, an empty list is returned
        """
        if '.' not in ensembl_transcript_id:
            # We must try to "up convert" the versionless ensembl ID to a versioned ID
            # because only versioned Ensembl IDs are in the uniprot IDMapping files
            # Yet we encountered non-versioned ENST0000124356 identifiers often in the pipeline
            # HOWEVER, there are a very few cases where non-versioned ENST transcript IDs are
            # in the Idmapping file
            if ensembl_transcript_id in PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID:
                ensembl_transcript_id = PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID[
                    ensembl_transcript_id
                ]

            # At this point, ensembl_transcript_id is most likely "versioned" but that does not mean
            # we have a cross-reference for it in the IdMapping file
            return PDBMapProtein.best_unps(
                PDBMapProtein._ID_to_uniprot['Ensembl_TRS'].get(ensembl_transcript_id, []))

    @classmethod
    def enst2unp(cls, ensembl_protein_id: str) -> List[str]:
        """
        Return the uniprot ID(s) cross-referenced to an ensembl protein ID
        Versionless ensembl protein IDs (lacking a .N suffix) will be converted
        to versioned.

        @param ensembl_protein_id: ENST0000124356 etc or with .N suffix
        @return: The cross-referenced uniprot IDs, isoform-specific when available
                 It should never happen that the list has more than one element
                 But often, an empty list is returned
        """
        if '.' not in ensembl_protein_id:
            # We must try to "up convert" the versionless ensembl ID to a versioned ID
            # because only versioned Ensembl IDs are in the uniprot IDMapping files
            # Yet we encountered non-versioned ENST0000124356 identifiers often in the pipeline
            # HOWEVER, there are a very few cases where non-versioned ENST protein IDs are
            # in the Idmapping file
            if ensembl_protein_id in PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID:
                ensembl_protein_id = PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID[
                    ensembl_protein_id
                ]

        # At this point, ensembl_protein_id is most likely "versioned" but that does not mean
        # we have a cross-reference for it in
        return PDBMapProtein.best_unps(
            PDBMapProtein._ID_to_uniprot['Ensembl_TRS'].get(ensembl_protein_id, []))

    @classmethod
    def unp2pdb(cls, uniprot_id: str):
        """
        Return Protein Data Bank ID associated with UniProt ID
        PDBs never have transcript identifiers in their depositions
        So we tend to find PDB lists from the canonical uniprot ID.

        @param uniprot_id:
        @return:
        """

        pdb_list = PDBMapProtein._uniprot_to_ID.get(uniprot_id, {}).get('PDB', [])
        if not pdb_list:
            base_unp = uniprot_id.split('-')[0]
            pdb_list = PDBMapProtein._uniprot_to_ID.get(base_unp, {}).get('PDB', [])

        # Chris Moth notes that on 2022-August the sec2prim seems hopelessly broken
        # We need to look at it again - but uniprot seems to not keep t going.
        # For PDB list, we probably should go directly to SIFTS API.
        # If still no joy - try to get a pdb based on a retired UniProt ID
        # (RCSB does NOT update .pdbs as uniprot IDentifiers change)
        # if not pdbs:
        #    currentUnp = PDBMapProtein._sec2prim.get(base_unp, None)
        #    if currentUnp:
        #        pdbs = PDBMapProtein._unp2pdb[currentUnp]

        return pdb_list

    @classmethod
    def isunp(cls, unp):
        # Return True if the provided ID is a UNP ID
        return unp in PDBMapProtein._uniprot_to_ID
        #
        # 2022-08-30 - Chris Moth, again, removing all these secondary uniprot ID lookups
        # or unp in PDBMapProtein._sec2prim

    @classmethod
    def ishgnc(cls, gene_name):
        # Return True if the provided ID is a known gene name
        return bool(PDBMapProtein._ID_to_uniprot.get('Gene_Name', {}).get(gene_name, []))

    @classmethod
    def load_idmapping(cls, id_mapping_file_name: str) -> None:
        # We check that ONLY our special HUMAN uniprot idmapping file
        # is included.  This is the special, swiss-curated-only file
        # that we made from the DOWNLOAD_uniprot.bash script
        #
        human_id_mapping_swiss_curated_filename = "HUMAN_9606_idmapping_sprot.dat.gz"
        if human_id_mapping_swiss_curated_filename not in id_mapping_file_name:
            LOGGER.critical("id_mapping_file_name %s does not include: %s",
                            id_mapping_file_name, human_id_mapping_swiss_curated_filename)
            sys.exit()

        LOGGER.info("Opening idmapping file: %s", id_mapping_file_name)
        # Load UniProt crossreferences, keyed on UniProt
        # I adopt the column header names
        # 'unp', 'ID_type', and 'ID' because these are used in
        _unp: str
        _id_type: str
        _id: str
        # our SQL load of the IDMapping file triples
        with gzip.open(id_mapping_file_name, 'rt', encoding='ascii') as fin:
            reader = csv.reader(fin, delimiter='\t')
            # Parse the comprehensive UniProt ID mapping resource
            # and extract all necessary information from the 3 columns
            for _unp, _id_type, _id in reader:
                """
                # Extract the UniProt isoform identifier (unpN), if present
                unpN = None
                if '-' in unp:
                    unpN = unp  # Retain the full UniprotAC-nn as the isoform
                    unpBest = unp
                    unp, iso = unp.split('-')
                else:
                    # Proteins with a single (thus canonical) isoform have no identifier
                    unpBest = unp
                """
                # This is necessary to avoid an unnecessary GRCh37/GRCh38 version mismatch
                # e.g. Some canonical isoforms in UniProt map to transcripts that only exist
                # in GRCh38. However, -002+ isoforms sometimes match previous transcripts
                # from GRCh37, so they need to be retained for completeness.
                # The best fix would be to update all genetic resources and datasets to GRCh38, but
                # so long as our collaborators and all major genetics cohorts continue to provide
                # all data in GRCh37, that is not a feasible solution.
                # # Only record EnsEMBL transcripts/proteins matching the UniProt canonical isoform
                # if iso != "1":
                #   continue

                # RESUME WORKING HERE
                if _id_type in PDBMapProtein.ID_TYPES_OF_INTEREST:
                    if _unp not in PDBMapProtein._uniprot_to_ID:
                        PDBMapProtein._uniprot_to_ID[_unp] = {}
                    _ids_list = PDBMapProtein._uniprot_to_ID[_unp].get(_id_type, [])
                    _ids_list.append(_id)
                    PDBMapProtein._uniprot_to_ID[_unp][_id_type] = _ids_list

                    _uniprot_id_list = PDBMapProtein._ID_to_uniprot[_id_type].get(_id, [])
                    _uniprot_id_list.append(_unp)
                    PDBMapProtein._ID_to_uniprot[_id_type][_id] = _uniprot_id_list

                    # In the subcase where the third column is a versioned ENSEMBL ID, we need to
                    # additionally create a map with the unversioned ID
                    # Note we are still in a block retaining only ID_TYPES_OF_INTEREST
                    if _id_type.startswith('Ensembl'):
                        _unversioned_ensembl_ID = _id.split('.')[0]  # Strip off the trailing version number
                        _existing_versioned_ensembl_ID = PDBMapProtein.\
                            _unversioned_ensembl_ID_to_versioned_ensembl_ID.get(
                                _unversioned_ensembl_ID, None)

                        # if we have already created a cross-ref from ENST1234 to ENST1234.9, then this is a second
                        # and that merits a warning about the idmapping file
                        if not _existing_versioned_ensembl_ID:
                            PDBMapProtein._unversioned_ensembl_ID_to_versioned_ensembl_ID[
                                _unversioned_ensembl_ID] = _id
                        else:  # We need to crash if we see 2 different versions of ID mapping to our base ENST
                            if _existing_versioned_ensembl_ID != _id:
                                LOGGER.critical("Ensembl ID %s is not uniquely cross-referenced in uniprot",
                                                _unversioned_ensembl_ID)
                                sys.exit("Parsing halted due to ambiguous Ensembl versionless->Versioned xref")

                    # In the subcase where the third column is a versioned ENSEMBL ID, we need to
                    # additionally create a map with the unversioned ID
                    # Note we are still in a block retaining only ID_TYPES_OF_INTEREST
                    if _id_type.startswith('RefSeq'):
                        _unversioned_refseq_ID = _id.split('.')[0]  # Strip off the trailing version number
                        _existing_versioned_refseq_ID = PDBMapProtein._unversioned_refseq_ID_to_versioned_refseq_ID.get(
                            _unversioned_refseq_ID, None)

                        # if we have already created a cross-ref from ENST1234 to ENST1234.9, then this is a second
                        # and that merits a warning about the idmapping file
                        if not _existing_versioned_refseq_ID:
                            PDBMapProtein._unversioned_refseq_ID_to_versioned_refseq_ID[
                                _unversioned_refseq_ID] = _id
                        else:  # We need to crash if we see 2 different versions of ID mapping to our base RefSeq
                            if _existing_versioned_refseq_ID != _id:
                                LOGGER.critical("RefSeq ID %s is not uniquely cross-referenced in uniprot",
                                                _unversioned_refseq_ID)
                                _id_version = int(_id.split('.')[1])
                                _existing_version = int(_existing_versioned_refseq_ID.split('.')[1])
                                # If we encounter NP_057260.3 after NO_057260.2, then install the larger one
                                if _id_version > _existing_version:
                                    PDBMapProtein._unversioned_refseq_ID_to_versioned_refseq_ID[
                                        _unversioned_refseq_ID] = _id
                                # sys.exit("Parsing halted due to ambiguous RefSeq versionless->Versioned xref")

        LOGGER.info("Cross-references for %d uniprot IDs loaded into RAM" % len(PDBMapProtein._uniprot_to_ID))

    @classmethod
    def load_sec2prim(cls, sec2prim_fname: str):
        # Method to load the UniProt secondary -> primary AC mapping
        # load_idmapping MUST be called before calling this function
        # I believe this secondary protein ID database is not longe rused anywhere.
        if not PDBMapProtein._uniprot_to_ID:
            raise Exception("ERROR (PDBMapProtein) load_sec2prim cannot be called before load_idmapping")
        with open(sec2prim_fname) as fin:
            reader = csv.reader(fin, delimiter='\t')
            sec2prim = {}
            for (sec, prim) in reader:
                # Ensure that this primary AC is human and mappable
                if prim in PDBMapProtein._uniprot_to_ID:
                    sec2prim[sec] = prim
        PDBMapProtein._sec2prim = sec2prim

    @classmethod
    def load_sprot(cls, sprot_fname):
        # Method to load all SwissProt IDs
        # I do not believe this is used anywhere in the pipeline
        with open(sprot_fname) as fin:
            sprot = []
            unp2species = {}
            for line in fin:
                # Extract the field ID, always [0:2]
                # Trim to AC list (AC field assumed)
                row = [line[0:2], line.rstrip()[2:-1].lstrip()]
                if row[0] == 'ID':
                    species = row[1].split()[0].split('_')[-1]
                elif row[0] == 'AC':
                    ac_list = row[1].split('; ')
                    sprot.append(ac_list[0])  # Most up to date AC
                    for ac in ac_list:
                        unp2species[ac] = species
        sprot.sort()
        PDBMapProtein.sprot = sprot
        PDBMapProtein.unp2species = unp2species

    @classmethod
    def check_loaded(cls):
        # Checks if external database ID mapping has been loaded
        if not PDBMapProtein._uniprot_to_ID or not PDBMapProtein._ID_to_uniprot:
            msg = "ERROR (UniProt) ID Mapping must be loaded before use."
            raise Exception(msg)
        # Checks if secondary to primary UniProt ID mapping has been loaded
        # Probably should remove this one.
        return bool(PDBMapProtein._sec2prim)


# Main check
if __name__ == "__main__":
    sys.stderr.write("Class definition. Should not be called from command line.\n")
    sys.exit(1)
