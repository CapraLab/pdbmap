#!/usr/bin/perl
# An example script demonstrating the use of BioMart API.
# This perl API representation is only available for configuration versions >=  0.5 
use strict;
use lib '/labs/twells/sivleyrm/pdbmap/biomart-perl/lib';
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "/labs/twells/sivleyrm/pdbmap/biomart.conf";
#
# NB: change action to 'clean' if you wish to start a fresh configuration  
# and to 'cached' if you want to skip configuration step on subsequent runs from the same registry
#

my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

		
	$query->setDataset("hsapiens_snp");
	$query->addFilter("variation_source", ["dbSNP"]);
	$query->addFilter("so_parent_name", ["coding_sequence_variant"]);
	$query->addAttribute("chr_name");
	$query->addAttribute("chrom_start");
	$query->addAttribute("refsnp_id");
	$query->addAttribute("refsnp_source");
	$query->addAttribute("ensembl_transcript_stable_id");
	$query->addAttribute("ensembl_gene_stable_id");
	$query->addAttribute("minor_allele_freq");

$query->formatter("TSV");

my $query_runner = BioMart::QueryRunner->new();
############################## GET COUNT ############################
# $query->count(1);
# $query_runner->execute($query);
# print $query_runner->getCount();
#####################################################################


############################## GET RESULTS ##########################
# to obtain unique rows only
# $query_runner->uniqueRowsOnly(1);

$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();
#####################################################################
