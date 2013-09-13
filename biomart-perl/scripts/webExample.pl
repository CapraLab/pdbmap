# an example script demonstrating the use of BioMart webservice
use strict;
use LWP::UserAgent;


open (FH,$ARGV[0]) || die ("\nUsage: perl webExample.pl Query.xml\n\n");

my $xml;
while (<FH>){
    $xml .= $_;
}
close(FH);


my $path="http://www.biomart.org/biomart/martservice?";
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
my $ua = LWP::UserAgent->new;

my $response;

$ua->request($request, 
	     sub{   
		 my($data, $response) = @_;
		 if ($response->is_success) {
		     print "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);


