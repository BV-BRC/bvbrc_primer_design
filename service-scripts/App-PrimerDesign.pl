#
# The PrimerDesign application.
use strict;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Slurp;
use File::Basename;
use IPC::Run 'run';
use JSON;
use File::Copy ('copy', 'move');
use P3DataAPI;
use Bio::KBase::AppService::AppConfig;
use Bio::KBase::AppService::AppScript;
use Cwd;
use URI::Escape;

our $global_ws;
our $global_token;

our $shock_cutoff = 10_000;

our $debug = 0;
$debug = $ENV{"PrimerDesignDebug"} if exists $ENV{"PrimerDesignDebug"};
print STDERR "debug = $debug\n" if $debug;
print STDERR "args = ", join("\n", @ARGV), "\n" if $debug;
our @analysis_step => ();# collect info on sequence of analysis steps
our @step_stack => (); # for nesting of child steps within parent steps

my $data_url = Bio::KBase::AppService::AppConfig->data_api_url;
#my $data_url = "http://www.alpha.patricbrc.org/api";
print STDERR "data_url=\n$data_url\n" if $debug;

my $script = Bio::KBase::AppService::AppScript->new(\&design_primers, \&preflight);
my $rc = $script->run(\@ARGV);

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;
    print STDERR "preflight: num params=", scalar keys %$params, "\n";

    my $time = 60 * 60 * 1;

    my $pf = {
	cpu => 1,
	memory => "16G",
    runtime => $time,
	storage => 0,
	is_control_task => 0,
    };
    return $pf;
}

sub design_primers {
    my ($app, $app_def, $raw_params, $params) = @_;
    my %primer3_params; # copy the params that need to be passed to primer3_core
    for my $name (keys %{$params}) {
        next if $name =~ /^[a-z]/; # parameters with lower-case keys are not meant for primer3_core
        $primer3_params{$name} = $params->{$name};
    }

    my $time1 = `date`;
    print "Proc DesignPrimers ", Dumper($app_def, $raw_params, $params);
    $global_token = $app->token()->token();
    # print STDERR "Global token = $global_token\n";
    my @outputs; # array of tuples of (filename, filetype)
    my $api = P3DataAPI->new;
    my $tmpdir = File::Temp->newdir( "/tmp/PrimerDesign_XXXXX", CLEANUP => !$debug );
    system("chmod", "755", "$tmpdir");
    print STDERR "created temp dir: $tmpdir, cleanup = ", !$debug, "\n";
   
    run("echo $tmpdir && ls -ltr $tmpdir");

    if ($params->{input_type} eq "sequence_text") {
        $primer3_params{SEQUENCE_TEMPLATE} = $params->{sequence_input};
    }
    elsif ($params->{input_type} eq "workspace_fasta") {
        my $fasta_input_file = $params->{sequence_input};
        $fasta_input_file =~ s/.*\///; #should yield basename
        my @cmd = ("p3-cp", "ws:" . $params->{sequence_input}, $fasta_input_file);
        print STDERR "@cmd\n";
        my $ok = IPC::Run::run(\@cmd);
        if (!$ok)
        {
            warn "Error $? copying output with @cmd\n";
        }
        open F, $fasta_input_file;
        my $seq = '';
        $_ = <F>;
        $_ =~ />(\S+)/ or die "$fasta_input_file does not look like fasta";
        $primer3_params{SEQUENCE_ID} = $1;
        print STDERR "Reading sequence from workspace file $params->{workspace_fasta}\n";
        print STDERR "setting SEQUENCE_ID to $primer3_params{SEQUENCE_ID} from workspace file\n";
        print STDERR "(was $params->{SEQUENCE_ID})\n" if $params->{SEQUENCE_ID} and $params->{SEQUENCE_ID} ne $primer3_params{SEQUENCE_ID};
        while (<F>) {
            last if /^>/;
            chomp;
            $seq .= $_;
        }
        $primer3_params{SEQUENCE_TEMPLATE} = $seq;
        close F;
    }
    elsif ($params->{input_type} eq "database_id") {
        print STDERR "Input type is 'database_id' which is not supported yet.";
        exit(1);
    }    
    else {
        print STDERR "Input type is '$params->{input_type}' which is not supported yet.";
        exit(1);
    }    

    if ($primer3_params{SEQUENCE_TEMPLATE} =~ /[[\]<>{}-]/) {
        print STDERR "Need to handle sequence markups\n" if $debug;
        handle_sequence_region_markup(\%primer3_params);
    }

    # write primer3 input to file
    my $p3params_file = "$params->{output_file}_Primer3_input.txt";
    open F, ">$tmpdir/$p3params_file";
    print STDERR "primer3_params = %primer3_params\n";
    for my $name (keys %primer3_params) {
        print F "$name=$primer3_params{$name}\n";
    }
    print F "=\n"; #marks end of input to primer3
    close F;

    my $cwd = getcwd();
    chdir($tmpdir);

    my $primer3_output_file = "$params->{output_file}_Primer3_output.txt";
    my @command = ("primer3_core", "--output", $primer3_output_file);
    print "run command: @command\n";
    my $out = run(\@command, "<", $p3params_file);
    #my ($out, $err) = run_cmd(", @command));
    #print STDERR "\nSTDOUT:\n$out\n";
    #print STDERR "\nSTDERR:\n$err\n";
    my %resultsHash;
    print STDERR "Now showing results from runPrimer3\n";
    open F, $primer3_output_file;
    while (<F>)
    {
        chomp;
        my ($key, $val) = split("=", $_, '2');
        $resultsHash{$key} = $val;
        print "results: $key\tvalue=$resultsHash{$key}\n" if $debug;
    }
    my $pair_count = $resultsHash{PRIMER_PAIR_NUM_RETURNED}; 

    if (0) {
        my $html .= generate_dynamic_html(\%resultsHash);
        my $html_file = "$tmpdir/$params->{output_file}_dynamic_report.html";
        open F, ">$html_file";
        print F $html;
        close F;
        push @outputs, [$html_file, "html"];
    }

    my $html = format_primer3_output_to_single_html_table(\%resultsHash);
    my $html_file = "$tmpdir/$params->{output_file}_table.html";
    open F, ">$html_file";
    print F $html;
    close F;
    push @outputs, [$html_file, "html"];

    if ($resultsHash{"PRIMER_PAIR_NUM_RETURNED"} > 0) {
        my $fasta = write_primers_to_fasta(\%resultsHash);
        my $fasta_file = "$tmpdir/$params->{output_file}_primers.fasta";
        open F, ">$fasta_file";
        print F $fasta;
        close F;
        push @outputs, [$fasta_file, "Feature_DNA_FASTA"];
    }

    push @outputs, [$primer3_output_file, "txt"];
    push @outputs, [$p3params_file, "txt"];

    print STDERR '\@outputs = '. Dumper(\@outputs);
    my $output_folder = $app->result_folder();
    for my $output (@outputs) {
        my($ofile, $type) = @$output;
        next if $type eq 'folder';
        
        if (! -f $ofile) {
            warn "Output file '$ofile' of type '$type' does not exist\n";
            next;
        }
        
        my $filename = basename($ofile);
        #print STDERR "Output folder = $output_folder\n";
        print STDERR "Saving $filename => $output_folder as $type\n";
        if (0) {
           $app->workspace->save_file_to_file($ofile, {}, "$output_folder/$filename", $type, 1,
               (-s $ofile > $shock_cutoff ? 1 : 0), # use shock for larger files
               $global_token);
        }
        else {
            my $ext = $1 if $ofile =~ /.*\.(\S+)$/;
            my @cmd = ("p3-cp", "-f", "-m", "${ext}=$type", $ofile, "ws:" . $app->result_folder);
            print STDERR "@cmd\n";
            my $ok = IPC::Run::run(\@cmd);
            if (!$ok)
            {
                warn "Error $? copying output with @cmd\n";
            }
        }
    }
    chdir($cwd);
    my $time2 = `date`;
    #print STDERR ("Start: $time1"."End:   $time2"r. "$tmpdir/DONE\n");
}

sub write_primers_to_fasta {
    my $results = shift;
    my $num_pairs = $results->{"PRIMER_PAIR_NUM_RETURNED"};
    my $retval = "";
    for my $i (0 .. $num_pairs-1) {
        my $display_index = $i+1;
        $retval .= ">Primer_" . $display_index . "_forward\n" . $results->{"PRIMER_LEFT_${i}_SEQUENCE"} . "\n";
        $retval .= ">Primer_" . $display_index . "_reverse\n" . $results->{"PRIMER_RIGHT_${i}_SEQUENCE"} . "\n";
        $retval .= ">Internal_oligo_" . $display_index . "\n" . $results->{"PRIMER_INTERNAL_${i}_SEQUENCE"} . "\n";
    }
    return $retval;
}

sub format_primer3_output_to_single_html_table {
    my $results = shift;
    print "format to table: results = $results\n";
    my $html = "";
    # Figure out if any primers were returned and
    # write a helping page if no primers are returned  
    unless ($results->{"PRIMER_PAIR_NUM_RETURNED"}) {
        $html = "No primers were returned.\n";
        return $html
    }
    $html .= qq(<head>
<style>
h3  {
    font-size: 1.2em;
    margin-bottom: .2em;
}
table, th, td {
  border: 1px solid black;
}
td.empty {
    background: lightgray;
    }

td.sep {
	background: gray;
	height: 0.5em;
	}

table {
    font-family:sans-serif; 
    margin: .2em;
    border-collapse: collapse;
}

td.num {
    text-align:right;
}

th {
    padding: .1em;
    text-align: left;
}

td {
    padding-left: .5em;
    padding-right: .5em;
}
</style>
</head>
<body>
);
    $html .= "Sequence ID: $results->{SEQUENCE_ID}<br>\n" if exists $results->{SEQUENCE_ID};
    $html .= "Sequence length: " . length($results->{SEQUENCE_TEMPLATE}) . "<br>\n";
    if (exists $results->{TARGET_REGION}) {
        $html .= "Target region:";
        for my $pair (@{$results->{TARGET_REGION}}) { 
            $html .= " " . join(",", @$pair);
        }
        $html .= "<br>\n";
    }  
    if (exists $results->{INCLUDED_REGION}) {
        $html .= "Included region:";
        for my $pair (@{$results->{INCLUDED_REGION}}) { 
            $html .= " " . join(",", @$pair);
            }
        $html .= "<br>\n";
    }  
    if (0 and exists $results->{PRIMER_PRODUCT_SIZE_RANGE}) {
        $html .= "Product size range: " . $results->{PRIMER_PRODUCT_SIZE_RANGE};
        $html .= "<br>\n";
    }  
    if (0 and exists $results->{PRIMER_MAX_SIZE}) {
        $html .= "Primer max length: " . $results->{PRIMER_MAX_SIZE};
        $html .= "<br>\n";
    }  
    if (0 and exists $results->{PRIMER_MIN_SIZE}) {
        $html .= "Primer min length: " . $results->{PRIMER_MIN_SIZE};
        $html .= "<br>\n";
    } 
    $html .= "<small>Note: sequence positions and pair indexes are 1-based here, but 0-based in Primer3 output.</small><br>\n"; 
    $html .= "<table>\n";
    $html .= "<tr><th>Pair#</th>\n";
    #<th>Region</th><th>Sequence (5'->3')</th><th class='primer'>Start</th><th class='primer'>End</th><th class='primer'>Length</th><th class='primer'>Tm</th><th class='primer'>GC%</th><th class='primer'>Compl. <br>any</th><th class='primer'>Compl. <br>end</th></tr>\n";
    $html .= "<th title=\"Class of region.\">Region</th>\n";
    $html .= "<th title=\"Sequence of the oligomer.\">Sequence (5&prime;&rarr;3&prime;)</th>\n";
    $html .= "<th title=\"Start position in target sequence. (1-based)\">Start</th>\n";
    $html .= "<th title=\"End position in target sequence. (1-based)\">End</th>\n";
    $html .= "<th title=\"Length in bases.\">Len</th>\n";
    $html .= "<th title=\"The melting temperature for the selected oligomer.\">&nbsp;T<sub>m</sub> &nbsp;</th>\n";
    $html .= "<th title=\"Percent GC for the selected oligo.\">GC%</th>\n";
    $html .= "<th title=\"Tendency of a primer to bind to itself, interfering with target sequence binding.\">Any <br>Compl.</th>\n";
    $html .= "<th title=\"The tendency of the 3'-END to bind an identical primer, allowing it to form a primer-dimer.\">End <br>Compl.</th>";
    $html .= "</tr>\n";
    my $pair_count = $results->{"PRIMER_PAIR_NUM_RETURNED"};
    for my $index (0..($pair_count-1)) {
	    $html .= "<tr>";
	    my ($start, $length) = split(",", $results->{"PRIMER_LEFT_${index}"});
	    my $end = $start + $length;
        $start += 1; # translate to 1-based coordinates for user
        my $product_start = $start;
	    my $num_rows = 3;
	    $num_rows++ if exists $results->{"PRIMER_INTERNAL_$index"};
	    $html .= "<td rowspan='$num_rows'>" . ($index + 1) . "</td>\t<td class='primer'>Forward primer</td>\t";
            $html .= "<td>" . $results->{"PRIMER_LEFT_${index}_SEQUENCE"} . "</td>\t";
            $html .= "<td class='num'>$start</td>\t";
            $html .= "<td class='num'>$end</td>\t";
            $html .= "<td class='num'>$length</td>\t";
            $html .= "<td class='num'>" . sprintf("%.2f", $results->{"PRIMER_LEFT_${index}_TM"}) . "</td>\t";
            $html .= "<td class='num'>" . sprintf("%d", $results->{"PRIMER_LEFT_${index}_GC_PERCENT"}) . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_LEFT_${index}_SELF_ANY_TH"} . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_LEFT_${index}_SELF_END_TH"} . "</td>\n";
	    $html .= "</tr>\n";

	    ($start, $length) = split(",", $results->{"PRIMER_RIGHT_${index}"});
        $start += 1; # translate to 1-based coordinates for user
        my $product_end = $start;
	    my $end = $start - $length + 1;
	    $html .= "<td>Reverse primer</td>\t";
            $html .= "<td>" . $results->{"PRIMER_RIGHT_${index}_SEQUENCE"} . "</td>\t";
            $html .= "<td class='num'>$start</td>\t";
            $html .= "<td class='num'>$end</td>\t";
            $html .= "<td class='num'>$length</td>\t";
            $html .= "<td class='num'>" . sprintf("%.2f", $results->{"PRIMER_RIGHT_${index}_TM"}) . "</td>\t";
            $html .= "<td class='num'>" . sprintf("%d", $results->{"PRIMER_RIGHT_${index}_GC_PERCENT"}) . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_RIGHT_${index}_SELF_ANY_TH"} . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_RIGHT_${index}_SELF_END_TH"} . "</td>\n";
	    $html .= "</tr>\n";

	    if (exists $results->{"PRIMER_INTERNAL_${index}"}) {
		($start, $length) = split(",", $results->{"PRIMER_INTERNAL_${index}"});
            my $end = $start + $length;
            $start += 1; # translate to 1-based coordinates for user
            $html .= "<td>Internal oligo</td>\t";
            $html .= "<td>" . $results->{"PRIMER_INTERNAL_${index}_SEQUENCE"} . "</td>\t";
            $html .= "<td class='num'>$start</td>\t";
            $html .= "<td class='num'>$end</td>\t";
            $html .= "<td class='num'>$length</td>\t";
            $html .= "<td class='num'>" . sprintf("%.2f", $results->{"PRIMER_INTERNAL_${index}_TM"}) . "</td>\t";
            $html .= "<td class='num'>" . sprintf("%d", $results->{"PRIMER_INTERNAL_${index}_GC_PERCENT"}) . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_INTERNAL_${index}_SELF_ANY_TH"} . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_INTERNAL_${index}_SELF_END_TH"} . "</td>\n";
            $html .= "</tr>\n";
	    }

	    $html .= "<td>Product/Primer pair</td>\t";
            $html .= "<td class='empty'></td>\t";
            $html .= "<td class='num'>$product_start</td>\t";
            $html .= "<td class='num'>$product_end</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_PAIR_${index}_PRODUCT_SIZE"} . "</td>\t";
            $html .= "<td class='empty'></td>\t";
            $html .= "<td class='empty'></td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_PAIR_${index}_COMPL_ANY_TH"} . "</td>\t";
            $html .= "<td class='num'>" . $results->{"PRIMER_PAIR_${index}_COMPL_END_TH"} . "</td>\n";
	    $html .= "</tr>\n";
	    if ($index < ($pair_count-1)) {
		    $html .= "<tr><td class='sep' colspan='10'></td></tr>\n";
	    }
    }
    $html .= "</table>\n";
    return $html;
}

sub generate_dynamic_html {
    my $resultsHash = shift;
    print STDERR "in generate_dynamic_html: resultsHash = $resultsHash\n";
    print STDERR "num primers = $resultsHash->{PRIMER_PAIR_NUM_RETURNED}\n";
    if ($resultsHash->{PRIMER_PAIR_NUM_RETURNED} == 0) {
        my $html = "Problem: no primer pairs returned.\n";
        
        my $error_messages = "";
        for my $key (sort keys %$resultsHash) {
            $error_messages .= "<li>$key:&nbsp;&nbsp;$resultsHash->{$key}\n" if $key =~ /ERROR/i;
        }
        if ($error_messages) {
            $html .= "<p>Error messages from Primer3:<ul>\n";
            $html .= $error_messages;
            $html .= "</ul>\n";
        }
        return $html
    }

    my %html_param  = (
    "right_primer_color" => "#ff9999",
    "left_primer_color" => "#ffff66",
    "target_sequence_color" => "#ccffcc",
    "row_alt_color" => "#dddddd",
    "row_highlight_color" => "#9999ff" );

    my $html .= qq(<html>\n<head>
<style>
h3  {
    font-size: 1.2em;
    margin-bottom: .2em;
}

table {
    font-family:sans-serif; 
    margin: .2em;
}

table.numbers {
    text-align:right;
}

th {
    padding: .1em;
    text-align: left;
}

td {
    padding-left: .5em;
    padding-right: .5em;
}
.primer3plus_left_primer { background-color: $html_param{left_primer_color} }
.primer3plus_right_primer { background-color: $html_param{right_primer_color} }
.primer3plus_target_sequence { background-color: $html_param{target_sequence_color} }

</style>
</head>
<body>
);

    $html .= generate_javascript(\%html_param);

    $html .= generate_html_all_pairs_view($resultsHash, \%html_param);	

    $html .= "</body></html>\n";
}

sub generate_javascript {
    my $html_param = shift;
    return qq(
<script type="text/javascript">

let cur_pair_index = 0;
let right_primer_color="$html_param->{right_primer_color}";
let left_primer_color="$html_param->{left_primer_color}";
let alternate_row_colors = ["transparent", "$html_param->{row_alt_color}"];
let row_highlight_color="$html_param->{row_highlight_color}";

function update_current_pair() {
new_index = document.getElementById("select_pair_control").value;

old_row_color = alternate_row_colors[cur_pair_index % 2]
old_row = document.getElementById("PRIMER_LEFT_"+cur_pair_index);
old_row.style.backgroundColor=old_row_color;
old_row.children[4].style.backgroundColor="transparent";

old_row = document.getElementById("PRIMER_RIGHT_"+cur_pair_index);	
old_row.style.backgroundColor=old_row_color;
old_row.children[3].style.backgroundColor="transparent";

new_row = document.getElementById("PRIMER_LEFT_"+new_index);
new_row.style.backgroundColor=row_highlight_color;
new_row.children[4].style.backgroundColor=left_primer_color;

new_row = document.getElementById("PRIMER_RIGHT_"+new_index);
new_row.style.backgroundColor=row_highlight_color;
new_row.children[3].style.backgroundColor=right_primer_color;

document.getElementById("sequence_for_pair_"+cur_pair_index).style.display="none";
document.getElementById("sequence_for_pair_"+new_index).style.display="block";

document.getElementById("primer_pair_id").innerHTML=new_index;
document.getElementById("pair_data_"+cur_pair_index).style.backgroundColor="transparent";
document.getElementById("pair_data_"+new_index).style.backgroundColor=row_highlight_color;

cur_pair_index = new_index;
}

</script>
)
}

sub generate_html_all_pairs_view {
    my $results = shift;
    my $html_param = shift;
    my $pair_count = $results->{PRIMER_PAIR_NUM_RETURNED};
    my $html .= "<div id=\"primer3_results_container\" onload='update_current_pair(0)'>\n";
    $html .= "results: $results <br>\n";
    $html .= "num primers: $pair_count == $results->{PRIMER_PAIR_NUM_RETURNED} <br>\n";
    $html .= "<h3>Select Primer Pair:</b> <select id=\"select_pair_control\" onchange='update_current_pair(this.value)'>\n";
    for (my $i=0; $i < $pair_count; $i++) {
        $html .= "<option value=\"$i\">$i</option>\n";
    }
    $html .= "</select></h3>\n";

    $html .= "<table id='PRIMER_PAIR_TABLE' class='numbers'>\n";
    $html .= "<tr><th>#</th>\n";
    $html .= "<th>Dir</th>\n";
    $html .= "<th title=\"0-based start position in target sequence.\">Start</th>\n";
    $html .= "<th title=\"Length of primer in bases.\">Length</th>\n";
    $html .= "<th title=\"Sequence of the primer.\">Sequence</th>\n";
    $html .= "<th title=\"The melting temperature for the selected oligo.\">&nbsp;T<sub>m</sub> &nbsp;</th>\n";
    $html .= "<th title=\"Percent GC for the selected oligo.\">GC%</th>\n";
    $html .= "<th title=\"Tendency of a primer to bind to itself, interfering with target sequence binding.\">Any <br>Compl.</th>\n";
    $html .= "<th title=\"The tendency of the 3'-END to bind an identical primer, allowing it to form a primer-dimer.\">End <br>Compl.</th>";
    $html .= "<th title=\"Delta G of disruption of the five 3' bases of the primer.\">3' <br>Stab</th>\n";
    $html .= "<th title=\"Contribution of primer to pair penalty. Lower is better. See primer3.org/manual.html.\">Penalty</th></tr>\n";
    for my $primer_index (0 .. $pair_count-1) {
        my $row_color = $primer_index % 2 ? $html_param->{row_alt_color} : "transparent";
        $row_color = $html_param->{row_highlight_color} unless $primer_index;
        for my $dir ("LEFT", "RIGHT") {
            my $key_prefix = "PRIMER_${dir}_$primer_index";
            $html .= "<tr id=\"$key_prefix\" style=\"background-color:$row_color\">";
            if ($dir eq 'LEFT') {
                $html .= "<td rowspan='2'>$primer_index</td>"; 
            }
            $html .= "<td>".lc($dir)."</td>";
            my ($start, $length) = split(",", $results->{$key_prefix});
            $html .= "<td>$start</td>";
            $html .= "<td>$length</td>";
            my $seq_bg_color = "transparent";
            if ($primer_index == 0) {
                $seq_bg_color = $dir eq "LEFT" ? $html_param->{left_primer_color} : $html_param->{right_primer_color};
            }
            $html .= "<td style='background-color:$seq_bg_color; font-family:Courier,monospace'>" . $results->{$key_prefix . "_SEQUENCE"} . "</td>";
            $html .= "<td>" . sprintf("&nbsp;%.1f", $results->{$key_prefix . "_TM"}) . "</td>";
            $html .= "<td>" . sprintf("&nbsp;%.1f", $results->{$key_prefix . "_GC_PERCENT"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_ANY_TH"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_END_TH"}) . "</td>";
            $html .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_END_STABILITY"}) . "</td>";
            $html .= "<td>" . sprintf("%.2f", $results->{$key_prefix . "_PENALTY"}) . "</td>";
            $html .= "</tr>\n";
        }
    }
    $html .= "</table>\n";
    $html .= "<span style=\"font-size:0.7em\">Note that primer3 reports zero-based sequence start positions.</span>";

    $html .= "<h3>Statistics for Primer Pair: <span id='primer_pair_id'>0</span></h3>\n";
    $html .= "<table  onload='update_current_pair()' class='numbers'>\n";
    $html .= "<tr><th>#</th><th title=\"Length of the PCR product.\">Product <br>Size</th>\n";
    $html .= "<th title=\"Tendency of the 3'-ENDs of a primer pair to bind to each other and extend, forming a primer-dimer.\">End <br>Compl.</th>\n";
    $html .= "<th title=\"Tendency of a primer pair to bind to each other, interfering with target sequence binding.\">Any <br>Compl.</th>\n";
    $html .= "<th title=\"Lower is better. See primer3 manual for details (primer3.org/manual.html).\"<th>Penalty</th></tr>\n";
    for my $pair_index (0 .. $pair_count-1) {
        my $row_color = $pair_index % 2 ? $html_param->{row_alt_color} : "transparent";
        $row_color = $html_param->{row_highlight_color} unless $pair_index;
        $html .= "<tr id='pair_data_$pair_index' style='background-color:$row_color'>";
        my $prefix = "PRIMER_PAIR_$pair_index";
        $html .= "<td>$pair_index</td>";
        $html .= "<td>".$results->{$prefix . "_PRODUCT_SIZE"}."</td>";
        $html .= "<td>". sprintf("%.1f", $results->{$prefix . "_COMPL_END_TH"}) ."</td>";
        $html .= "<td>". sprintf("%.1f", $results->{$prefix . "_COMPL_ANY_TH"}) ."</td>";
        $html .= "<td>". sprintf("%.2f", $results->{$prefix . "_PENALTY"}) ."</td></tr>\n";
    }
    $html .= "</table>\n";
    $html .= "<h3>Primers in Template Sequence</h3>\n";
    $html .= "<div id=\"sequence_container\">\n";
    $html .= "<span style=\"font-size:0.7em\">Color key:\n";
    $html .= "&nbsp;&nbsp<a class='primer3plus_left_primer'>Left primer</a>";
    if ($results->{SEQUENCE_TARGET}) {
        $html .= "&nbsp;&nbsp;<a class='primer3plus_target_sequence'>Target sequence</a>";
    }
    $html .= "&nbsp;&nbsp;<a class='primer3plus_right_primer'>Right primer</a>";
    $html .= "</span>\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
        $html .= create_primer_pair_on_sequence_html($results, $index, $index == 0, $html_param); 
    } 
    $html .= "</div></div>\n"; 
    return $html;
}


sub generate_html_one_pair_view {
    my $results = shift;
    my $pair_count = $results->{PRIMER_PAIR_NUM_RETURNED};
    my $html .= "<b>Show Primer Pair:</b> <select id=\"select_pair_control\" onchange=\"update_current_pair(this.value)\">\n";
    for (my $i=0; $i < $pair_count; $i++) {
        $html .= "<option value=\"$i\">$i</option>\n";
    }
    $html .= "</select><p>\n";
    #$html .= "<span id=\"primer_pair_label\" style=\"font-weight: bold\">Showing Primer Pair 0</span>\n";
    $html .= "<div id=\"primer_pair_container\">\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
      $html .= create_primer_pair_table_html($results, $index, $index == 0); 
    } 
    $html .= "</div>\n"; 

    $html .= "<p>Primers in template sequence:<div id=\"sequence_container\" style=\"font-family:Courier,monospace\">\n";
    #now add html table for each pair to javascript variable
    for (my $index = 0; $index < $pair_count; $index++) {
      $html .= create_primer_pair_on_sequence_html($results, $index, $index == 0); 
    } 
    $html .= "</div>\n"; 
    return $html
}

    

sub create_primer_pair_table_html {
  my ($results, $primer_index, $visibility) = @_; 

  my $tableHTML .= "<table class='bvbrc_primer_table' id='PRIMER_PAIR_TABLE_$primer_index' style=\"display:";
  $tableHTML .= $primer_index eq 0 ? "block" : "none";  # first one is visible to start (changed by javascript on user interaction)
  $tableHTML .= "\">\n";
  $tableHTML .= "<tr><th>Dir</th><th>Start</th><th>Length</th><th>Sequence</th><th>Tm</th><th>GC</th><th>Any</th><th>End</th><th>3' Stab</th><th>Penalty</th></tr>\n";
  for my $dir ("LEFT", "RIGHT") {
      my $key_prefix = "PRIMER_${dir}_$primer_index";
      $tableHTML .= "<tr>";
      $tableHTML .= "<td>$dir</td>";
      my ($start, $length) = split(",", $results->{$key_prefix});
      $tableHTML .= "<td>$start</td>";
      $tableHTML .= "<td>$length</td>";
      my $color = $dir eq "LEFT" ? "#ccccff" : "rgb(250,240,75)";
      $tableHTML .= "<td style=\"background-color:$color\">" . $results->{$key_prefix . "_SEQUENCE"} . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_TM"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_GC_PERCENT"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_ANY"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_SELF_END"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_END_STABILITY"}) . "</td>";
      $tableHTML .= "<td>" . sprintf("%.1f", $results->{$key_prefix . "_PENALTY"}) . "</td>";
      $tableHTML .= "</tr>\n";
  }
  my $product_length = $results->{"PRIMER_PAIR_${primer_index}_PRODUCT_SIZE"};
  $tableHTML .= "<tr><td colspan=\"10\"><span style=\"font-weight:bold\">Product Length:</span> $product_length</td></tr>\n";
  $tableHTML .= "</table>\n";
  return $tableHTML;
}

###################################################
# divHTMLsequence: Prints out the sequence nicely # 
###################################################
sub create_primer_pair_on_sequence_html {
  my ($results, $sequence, $seqLength, $firstBase);
  my ($base, $count, $preCount, $postCount, $printBase) ;
  my ($format, $baseFormat, $preFormat, $pair_to_show);
  my (@targets, $region, $run, $counter, $madeRegion);
  my $tableHTML;
   
  $results = shift;
  $pair_to_show = shift;
  my $is_visible = shift;
  my $html_param = shift;

  $sequence = $results->{"SEQUENCE_TEMPLATE"};
  $format = $sequence;
  $format =~ s/\w/N/g;

  $seqLength = length ($sequence);
  $firstBase = 1; #$results->{"PRIMER_FIRST_BASE_INDEX"};
  
  if ((defined $results->{"SEQUENCE_TEMPLATE"})and ($results->{"SEQUENCE_TEMPLATE"} ne "")) {

  if (defined ($results->{"SEQUENCE_EXCLUDED_REGION"}) and (($results->{"SEQUENCE_EXCLUDED_REGION"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_EXCLUDED_REGION"};
      foreach $region (@targets) {
          $format = addRegion($format,$region,$firstBase,"E");
      }
  }
  if (defined ($results->{"SEQUENCE_TARGET"}) and (($results->{"SEQUENCE_TARGET"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_TARGET"};
      foreach $region (@targets) {
          $format = addRegion($format,$region,$firstBase,"T");
      }
  }
  if (defined ($results->{"SEQUENCE_INCLUDED_REGION"}) and (($results->{"SEQUENCE_INCLUDED_REGION"}) ne "")) {
      $format = addRegion($format,$results->{"SEQUENCE_INCLUDED_REGION"},$firstBase,"I");
  }
  
  # Add only the specified primer pair e.g. Detection
  if ($pair_to_show >= 0) {
      if (defined ($results->{"PRIMER_LEFT_" . $pair_to_show}) and (($results->{"PRIMER_LEFT_" . $pair_to_show}) ne "")) {
           $format = addRegion($format,$results->{"PRIMER_LEFT_" . $pair_to_show},$firstBase,"F");
      }
      if (defined ($results->{"PRIMER_INTERNAL_" . $pair_to_show}) and (($results->{"PRIMER_INTERNAL_" . $pair_to_show}) ne "")) {
     	   $format = addRegion($format,$results->{"PRIMER_INTERNAL_" . $pair_to_show},$firstBase,"O");
      }
      if (defined ($results->{"PRIMER_RIGHT_" . $pair_to_show}) and (($results->{"PRIMER_RIGHT_" . $pair_to_show}) ne "")) {
           $format = addRegion($format,$results->{"PRIMER_RIGHT_" . $pair_to_show},$firstBase,"R");
      }
  }
  # Mark no primers on the sequence e.g. Primer list
  elsif ($pair_to_show eq -1) {
      
  }
  if (defined ($results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"}) and (($results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"}) ne "")) {
      @targets = split ' ', $results->{"SEQUENCE_OVERLAP_JUNCTION_LIST"};
      foreach $region (@targets) {
          $madeRegion = $region - 1;
          $madeRegion .= ",2";
          $format = addRegion($format,$madeRegion,$firstBase,"-");
      }
  }

  ## Handy for testing:
  #  $sequence = $format;

  my $display_style = $is_visible ? "block" : "none";
  $tableHTML = "<table id=\"sequence_for_pair_$pair_to_show\" style=\"display: $display_style; font-family:Courier,monospace\">";
  $tableHTML .= qq{
     <colgroup>
       <col width="13%" style="text-align: right;">
       <col width="87%">
     </colgroup>
};
  $preFormat = "N";
  for (my $i=0; $i<$seqLength; $i++) {
     $count = $i;
     $preCount = $i - 1;
     $postCount = $i + 1;
     $base = substr($sequence,$i,1);
     $baseFormat = substr($format,$i,1);

     if (($count % 50) eq 0) {
         $printBase = $i + $firstBase;
         $tableHTML .= qq{     <tr>
       <td>$printBase&nbsp;&nbsp;</td>
       <td>};
     }

     if ((($count % 10) eq 0) and !(($count % 50) eq 0)) {
         $tableHTML .= qq{&nbsp;&nbsp;};
     }

     if ($preFormat ne $baseFormat) {
         if ($preFormat ne "J") {
             $tableHTML .= qq{</a>};
         }
         if ($baseFormat eq "N") {
             $tableHTML .= qq{<a>};
         }
         if ($baseFormat eq "E") {
             $tableHTML .= qq{<a class="primer3plus_excluded_region">};
         }
         if ($baseFormat eq "T") {
             $tableHTML .= qq{<a class="primer3plus_target_sequence">};
         }
         if ($baseFormat eq "I") {
             $tableHTML .= qq{<a class="primer3plus_included_region">};
         }
         if ($baseFormat eq "F") {
             $tableHTML .= qq{<a class="primer3plus_left_primer">};
         }
         if ($baseFormat eq "O") {
             $tableHTML .= qq{<a class="primer3plus_internal_oligo">};
         }
         if ($baseFormat eq "R") {
             $tableHTML .= qq{<a class="primer3plus_right_primer">};
         }
         if ($baseFormat eq "B") {
             $tableHTML .= qq{<a class="primer3plus_left_right_primer">};
         }
         if ($baseFormat eq "-") {
             $tableHTML .= qq{<a class="primer3plus_primer_overlap_pos">};
         }

     }

     $tableHTML .= qq{$base};

     if (($postCount % 50) eq 0) {
         $tableHTML .= qq{</a></td>
     </tr>
};
         $baseFormat = "J";
     }
     $preFormat = $baseFormat;
  }

  if (($postCount % 50) ne 0) {
      $tableHTML .= qq{</a></td>
     </tr>
};
  }

  $tableHTML .= "  </table>\n";
}

  return $tableHTML;
}

sub addRegion {
  my ($formatString, $region, $firstBase, $letter);
  my ($regionStart, $regionLength, $regionEnd);
  my ($stingStart, $stringRegion, $stringEnd, $stringLength);
  $formatString = shift;
  $region = shift;
  $firstBase = shift;
  $letter = shift;
  ($regionStart, $regionLength) = split "," , $region;
  $regionStart =~ s/\s//g;
  $regionLength =~ s/\s//g;
  $regionStart = $regionStart - $firstBase;
  if ($regionStart < 0) {
      $regionStart = 0;
  }
  if ($letter eq "R") {
  $regionStart = $regionStart - $regionLength + 1;  
  }
  $regionEnd = $regionStart + $regionLength;
  $stringLength = length $formatString;
  if ($regionEnd > $stringLength) {
      $regionEnd = $stringLength;
  }
  $stingStart = substr($formatString,0,$regionStart);
  $stringRegion = substr($formatString,$regionStart,$regionLength);
  $stringEnd = substr($formatString,$regionEnd,$stringLength);
  if (($letter ne "F") and ($letter ne "R")) {
      $stringRegion =~ s/\w/$letter/g;
  }
  if ($letter eq "F") {
      $stringRegion =~ tr/NETIFROB/FFFFFBFB/;
  }
  if  ($letter eq "R") {
      $stringRegion =~ tr/NETIFROB/RRRRBRRB/;
  }
  $formatString = $stingStart;
  $formatString .= $stringRegion;
  $formatString .= $stringEnd;

  return $formatString;
}

sub handle_sequence_region_markup {
    print STDERR "in handle_sequence_region_markup\n" if $debug;
    my ($primer3_params) = @_;
    my %region_types;
    my $seqpos = 0;
    my $seq_target_start = -1;
    my $included_region_start = -1;
    my $excluded_region_start = -1;
    for my $char (split '', $primer3_params->{SEQUENCE_TEMPLATE}) {
        if ($char eq '[') {
            $seq_target_start = $seqpos;
        }
        elsif ($char eq ']') {
            if ($seq_target_start >= 0) {
                my $length = $seqpos - $seq_target_start + 1;
                $primer3_params->{SEQUENCE_TARGET} .= " $seq_target_start,$length ";
                print STDERR "updated SEQUENCE_TARGET to $primer3_params->{SEQUENCE_TARGET}\n" if $debug;
            }
            $seq_target_start = -1;
        }
        if ($char eq '{') {
            $included_region_start = $seqpos;
        }
        elsif ($char eq '}') {
            if ($included_region_start >= 0) {
                my $length = $seqpos - $included_region_start + 1;
                $primer3_params->{SEQUENCE_INCLUDED_REGION} = "$included_region_start,$length"; #there can be only one
                print STDERR "updated SEQUENCE_INCLUDED_REGION to $primer3_params->{SEQUENCE_INCLUDED_REGION}\n" if $debug;
            }
            $included_region_start = -1;
        } 
        elsif ($char eq '<') {
            $excluded_region_start = $seqpos;
        }
        elsif ($char eq '>') {
            if ($excluded_region_start >= 0) {
                my $length = $seqpos - $excluded_region_start + 1;
                $primer3_params->{SEQUENCE_EXCLUDED_REGION} .= " $excluded_region_start,$length ";
                print STDERR "updated SEQUENCE_EXCLUDED_REGION to $primer3_params->{SEQUENCE_EXCLUDED_REGION}\n" if $debug;
            }
            $excluded_region_start = -1;
        }
        elsif ($char eq '-') {
            $primer3_params->{SEQUENCE_OVERLAP_JUNCTION_LIST} .= " $seqpos ";
            print STDERR "updated SEQUENCE_OVERLAP_JUNCTION_LIST to $primer3_params->{SEQUENCE_OVERLAP_JUNCTION_LIST}\n" if $debug;
        }
        elsif ($char =~ /[ACGT]/i) {
            $seqpos++;
        } 
    }
    # now remove all non-nucleotide characters
    $primer3_params->{SEQUENCE_TEMPLATE} =~ s/[^ACGTacgt]//g;
}

sub run_cmd {
    my ($cmd) = @_;
    my ($out, $err);
    run($cmd, '>', \$out, '2>', \$err)
        or die "Error running cmd=@$cmd, stdout:\n$out\nstderr:\n$err\n";
    # print STDERR "STDOUT:\n$out\n";
    # print STDERR "STDERR:\n$err\n";
    return ($out, $err);
}
