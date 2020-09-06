#Main program for pseudoknot removal
#It accepts input RNA secondary structure as BPSEQ format or dot-parentheses format
#There are four choices of scoring functions: No. of base pairs, no. of stems, no. of hydrogen
#bonds, and Turner free energy (Turner, D. H. & Mathews, D. H., NAR 2009)). The optimization goal
#for the first three options is to maximize the score as all the choices only give positive values.
#For the last option, the goal is to minimize the score (i.e. free energy).

#!/usr/bin/perl

use BpseqParser;
use BpseqWriter;
use BracketPairs;
use ChordModel;
use CircleGraph;
use DPParser;
use DPWriter;
use MIS;
use MWIS;
use PrimitivePseudoknotExtractor;
use ScoringFunctions;
use VertexSubset;
use strict;

#OS_BIT specifies the length of a bitstring used in the circle graph
use constant OS_BIT => 32;

if (@ARGV < 5) {
    print "Usage: perl $0 -i <input file format> -s <scoring function> <input file path> [-a : report all optimal solutions]\n";
    exit;
}

my ($input_file_path, $input_file_format, $scoring_fx_option);
my $is_report_all = 0;

for (my $i = 0; $i < @ARGV; $i++) {
    if ($ARGV[$i] eq '-i') {
	if (defined($input_file_format)) {
	    print "Duplicated input file format specification\n";
	    exit;
	}
	else {
	    $input_file_format = $ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '-s') {
	if (defined($scoring_fx_option)) {
	    print "Duplicated scoring function specification\n";
	    exit;
	}
	else {
	    $scoring_fx_option = $ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '-a') {
	$is_report_all = 1;
    }
    elsif (substr($ARGV[$i], 0, 1) eq '-') {
	print "Unknown parameter $ARGV[$i]\n";
	exit;
    }
    elsif (!defined($input_file_path)) {
	$input_file_path = $ARGV[$i];
    }
}

if (!defined($input_file_path)) {
    print "No input file path specified\n";
    exit;
}

#Select the scoring function according to the user option. It will be used to calculate the score of
#each stem in the MWIS algorithm
my ($scoring_function, $criteria, $is_fe) = ScoringFunctions->get_scoring_function($scoring_fx_option);
if (!defined($scoring_function)) {
    print "Unknown scoring function specified: $scoring_fx_option\n";
    exit;
}

my ($primitive_pseudoknots, $base_seq, $paired_pos_ptrs, $base_count, $structure_symbols, $base_seq_str);

#Parse the input structure file to generate pseudoknot objects
if ($input_file_format eq 'bpseq') {
    ($primitive_pseudoknots, $base_seq, $paired_pos_ptrs, $base_count) = BpseqParser->parse($input_file_path);
}
elsif ($input_file_format eq 'dp') {
    ($primitive_pseudoknots, $base_seq, $structure_symbols, $base_seq_str) = DPParser->parse($input_file_path);
}
else {
    print "Unknown input file format: $input_file_format\n";
    exit;
}

my $pseudoknot_base_pair_removal_pos = [];
my $prim_pseudoknot_count = 0;

#If free energy is selected as the scoring function, then MIS algorithm is applied to generate
#all MISs of the circle graph, and evaluated the free energy for each of them
if ($is_fe) {
    foreach (@{$primitive_pseudoknots}) {
	my $circle_graph = CircleGraph->new($_, OS_BIT);
	my $miss = MIS->get_mis($circle_graph, $criteria);
	my $base_pair_removal_pos = convert_to_base_pair_removal_pos_circle_graph($circle_graph, $miss);
	push @{$pseudoknot_base_pair_removal_pos}, $base_pair_removal_pos;
	$prim_pseudoknot_count++;
    }
}
#For other scoring function options, MWIS algorithm is applied to generate one/all MWISs from
#the chord model of the circle graph
else{
    foreach (@{$primitive_pseudoknots}) {
	my $chord_model = ChordModel->new($_);
	my $mwiss = MWIS->get_mwis($chord_model, $base_seq, $scoring_function, $criteria, $is_report_all);
	my $base_pair_removal_pos = convert_to_base_pair_removal_pos($chord_model, $mwiss);
	push @{$pseudoknot_base_pair_removal_pos}, $base_pair_removal_pos;
	$prim_pseudoknot_count++;
    }
}

#Combine the possible removal positions sets for all primitive pseudoknots
my $combined_base_pair_removal_pos = combine_base_pair_removal_pos($pseudoknot_base_pair_removal_pos, []);

#Determine the free energy of every structure converted from the MISs combinations of different
#primitive pseudoknots in the structure. It writes the structure to a temporary file and call
#RNAeval in ViennaRNA package to calculate its free energy
if ($is_fe) {
    my $mfe;
    my $mfe_base_pair_models = [];

    if (!defined($base_seq_str)) {
	$base_seq_str = join('', @{$base_seq});
    }

    foreach (@{$combined_base_pair_removal_pos}) {
	DPWriter->output_mfe_candidate($_, $paired_pos_ptrs, $structure_symbols, $base_seq_str);
	my $rna_eval_output = `RNAeval < MWIS_temp.dp`;
	$rna_eval_output =~ /(-?\d+\.\d+)/;
	if ($1 < $mfe || !defined($mfe)) {
	    $mfe_base_pair_models = [$_];
	    $mfe = $1;
	}
	elsif ($1 == $mfe) {
	    push @{$mfe_base_pair_models}, $_;
	}
    }

    $combined_base_pair_removal_pos = $mfe_base_pair_models;
}

if ($input_file_format eq 'bpseq') {
    BpseqWriter->output_results($combined_base_pair_removal_pos, $base_seq, $paired_pos_ptrs, $base_count);
}
elsif ($input_file_format eq 'dp') {
    DPWriter->output_results($combined_base_pair_removal_pos, $structure_symbols, $base_seq_str);
}

sub convert_to_base_pair_removal_pos_circle_graph {
    my ($circle_graph, $miss) = @_;

    my $base_pair_removal_pos = [];

    foreach my $mis (@{$miss}) {
	my $removed_vertex_nums = [];
	for (my $i = 0; $i < $mis->[0]; $i++) {
	    push @{$removed_vertex_nums}, $i;
	}

	for (my $i = 1; $i < @{$mis}; $i++) {
	    for (my $j = $mis->[$i - 1] + 1; $j < $mis->[$i]; $j++) {
		push @{$removed_vertex_nums}, $j;
	    }
	}

	for (my $i = $mis->[-1] + 1; $i < $circle_graph->get_vertex_count(); $i++) {
	    push @{$removed_vertex_nums}, $i;
	}

	my $removal_pos = {};
	foreach (@{$removed_vertex_nums}) {
	    my $vertex_attrs = $circle_graph->get_vertex_attrs_at($_);
	    my $stem_pairs = $vertex_attrs->{stem_pairs};
	    foreach (@{$stem_pairs}) {
		my ($pair_upstream_pos, $pair_downstream_pos) = @{$_};
		$removal_pos->{$pair_upstream_pos} = 1;
		$removal_pos->{$pair_downstream_pos} = 1;
	    }
	}

	push @{$base_pair_removal_pos}, $removal_pos;
    }

    return $base_pair_removal_pos;
}

sub convert_to_base_pair_removal_pos {
    my ($chord_model, $mwiss) = @_;

    my $base_pair_removal_pos = [];

    foreach my $mwis (@{$mwiss}) {
	my %removed_chord_edges = %{$chord_model->get_chord_edges()};
	foreach (@{$mwis}) {
	    delete $removed_chord_edges{$_->[0] . '-' . $_->[1]};
	}

	my $removal_pos = {};
	foreach my $removed_chord_edge (values %removed_chord_edges) {
	    my $removed_chord_base_pairs = $chord_model->get_chord_base_pairs($removed_chord_edge->[0], $removed_chord_edge->[1]);
	    foreach (@{$removed_chord_base_pairs}) {
		$removal_pos->{$_->[0]} = 1;
		$removal_pos->{$_->[1]} = 1;
	    }
	}

	push @{$base_pair_removal_pos}, $removal_pos;
    }

    return $base_pair_removal_pos;
}

sub combine_base_pair_removal_pos {
    my ($pseudoknot_base_pair_removal_pos, $combined_base_pair_removal_pos) = @_;

    my $expanded_base_pair_removal_pos = [];
    my $base_pair_removal_pos = pop @{$pseudoknot_base_pair_removal_pos};
    foreach my $removal_pos (@{$base_pair_removal_pos}) {
	if (defined($combined_base_pair_removal_pos->[0])) {
	    foreach (@{$combined_base_pair_removal_pos}) {
		my %expanded_removal_pos = (%{$removal_pos}, %{$_});
		push @{$expanded_base_pair_removal_pos}, \%expanded_removal_pos;
	    }
	}
	else {
	    push @{$expanded_base_pair_removal_pos}, $removal_pos;
	}
    }

    if (defined($pseudoknot_base_pair_removal_pos->[0])) {
	$expanded_base_pair_removal_pos = combine_base_pair_removal_pos($pseudoknot_base_pair_removal_pos, $expanded_base_pair_removal_pos);
    }

    return $expanded_base_pair_removal_pos;
}
