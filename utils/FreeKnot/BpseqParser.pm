#Parser for BPSEQ format
#It returns primitive pseudoknot objects, base sequence and paired positions

package BpseqParser;

use strict;

sub parse {
    my (undef, $bpseq_file_path) = @_;

    my ($base_seq, $paired_pos_ptrs) = ([], []);
    my ($next_paired_pos, $prev_paired_pos) = ({}, {});
    my $matched_pos = {};
    my $last_paired_pos = 0;
    my $base_count = 0;

    open (BPSEQ, "<$bpseq_file_path") or die "Cannot open file at $bpseq_file_path";

    while (<BPSEQ>) {
	if ($_ =~ /^([0-9]+) ([A-Za-z]{1}) ([0-9]+)[\r\n]*$/) {
	    my ($pos, $base, $paired_pos) = ($1, $2, $3);
	    if ($pos != ++$base_count) {
		die "Base position $base_count is missing";
	    }

	    if ($paired_pos > 0) {
		if ($pos < $paired_pos) {
		    $matched_pos->{$pos} = $paired_pos;
		}
		else {
		    if ($matched_pos->{$paired_pos} != $pos) {
			die "Unmatched pair position $pos and $paired_pos";
		    }
		}

		$next_paired_pos->{$last_paired_pos} = $pos;
		$prev_paired_pos->{$pos} = $last_paired_pos;
		$last_paired_pos = $pos;
	    }

	    $paired_pos_ptrs->[$pos] = $paired_pos;
	    $base_seq->[$pos - 1] = $base;
	}
	elsif ($_ !~ /^#.*/ && $_ !~ /^\s+/) {
	    die "Unknown input: $_";
	}
    }

    $next_paired_pos->{$last_paired_pos} = 0;
    $prev_paired_pos->{0} = $last_paired_pos;

    close BPSEQ or die "Cannot close file at $bpseq_file_path";

    #Group the base pairs into base pair stems
    my ($stem_outermost_pairs, $stems) = _group_to_stems($next_paired_pos, $prev_paired_pos, $paired_pos_ptrs);
    #Extract primitive pseudoknots from the base pair stems
    my $primitive_pseudoknots = PrimitivePseudoknotExtractor->extract($stem_outermost_pairs, $stems, $paired_pos_ptrs);

    return ($primitive_pseudoknots, $base_seq, $paired_pos_ptrs, $base_count);
}

sub _group_to_stems {
    my ($next_paired_pos, $prev_paired_pos, $paired_pos_ptrs) = @_;

    my $stems = {};
    my $stem_outermost_pairs = [];
    my $stem;
    my $last_pair;

    my $curr_pos = $next_paired_pos->{0};
    while ($curr_pos > 0) {
	my $paired_pos = $paired_pos_ptrs->[$curr_pos];
	if ($paired_pos < $curr_pos) {
	    undef $last_pair;
	    $curr_pos = $next_paired_pos->{$curr_pos};
	    next;
	}

	my $curr_pair = [$curr_pos, $paired_pos];

	if (defined($last_pair) && $prev_paired_pos->{$last_pair->[1]} == $paired_pos) {
	    push @{$stem}, $curr_pair;
	}
	else {
	    $stem = [$curr_pair];
	    $stems->{$curr_pos} = $stem;
	    push @{$stem_outermost_pairs}, $curr_pair;
	}

	$last_pair = $curr_pair;
	$curr_pos = $next_paired_pos->{$curr_pos};
    }

    return ($stem_outermost_pairs, $stems);
}

1;
