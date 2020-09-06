#Parser for dot-parentheses format
#It returns primitive pseudoknot objects, base sequence and dot-parentheses array

package DPParser;

use strict;

use constant DOT => '.';

sub parse{
    my (undef, $dp_file_path) = @_;

    my $primitive_pseudoknots = [];
    my ($base_seq_str, $secondary_structure) = ('', '');

    open (DP, "<$dp_file_path") or die "Cannot open file at $dp_file_path";
    while (<DP>) {
	if ($_ =~ /^([A-Za-z]+)[\r\n]*$/) {
	    $base_seq_str = $base_seq_str . $1;
	}
	elsif ($_ =~ /^([\.\(\)\[\]\{\}<>A-Za-z]+)[\r\n]*$/) {
	    $secondary_structure = $secondary_structure . $1;
	}
	elsif ($_ !~ /^#.*/ && $_ !~ /^\s+/) {
	    die "Unknown input: $_";
	}
    }

    close DP or die "Cannot close file at $dp_file_path";

    if ($base_seq_str eq '') {
	die 'Base sequence is missing';
    }

    if ($secondary_structure eq '') {
	die 'Secondary structure is missing';
    }

    if (length($base_seq_str) != length($secondary_structure)) {
	die 'Base sequence length not equal to secondary structure length';
    }

    #Group the base pairs into base pair stems
    my ($stem_outermost_pairs, $stems, $paired_pos_ptrs, $structure_symbols) = _group_to_stems($secondary_structure);
    #Extract primitive pseudoknots from the base pair stems
    my $primitive_pseudoknots = PrimitivePseudoknotExtractor->extract($stem_outermost_pairs, $stems, $paired_pos_ptrs);
    my @base_seq = split(//, $base_seq_str);

    return $primitive_pseudoknots, \@base_seq, $structure_symbols, $base_seq_str;
}

sub _group_to_stems {
    my $secondary_structure = shift;
    my $stems = {};
    my ($stem_outermost_pairs, $stem, $outermost_base_pair) = ([], [], []);
    my $paired_pos_ptrs = [];
    my $unsettled_bracket_upstream_pos = {};
    my $next_paired_pos = {};
    my $last_paired_pos = 0;

    my @structure_symbols = split(//, $secondary_structure);
    my $structure_length = scalar @structure_symbols;

    for (my $i = 0; $i < $structure_length; $i++) {
	my $symbol = $structure_symbols[$i];
	if ($symbol eq DOT) {
	    next;
	}
	elsif (BracketPairs->is_open_bracket($symbol)) {
	    my $unsettled_upstream_pos = $unsettled_bracket_upstream_pos->{$symbol};
	    if (!defined($unsettled_upstream_pos)) {
		$unsettled_upstream_pos = [];
		$unsettled_bracket_upstream_pos->{$symbol} = $unsettled_upstream_pos;
	    }

	    my $curr_upstream_pos = $i + 1;
	    push @{$unsettled_upstream_pos}, $curr_upstream_pos;

	    if (defined($outermost_base_pair->[0])) {
		($stem_outermost_pairs, $stems, $outermost_base_pair, $stem) = _add_to_stems($stem_outermost_pairs, $stems, $outermost_base_pair, $stem);
	    }

	    $next_paired_pos->{$last_paired_pos} = $curr_upstream_pos;
	    $last_paired_pos = $curr_upstream_pos;
	}
	else {
	    my $pair_open_bracket = BracketPairs->get_open_bracket($symbol);
	    my $unsettled_upstream_pos = $unsettled_bracket_upstream_pos->{$pair_open_bracket};
	    if (defined($unsettled_upstream_pos) && defined($unsettled_upstream_pos->[0])) {
		my $paired_upstream_pos = pop @{$unsettled_upstream_pos};
		my $curr_downstream_pos = $i + 1;

		if (defined($outermost_base_pair->[0])) {
		    if ($next_paired_pos->{$paired_upstream_pos} != $outermost_base_pair->[0]) {
			($stem_outermost_pairs, $stems, $outermost_base_pair, $stem) = _add_to_stems($stem_outermost_pairs, $stems, $outermost_base_pair, $stem);
		    }

		    $outermost_base_pair = [$paired_upstream_pos, $curr_downstream_pos];
		    unshift @{$stem}, $outermost_base_pair;
		}
		else {
		    $outermost_base_pair = [$paired_upstream_pos, $curr_downstream_pos];
		    $stem = [$outermost_base_pair];
		}

		$paired_pos_ptrs->[$paired_upstream_pos] = $curr_downstream_pos;
		$paired_pos_ptrs->[$curr_downstream_pos] = $paired_upstream_pos;

		$next_paired_pos->{$last_paired_pos} = $curr_downstream_pos;
		$last_paired_pos = $curr_downstream_pos;
	    }
	    else {
		die "Closing bracket $symbol not paired\n";
	    }
	}
    }

    if (!_is_all_open_bracket_settled($unsettled_bracket_upstream_pos)) {
	die "Unpaired open bracket remains\n";
    }

    if (defined($outermost_base_pair->[0])) {
	($stem_outermost_pairs, $stems, undef, undef) = _add_to_stems($stem_outermost_pairs, $stems, $outermost_base_pair, $stem);
    }

    my @sorted_outermost_pairs = sort {$a->[0] <=> $b->[0]} @{$stem_outermost_pairs};

    return (\@sorted_outermost_pairs, $stems, $paired_pos_ptrs, \@structure_symbols);
}

sub _add_to_stems {
    my ($stem_outermost_pairs, $stems, $stem_outermost_pair, $stem) = @_;

    $stems->{$stem_outermost_pair->[0]} = $stem;
    push @{$stem_outermost_pairs}, $stem_outermost_pair;

    return ($stem_outermost_pairs, $stems, [], []);
}

sub _is_all_open_bracket_settled {
    my $unsettled_open_bracket_pos = shift;

    foreach (values %{$unsettled_open_bracket_pos}) {
	if (defined($_->[0])) {
	    return 0;
	}
    }

    return 1;
}

1;
