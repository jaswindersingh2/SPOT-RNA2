#Writer for dot-parentheses format

package DPWriter;

use strict;

use constant DOT => '.';
use constant OPEN_BRACKET => '(';
use constant CLOSE_BRACKET => ')';
use constant TEMP_DP_FILE => 'MWIS_temp.dp';

sub output_results {
    my (undef, $combined_base_pair_removal_pos, $structure_symbols, $base_seq_str) = @_;

    if (@{$combined_base_pair_removal_pos} == 0) {
	my $output_structure = join('', @{$structure_symbols});
	print "$base_seq_str\n$output_structure\n";
    }

    foreach (@{$combined_base_pair_removal_pos}) {
	my $output_structure = '';
	for (my $i = 0; $i < @{$structure_symbols}; $i++) {
	    if (exists($_->{$i + 1})) {
		$output_structure = $output_structure . DOT;
	    }
	    else {
		$output_structure = $output_structure . $structure_symbols->[$i];
	    }
	}

	print "$base_seq_str\n$output_structure\n";
    }
}

sub output_mfe_candidate {
    my (undef, $base_pair_removal_pos, $paired_pos_ptrs, $structure_symbols, $base_seq_str) = @_;

    my $base_seq_len = length($base_seq_str);
    my $output_structure = '';
    if (defined($paired_pos_ptrs)) {
	for (my $i = 1; $i <= $base_seq_len; $i++) {
	    if (exists($base_pair_removal_pos->{$i})) {
		$output_structure = $output_structure . DOT;
	    }
	    else {
		my $paired_pos = $paired_pos_ptrs->[$i];
		if ($paired_pos == 0) {
		    $output_structure = $output_structure . DOT;
		}
		elsif ($i < $paired_pos) {
		    $output_structure = $output_structure . OPEN_BRACKET;
		}
		else {
		    $output_structure = $output_structure . CLOSE_BRACKET;
		}
	    }
	}
    }
    elsif (defined($structure_symbols)) {
	for (my $i = 1; $i <= $base_seq_len; $i++) {
	    if (exists($base_pair_removal_pos->{$i})) {
		$output_structure = $output_structure . DOT;
	    }
	    else {
		$output_structure = $output_structure . $structure_symbols->[$i - 1];
	    }
	}
    }

    $output_structure =~ s/[\[\{<A-Z]/\(/g;
    $output_structure =~ s/[\]\}>a-z]/\)/g;

    open (DP, ">" . TEMP_DP_FILE) or die "Cannot open file at " . TEMP_DP_FILE;
    print DP "$base_seq_str\n$output_structure\n";
    close DP or die "Cannot close file at " . TEMP_DP_FILE;
}

1;
