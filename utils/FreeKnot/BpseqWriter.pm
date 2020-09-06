#Writer for BPSEQ format

package BpseqWriter;

use strict;

sub output_results {
    my (undef, $combined_base_pair_removal_pos, $base_seq, $paired_pos_ptrs, $base_count) = @_;

    if (@{$combined_base_pair_removal_pos} == 0) {
	for (my $i = 1; $i <= $base_count; $i++) {
	    print $i . ' ' . $base_seq->[$i - 1] . ' ' . $paired_pos_ptrs->[$i] . "\n";
	}
    }

    foreach (@{$combined_base_pair_removal_pos}) {
	for (my $i = 1; $i <= $base_count; $i++) {
	    print $i . ' ' . $base_seq->[$i - 1] . ' ';
	    if (exists($_->{$i})) {
		print "0\n";
	    }
	    else {
		print $paired_pos_ptrs->[$i] . "\n";
	    }
	}
    }
}

1;
