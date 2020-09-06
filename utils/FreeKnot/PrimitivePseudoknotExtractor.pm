#Module that extracts primitive pseudoknots from all the base pair stems of the RNA secondary structure

package PrimitivePseudoknotExtractor;

use strict;

sub extract {
    my (undef, $stem_outermost_pairs, $stems, $paired_pos_ptrs) = @_;

    #Group together the crossing stems of a pseudoknot
    my ($knotted_pair_pos_groups, $outermost_pair_crossings) = _group_knotted_outermost_pairs($stem_outermost_pairs);
    #Create the pseudoknot objects
    my $primitive_pseudoknots = _get_prim_pseudoknots($stems, $knotted_pair_pos_groups, $outermost_pair_crossings, $paired_pos_ptrs);

    return $primitive_pseudoknots;
}

sub _group_knotted_outermost_pairs {
    my $stem_outermost_pairs = shift;

    my $knotted_pair_pos_groups = [];
    my $outermost_pair_crossings = {};
    my $paired_pos_to_group_id = {};
    my $max_group_id;

    my $outermost_pair_count = @{$stem_outermost_pairs};

    for (my $i = 0; $i < $outermost_pair_count; $i++) {
	my ($curr_pair_upstream_pos, $curr_pair_downstream_pos) = @{$stem_outermost_pairs->[$i]};
	my $curr_pair_group_id = $paired_pos_to_group_id->{$curr_pair_upstream_pos};

	my $succ_pair_crossings = [];

	for (my $j = $i + 1; $j < $outermost_pair_count; $j++) {
	    my ($candidate_pair_upstream_pos, $candidate_pair_downstream_pos) = @{$stem_outermost_pairs->[$j]};
	    if ($candidate_pair_upstream_pos > $curr_pair_downstream_pos) {
		last;
	    }

	    if ($candidate_pair_downstream_pos > $curr_pair_downstream_pos) {
		my $crossing_pair_group_id = $paired_pos_to_group_id->{$candidate_pair_upstream_pos};
		if (defined($curr_pair_group_id)) {
		    if (!defined($crossing_pair_group_id)) {
			push @{$knotted_pair_pos_groups->[$curr_pair_group_id]}, $candidate_pair_upstream_pos;
			push @{$knotted_pair_pos_groups->[$curr_pair_group_id]}, $candidate_pair_downstream_pos;
			$paired_pos_to_group_id->{$candidate_pair_upstream_pos} = $curr_pair_group_id;
		    }
		    elsif ($crossing_pair_group_id != $curr_pair_group_id) {
			my @merged_pos_group = (@{$knotted_pair_pos_groups->[$curr_pair_group_id]}, @{$knotted_pair_pos_groups->[$crossing_pair_group_id]});
			$knotted_pair_pos_groups->[$curr_pair_group_id] = \@merged_pos_group;

			foreach (@{$knotted_pair_pos_groups->[$crossing_pair_group_id]}) {
			    if (exists($paired_pos_to_group_id->{$_})) {
				$paired_pos_to_group_id->{$_} = $curr_pair_group_id;
			    }
			}

			delete $knotted_pair_pos_groups->[$crossing_pair_group_id];
		    }
		}
		else {
		    if (defined($crossing_pair_group_id)) {
			$curr_pair_group_id = $crossing_pair_group_id;
			push @{$knotted_pair_pos_groups->[$curr_pair_group_id]}, $curr_pair_upstream_pos;
			push @{$knotted_pair_pos_groups->[$curr_pair_group_id]}, $curr_pair_downstream_pos;
		    }
		    else {
			$curr_pair_group_id = $max_group_id++;
			$knotted_pair_pos_groups->[$curr_pair_group_id] = [$curr_pair_upstream_pos, $curr_pair_downstream_pos, $candidate_pair_upstream_pos, $candidate_pair_downstream_pos];
			$paired_pos_to_group_id->{$candidate_pair_upstream_pos} = $curr_pair_group_id;
		    }
		}

		push @{$succ_pair_crossings}, $candidate_pair_upstream_pos;
	    }
	}

	$outermost_pair_crossings->{$curr_pair_upstream_pos} = $succ_pair_crossings;
    }

    return ($knotted_pair_pos_groups, $outermost_pair_crossings);
}

sub _get_prim_pseudoknots {
    my ($stems, $knotted_pair_pos_groups, $outermost_pair_crossings, $paired_pos_ptrs) = @_;

    my $primitive_pseudoknots = [];

    for (my $i = 0; $i < @{$knotted_pair_pos_groups}; $i++) {
	if (!defined($knotted_pair_pos_groups->[$i])) {
	    next;
	}

	my @sorted_knot_pair_pos = sort {$a <=> $b} @{$knotted_pair_pos_groups->[$i]};
	my $prev_knot_pair_pos = {};
	for (my $j = 1; $j < @sorted_knot_pair_pos; $j++) {
	    $prev_knot_pair_pos->{$sorted_knot_pair_pos[$j]} = $sorted_knot_pair_pos[$j - 1];
	}

	my ($prim_pseudoknot_stems, $prim_pseudoknot_stem) = ([], []);
	my $knot_pair_pos_to_stem_id = {};
	my $max_stem_id = 0;

	for (my $j = 0; $j < (@sorted_knot_pair_pos - 1); $j++) {
	    my $curr_pos = $sorted_knot_pair_pos[$j];
	    my $curr_paired_pos = $paired_pos_ptrs->[$curr_pos];
	    if ($curr_pos > $curr_paired_pos) {
		next;
	    }

	    my @merged_stem = (@{$prim_pseudoknot_stem}, @{$stems->{$curr_pos}});
	    $prim_pseudoknot_stem = \@merged_stem;

	    my $next_pos = $sorted_knot_pair_pos[$j + 1];
	    my $next_paired_pos = $paired_pos_ptrs->[$next_pos];
	    if ($prev_knot_pair_pos->{$curr_paired_pos} != $next_paired_pos) {
		push @{$prim_pseudoknot_stems}, $prim_pseudoknot_stem;
		$knot_pair_pos_to_stem_id->{$curr_pos} = $max_stem_id++;
		$prim_pseudoknot_stem = [];
	    }
	}

	my $prim_pseudoknot_stem_crossings = [];
	while (my ($knot_pair_upstream_pos, $stem_id) = each %{$knot_pair_pos_to_stem_id}) {
	    my $stem_crossings = [];
	    my $knot_pair_crossings = $outermost_pair_crossings->{$knot_pair_upstream_pos};
	    foreach (@{$knot_pair_crossings}) {
		if (exists($knot_pair_pos_to_stem_id->{$_})) {
		    push @{$stem_crossings}, $knot_pair_pos_to_stem_id->{$_};
		}
	    }

	    $prim_pseudoknot_stem_crossings->[$stem_id] = $stem_crossings;
	}

	push @{$primitive_pseudoknots}, [$prim_pseudoknot_stems, $prim_pseudoknot_stem_crossings];
    }

    return $primitive_pseudoknots;
}

1;
