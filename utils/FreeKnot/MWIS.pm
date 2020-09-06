#Modified circle graph MWIS algorithm based on that proposed by Valiente (Valiente, G., 2003), with
#enhancement suggested by Nash et al. (Nash, N., Lelait, S., and Gregg, D., 2009). It operates with
#the chord model and reports either single solution or all solutions according to the user option.

package MWIS;

use strict;

sub get_mwis {
    my (undef, $chord_model, $base_seq, $scoring_function, $criteria, $is_report_all) = @_;

    my $chord_weights = _get_chord_weights($chord_model, $base_seq, $scoring_function);

    my $end_point_count = $chord_model->get_chord_edge_count() * 2;

    #Enhancement by Nash et al. to get MWISs (in variable c) and the scores (in variable cmis) in
    #every region bounded by the endpoints of each chord.
    my ($m, $p) = ([], []);
    my ($cmis, $c) = ({}, {});

    for (my $i = 1; $i <= $end_point_count + 1; $i++) {
	$m->[$i] = 0;
	$p->[$i] = [0];
    }

    my $last = 1;

    for (my $i = 1; $i <= $end_point_count; $i++) {
	if ($chord_model->is_left_end_point($i)) {
	    next;
	}

	my ($left_end_point, $right_end_point) = @{$chord_model->get_chord_edge_by_end_point($i)};

	for (my $j = $last; $j > $left_end_point; $j--) {
	    $m->[$j] = $m->[$j + 1];
	    $p->[$j] = $p->[$j + 1];

	    if ($chord_model->is_left_end_point($j)) {
		my (undef, $inner_right_end_point) = @{$chord_model->get_chord_edge_by_end_point($j)};
		my $candidate_m = $m->[$inner_right_end_point + 1] + $cmis->{$j . '-' . $inner_right_end_point};

		if (($criteria eq 'max' && $candidate_m > $m->[$j]) ||
		    ($criteria eq 'min' && $candidate_m < $m->[$j])) {
		    $m->[$j] = $candidate_m;
		    $p->[$j] = [$inner_right_end_point];
		}
		elsif ($is_report_all && $candidate_m == $m->[$j]) {
		    my @arr_clone = @{$p->[$j + 1]};
		    $p->[$j] = [$inner_right_end_point];
		    push @{$p->[$j]}, @arr_clone;
		}
	    }
	}

	$cmis->{$left_end_point . '-' . $right_end_point} = $m->[$left_end_point + 1] + $chord_weights->{$left_end_point . '-' . $right_end_point};
	$c->{$left_end_point . '-' . $right_end_point} = _add_front($p, $left_end_point + 1, $chord_model, []);
	$last = $left_end_point;
    }

    #Algorithm proposed by Valiente to obtain MWISs starting at each endpoint. Only those chords
    #in the MWIS that are not bounded by other chords in the same MWIS set are stored.
    my ($t_structures, $t_struct_weights) = ([], []);

    foreach (@{$chord_model->get_chord_end_point_nums()}) {
	$t_structures->[$_] = [[]];

	if (!$chord_model->is_left_end_point($_)) {
	    if ($_ < $end_point_count) {
		@{$t_structures->[$_]} = @{$t_structures->[$_ + 1]};
		$t_struct_weights->[$_] = $t_struct_weights->[$_ + 1];
	    }
	    else {
		$t_struct_weights->[$_] = 0;
	    }
	}
	else {
	    my $chord_edge = $chord_model->get_chord_edge_by_end_point($_);
	    my $candidate_total_chord_weight = $cmis->{$chord_edge->[0] . '-' . $chord_edge->[1]};

	    if ($chord_edge->[1] < $end_point_count) {
		$candidate_total_chord_weight += $t_struct_weights->[$chord_edge->[1] + 1];
	    }

	    if (($criteria eq 'max' && $candidate_total_chord_weight > $t_struct_weights->[$_ + 1]) ||
		($criteria eq 'min' && $candidate_total_chord_weight < $t_struct_weights->[$_ + 1]) ||
		($candidate_total_chord_weight == $t_struct_weights->[$_ + 1] && $is_report_all)) {
		my $generated_new_t_structures;

		if ($candidate_total_chord_weight == $t_struct_weights->[$_ + 1]) {
		    @{$generated_new_t_structures} = @{$t_structures->[$_ + 1]};
		}
		else {
		    $generated_new_t_structures = [];
		}

		if ($chord_edge->[1] < $end_point_count) {
		    foreach my $t_structure (@{$t_structures->[$chord_edge->[1] + 1]}) {
			my @new_t_structure = @{$t_structure};
			unshift @new_t_structure, $chord_edge;
			push @{$generated_new_t_structures}, \@new_t_structure;
		    }
		}
		else {
		    push @{$generated_new_t_structures}, [$chord_edge];
		}

		$t_structures->[$_] = $generated_new_t_structures;
		$t_struct_weights->[$_] = $candidate_total_chord_weight;
	    }
	    else {
		$t_structures->[$_] = $t_structures->[$_ + 1];
		$t_struct_weights->[$_] = $t_struct_weights->[$_ + 1];
	    }
	}
    }

    my $mwiss = _restore_chord_mwiss($t_structures->[1], $c);

    return $mwiss;
}

#Generate all the MWISs in the region bounded by the endpoints of a single chord
sub _add_front {
    my ($p, $start_pos, $chord_model, $org_c_element) = @_;

    my $p_element = $p->[$start_pos];

    if ($p_element->[0] > 0) {
	my $new_c_element = [];

	foreach (@{$p_element}) {
	    my $chord_edge = $chord_model->get_chord_edge_by_end_point($_);
	    my $expanded_c_element = [];

	    if (!defined($org_c_element->[0])) {
		push @{$expanded_c_element}, [$chord_edge];
	    }
	    else {
		foreach my $element_value (@{$org_c_element}) {
		    my @arr_clone = @{$element_value};
		    push @arr_clone, $chord_edge;
		    push @{$expanded_c_element}, \@arr_clone;
		}
	    }

	    my $new_values = _add_front($p, $_, $chord_model, $expanded_c_element);
	    push @{$new_c_element}, @{$new_values};
	}

	return $new_c_element;
    }

    return $org_c_element;
}

sub _get_chord_weights {
    my ($chord_model, $base_seq, $scoring_function) = @_;

    my $chord_weights = {};

    foreach (values %{$chord_model->get_chord_edges()}) {
	my $chord_base_pairs = $chord_model->get_chord_base_pairs($_->[0], $_->[1]);
	my $chord_attrs = {};
	$chord_attrs->{base_pairs} = $chord_base_pairs;
	$chord_attrs->{pair_count} = @{$chord_base_pairs};
	$chord_weights->{$_->[0] . '-' . $_->[1]} = $scoring_function->($chord_attrs, $base_seq);
    }

    return $chord_weights;
}

#Recover the MWISs from the chord sets in variable c
sub _restore_chord_mwiss {
    my ($chord_edge_sets, $c) = @_;

    my $chord_mwiss = [];

    foreach my $chord_edge_set (@{$chord_edge_sets}) {
	my $single_chord_edge_set_mwiss = [$chord_edge_set];

	foreach my $chord_edge (@{$chord_edge_set}) {
	    my $inner_chord_edge_sets = $c->{$chord_edge->[0] . '-' . $chord_edge->[1]};
	    if (!defined($inner_chord_edge_sets->[0])) {
		next;
	    }

	    my $inner_chord_mwiss = _restore_chord_mwiss($inner_chord_edge_sets, $c);
	    my @org_single_chord_edge_set_mwiss = @{$single_chord_edge_set_mwiss};
	    $single_chord_edge_set_mwiss = [];

	    foreach my $single_chord_edge_set_mwis (@org_single_chord_edge_set_mwiss) {
		foreach my $inner_chord_mwis (@{$inner_chord_mwiss}) {
		    my @merged_mwis = (@{$single_chord_edge_set_mwis}, @{$inner_chord_mwis});
		    push @{$single_chord_edge_set_mwiss}, \@merged_mwis;
		}
	    }
	}

	push @{$chord_mwiss}, @{$single_chord_edge_set_mwiss};

    }

    return $chord_mwiss;
}

1;
