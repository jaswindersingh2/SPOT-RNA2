#The MIS algorithm module. It is an extension of the k-MIS algorithm proposed by Byskov (Byskov, J., 2004)

package MIS;

use strict;

use constant D => 3;

my $miss;
my $checked_sets;

sub get_mis {
    my (undef, $circle_graph) = @_;

    $miss = [];
    $checked_sets = {};

    #Initialize the vertex set with goal opposing vertices filtered
    my $vertex_set = VertexSubset->new($circle_graph);
    #Call the branching algorithm _search_mis, all the MWISs will be stored in $mwiss
    _search_mis($vertex_set, [], $circle_graph);

    undef $checked_sets;

    return $miss;
}

sub _search_mis {
    my ($vertex_subset, $candidate_set, $circle_graph) = @_;

    if ($vertex_subset->get_size() == 0) {
	#If the vertex subset is empty, check whether the $candidate_set is an independent set. If so then it is
	#an MIS and the toal vertex weight is evaluated. Those with the best overall weight (according to the
	#goal specified by $criteria) are put in $miss. Since the same subset may appear more than once,
	#$checked_sets stores all the subset verified before to avoid unnecessary checking.
	@{$candidate_set} = sort {$a <=> $b} @{$candidate_set};
	my $candidate_set_id = join('-', @{$candidate_set});
	if (!exists($checked_sets->{$candidate_set_id}) && _is_independent_set($candidate_set, $circle_graph)) {
	    push @{$miss}, $candidate_set;
	    $checked_sets->{$candidate_set_id} = $candidate_set;
	}
    }
    else {
	my ($highest_degree_vertices, $highest_vertex_degree) = $vertex_subset->get_highest_degree_vertex_info();
	#If the highest vertex degree is at least D, select a vertex with such degree to branch
	if ($highest_vertex_degree >= D) {
	    my @self_adj_vertices = (@{$vertex_subset->get_adjacent_vertices_at($highest_degree_vertices->[0])}, $highest_degree_vertices->[0]);
	    my @expanded_candidate_set = (@{$candidate_set}, $highest_degree_vertices->[0]);
	    #Branch on by including the selected vertex in $candidate_set
	    _search_mis($vertex_subset->get_subset(\@self_adj_vertices), \@expanded_candidate_set, $circle_graph);

	    #Branch on by just excluding the selected vertex in $candidate_set
	    _search_mis($vertex_subset->get_subset([$highest_degree_vertices->[0]]), $candidate_set, $circle_graph);
	}
	#If the highest vertex degree is lower than D, select a vertex with the lowest vertex degree to branch instead
	else {
	    my ($lowest_degree_vertices, undef) = $vertex_subset->get_lowest_degree_vertex_info();
	    my $adj_vertices = $vertex_subset->get_adjacent_vertices_at($lowest_degree_vertices->[0]);
	    my @self_adj_vertices1 = (@{$adj_vertices}, $lowest_degree_vertices->[0]);
	    my @expanded_candidate_set1 = (@{$candidate_set}, $lowest_degree_vertices->[0]);
	    #Branch on by including the selected vertex in $candidate_set
	    _search_mis($vertex_subset->get_subset(\@self_adj_vertices1), \@expanded_candidate_set1, $circle_graph);

	    #Branch on by enumerating and including each adjacent vertex of the selected vertex in $candidate_set
	    foreach (@{$adj_vertices}) {
		my @expanded_candidate_set2 = (@{$candidate_set}, $_);
		my @self_adj_vertices2 = (@{$vertex_subset->get_adjacent_vertices_at($_)}, $_);
		_search_mis($vertex_subset->get_subset(\@self_adj_vertices2), \@expanded_candidate_set2, $circle_graph);
	    }
	}
    }
}

sub _is_independent_set {
    my ($candidate_set, $circle_graph) = @_;

    my ($all_non_adj_vertex_mask, $candidate_set_bitstrings) = ([], []);

    for (my $i = @{$candidate_set} - 1; $i >= 0; $i--) {
	my $non_adj_vertex_mask = $circle_graph->get_non_adj_vertex_mask_at($candidate_set->[$i]);
	for (my $j = 0; $j < @{$candidate_set_bitstrings}; $j++) {
	    if (($candidate_set_bitstrings->[$j] & $non_adj_vertex_mask->[$j]) != $candidate_set_bitstrings->[$j]) {
		return 0;
	    }
	}

	my ($vertex_bitstring_segment_num, $vertex_bitstring) = @{$circle_graph->get_vertex_bitstring_segment_at($candidate_set->[$i])};
	$candidate_set_bitstrings->[$vertex_bitstring_segment_num] = $candidate_set_bitstrings->[$vertex_bitstring_segment_num] | $vertex_bitstring;
    }

    return 1;
}

1;
