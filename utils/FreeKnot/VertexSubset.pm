#Module that represents the vertex subset in the MWIS algorithm. All the vertices of the knot-stem
#graph are added to this subset (with the goal opposing vertices filtered) at initialization. When
#the MWIS algorithm proceeds, vertices are gradually removed from this subset and the algorithm
#stops when this subset is empty.
#
#This subset also keeps the adjacent vertices for each vertex in it, as well as the vertex degrees.
#It enables the MWIS algorithm to select the highest degree and lowest degree vertices, and to
#further generate a new subset of it while updating the adjacent vertices and vertex degrees.

package VertexSubset;

use strict;

sub new {
#    my (undef, $circle_graph, $stem_scores, $criteria) = @_;
    my (undef, $circle_graph) = @_;

    my ($vertex_degrees, $adj_vertex_sets) = ({}, {});

    my $subset_size = 0;

    for (my $i = $circle_graph->get_vertex_count() - 1; $i >= 0; $i--) {
	$vertex_degrees->{$i} = 0;

	foreach (@{$circle_graph->get_edges_at($i)}) {
	    $vertex_degrees->{$i}++;
	    $vertex_degrees->{$_}++;
	    $adj_vertex_sets->{$i}{$_} = 1;
	    $adj_vertex_sets->{$_}{$i} = 1;
	}

	$subset_size++;
    }

    my ($highest_degree_vertices, $lowest_degree_vertices, $highest_vertex_degree, $lowest_vertex_degree) = _get_highest_and_lowest_degree_vertices($vertex_degrees);

    my $self = {};
    $self->{subset_size} = $subset_size;
    $self->{vertex_degrees} = $vertex_degrees;
    $self->{adj_vertex_sets} = $adj_vertex_sets;
    $self->{highest_degree_vertices} = $highest_degree_vertices;
    $self->{lowest_degree_vertices} = $lowest_degree_vertices;
    $self->{highest_vertex_degree} = $highest_vertex_degree;
    $self->{lowest_vertex_degree} = $lowest_vertex_degree;

    bless $self;

    return $self;
}

#Generate a new subset instance by removing the vertices specified in the input
sub get_subset {
    my ($self, $vertices_to_remove) = @_;

    my $subset_size = 0;
    my ($subset_vertex_degrees, $subset_adj_vertex_sets) = ({}, {});

    my %delete_vertices = map {$_ => 1} @{$vertices_to_remove};
    my $vertex_degrees = $self->{vertex_degrees};
    foreach (keys %{$vertex_degrees}) {
	if (!exists($delete_vertices{$_})) {
	    $subset_vertex_degrees->{$_} = 0;
	    $subset_adj_vertex_sets->{$_} = {};
	    $subset_size++;
	}
    }

    my $adj_vertex_sets = $self->{adj_vertex_sets};
    while (my ($vertex, $adj_vertices) = each %{$adj_vertex_sets}) {
	if (!exists($delete_vertices{$vertex})) {
	    foreach (keys %{$adj_vertices}) {
		if ($vertex < $_ && !exists($delete_vertices{$_})) {
		    $subset_adj_vertex_sets->{$vertex}{$_} = 1;
		    $subset_adj_vertex_sets->{$_}{$vertex} = 1;
		    $subset_vertex_degrees->{$vertex}++;
		    $subset_vertex_degrees->{$_}++;
		}
	    }
	}
    }

    my ($highest_degree_vertices, $lowest_degree_vertices, $highest_vertex_degree, $lowest_vertex_degree) = _get_highest_and_lowest_degree_vertices($subset_vertex_degrees);

    my $subset_self = {};
    $subset_self->{subset_size} = $subset_size;
    $subset_self->{vertex_degrees} = $subset_vertex_degrees;
    $subset_self->{adj_vertex_sets} = $subset_adj_vertex_sets;
    $subset_self->{highest_degree_vertices} = $highest_degree_vertices;
    $subset_self->{lowest_degree_vertices} = $lowest_degree_vertices;
    $subset_self->{highest_vertex_degree} = $highest_vertex_degree;
    $subset_self->{lowest_vertex_degree} = $lowest_vertex_degree;

    bless $subset_self;

    return $subset_self;
}

sub _get_highest_and_lowest_degree_vertices {
    my $vertex_degrees = shift;

    my ($highest_degree_vertices, $lowest_degree_vertices) = ([], []);
    my ($highest_vertex_degree, $lowest_vertex_degree) = (-1, -1);

    while (my ($vertex, $vertex_degree) = each %{$vertex_degrees}) {
	if ($vertex_degree > $highest_vertex_degree) {
	    $highest_degree_vertices = [$vertex];
	    $highest_vertex_degree = $vertex_degree;
	}
	elsif ($vertex_degree == $highest_vertex_degree) {
	    push @{$highest_degree_vertices}, $vertex;
	}

	if ($vertex_degree < $lowest_vertex_degree || $lowest_vertex_degree < 0) {
	    $lowest_degree_vertices = [$vertex];
	    $lowest_vertex_degree = $vertex_degree;
	}
	elsif ($vertex_degree == $lowest_vertex_degree) {
	    push @{$lowest_degree_vertices}, $vertex;
	}
    }

    my @sorted_highest_degree_vertices = sort {$a <=> $b} @{$highest_degree_vertices};
    my @sorted_lowest_degree_vertices = sort {$a <=> $b} @{$lowest_degree_vertices};

    return \@sorted_highest_degree_vertices, \@sorted_lowest_degree_vertices, $highest_vertex_degree, $lowest_vertex_degree;
}

sub get_size {
    my $self = shift;

    return $self->{subset_size};
}

sub get_vertices {
    my $self = shift;

    my @vertices = sort {$a <=> $b} keys %{$self->{vertex_degrees}};

    return \@vertices;
}

sub get_adjacent_vertices_at {
    my ($self, $vertex) = @_;

    my $adj_vertex_sets = $self->{adj_vertex_sets};
    if (exists($adj_vertex_sets->{$vertex})) {
	my @adj_vertices = sort {$a <=> $b} keys %{$adj_vertex_sets->{$vertex}};
	return \@adj_vertices;
    }

    return [];
}

sub get_highest_degree_vertex_info {
    my $self = shift;

    return $self->{highest_degree_vertices}, $self->{highest_vertex_degree};
}

sub get_lowest_degree_vertex_info {
    my $self = shift;

    return $self->{lowest_degree_vertices}, $self->{lowest_vertex_degree};
}

1;
