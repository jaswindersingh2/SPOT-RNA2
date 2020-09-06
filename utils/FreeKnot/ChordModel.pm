#Chord model of the circle graph representing a primitive pseudoknot. Each chord denotes a unique
#crossing base pair stem in the primitive pseudoknot. If two stems cross, then their corresponding
#chords also cross. Each chord is associated with its underlying base pairs.

package ChordModel;

use strict;

sub new {
    my (undef, $primitive_pseudoknot) = @_;

    my $prim_pseudoknot_stems = $primitive_pseudoknot->[0];
    my $chord_end_point_num_map = _get_chord_end_point_num_map($prim_pseudoknot_stems);

    my ($chord_edges, $all_chord_base_pairs) = ({}, {});
    my ($chord_end_point_nums, $end_point_to_edge_map, $is_left_end_points) = ([], [], []);

    foreach (@{$prim_pseudoknot_stems}) {
	my $chord_left_end_point_num = $chord_end_point_num_map->{$_->[0][0]};
	my $chord_right_end_point_num = $chord_end_point_num_map->{$_->[0][1]};
	push @{$chord_end_point_nums}, ($chord_left_end_point_num, $chord_right_end_point_num);
	$is_left_end_points->[$chord_left_end_point_num] = 1;
	$is_left_end_points->[$chord_right_end_point_num] = 0;
	$all_chord_base_pairs->{$chord_left_end_point_num . '-' . $chord_right_end_point_num} = $_;

	my $chord_edge = [$chord_left_end_point_num, $chord_right_end_point_num];
	$chord_edges->{$chord_left_end_point_num . '-' . $chord_right_end_point_num} = $chord_edge;
	$end_point_to_edge_map->[$chord_left_end_point_num] = $chord_edge;
	$end_point_to_edge_map->[$chord_right_end_point_num] = $chord_edge;
    }

    my @sorted_chord_end_point_nums = sort {$b <=> $a} @{$chord_end_point_nums};

    my $self = {};
    $self->{chord_end_point_nums} = \@sorted_chord_end_point_nums;
    $self->{chord_edges} = $chord_edges;
    $self->{end_point_to_edge_map} = $end_point_to_edge_map;
    $self->{is_left_end_points} = $is_left_end_points;
    $self->{all_chord_base_pairs} = $all_chord_base_pairs;

    bless $self;

    return $self;
}

sub _get_chord_end_point_num_map {
    my $prim_pseudoknot_stems = shift;

    my $stem_end_points = [];

    foreach (@{$prim_pseudoknot_stems}) {
	push @{$stem_end_points}, $_->[0][0];
	push @{$stem_end_points}, $_->[0][1];
    }

    my @sorted_stem_end_points = sort {$a <=> $b} @{$stem_end_points};

    my $chord_end_point_num_map = {};
    for (my $i = 0; $i < @sorted_stem_end_points; $i++) {
	$chord_end_point_num_map->{$sorted_stem_end_points[$i]} = $i + 1;
    }

    return $chord_end_point_num_map;
}

sub get_chord_end_point_nums {
    my $self = shift;

    return $self->{chord_end_point_nums};
}

sub get_chord_edges {
    my $self = shift;

    return $self->{chord_edges};
}

sub get_chord_edge_count {
    my $self = shift;

    return scalar(keys %{$self->{chord_edges}});
}

sub get_chord_edge_by_end_point {
    my ($self, $end_point_num) = @_;

    my $end_point_to_edge_map = $self->{end_point_to_edge_map};

    return $end_point_to_edge_map->[$end_point_num];
}

sub is_left_end_point {
    my ($self, $end_point_num) = @_;

    my $is_left_end_points = $self->{is_left_end_points};

    return $is_left_end_points->[$end_point_num];
}

sub get_chord_base_pairs {
    my ($self, $chord_left_end_point, $chord_right_end_point) = @_;

    my $all_chord_base_pairs = $self->{all_chord_base_pairs};

    return $all_chord_base_pairs->{$chord_left_end_point . '-' . $chord_right_end_point};
}

1;
