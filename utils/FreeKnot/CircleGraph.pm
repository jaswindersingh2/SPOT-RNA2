#Circle graph is the graphical model for a primitive pseudoknot. Each vertex represents a crossing stem of
#the pseudoknot, and each edge represents a crossing between two stems. The vertex attributes store information
#such as number of base pairs, paired positions.
#
#Every vertex is represented by a unique bitstring and its adjacent vertices are represented by a bitstring
#mask. The least significant bit (LSB) represents the most preceding vertex and the most significant bit (MSB)
#represents the least preceding vertex of the knot-stem graph. However, the no. of vertices may exceed the
#length of one bitstring. To solve this problem, multiple bitstrings are required to form a bitstring long
#enough for each bit position to uniquely identify a vertex. This 'long' bitstring is disassembled into an
#array of bitstrings and each array element is called a bitstring segment. Every bit position of the 'long'
#bitstring is then transformed by a (segment no., segment bitstring) pair.

package CircleGraph;

use strict;

sub new {
    my (undef, $primitive_pseudoknot, $os_bit) = @_;

    my $vertex_attrs = [];
#    my ($stem_pair_counts, $gains) = ([], []);
    my ($vertex_bitstring_segments, $non_adj_vertex_masks) = ([], []);

    my ($prim_pseudoknot_stems, $prim_pseudoknot_stem_crossings) = @{$primitive_pseudoknot};
    my $vertex_count = @{$prim_pseudoknot_stems};

    my ($bitstring_segment_num, $vertex_bit) = (0, 0);

    for (my $i = $vertex_count - 1; $i >= 0; $i--) {
	my $prim_pseudoknot_stem = $prim_pseudoknot_stems->[$i];
#	$stem_pair_counts->[$i] = @{$prim_pseudoknot_stem};
#	$gains->[$i] = $stem_pair_counts->[$i];

	$vertex_bitstring_segments->[$i] = [$bitstring_segment_num, 1 << $vertex_bit];
	my $non_adj_vertex_mask_bitstrings = [];

	my $stem_crossings = $prim_pseudoknot_stem_crossings->[$i];
	my $next_crossing_index = @{$stem_crossings} - 1;
	my $next_crossing_stem_id;
	if ($next_crossing_index >= 0) {
	    $next_crossing_stem_id = $stem_crossings->[$next_crossing_index];
	}

	for (my $j = $vertex_count - 1; $j > $i; $j--) {
	    if ($next_crossing_index >= 0 && $j == $next_crossing_stem_id) {
#		$gains->[$i] -= $stem_pair_counts->[$j];
#		$gains->[$j] -= $stem_pair_counts->[$i];
		if (--$next_crossing_index >= 0) {
		    $next_crossing_stem_id = $stem_crossings->[$next_crossing_index];
		}
	    }
	    else {
		my $non_adj_vertex_bitstring_segment_num = $vertex_bitstring_segments->[$j][0];
		$non_adj_vertex_mask_bitstrings->[$non_adj_vertex_bitstring_segment_num] = $non_adj_vertex_mask_bitstrings->[$non_adj_vertex_bitstring_segment_num] | $vertex_bitstring_segments->[$j][1];
	    }
	}
	
	$non_adj_vertex_masks->[$i] = $non_adj_vertex_mask_bitstrings;

	if (++$vertex_bit == $os_bit) {
	    $bitstring_segment_num++;
	    $vertex_bit = 0;
	}
    }

    for (my $i = 0; $i < $vertex_count; $i++) {
	my $attrs = {};
#	$attrs->{pair_count} = $stem_pair_counts->[$i];
#	$attrs->{gain} = $gains->[$i];
	$attrs->{stem_pairs} = $prim_pseudoknot_stems->[$i];
	$vertex_attrs->[$i] = $attrs;
    }

    my $self = {};
    $self->{vertex_count} = $vertex_count;
    $self->{vertex_attrs} = $vertex_attrs;
    $self->{edges} = $prim_pseudoknot_stem_crossings;
    $self->{vertex_bitstring_segments} = $vertex_bitstring_segments;
    $self->{non_adj_vertex_masks} = $non_adj_vertex_masks;

    bless $self;

    return $self;
}

sub get_vertex_count {
    my $self = shift;

    return $self->{vertex_count};
}

sub get_vertex_attrs_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return [];
    }

    my $vertex_attrs = $self->{vertex_attrs};

    return $vertex_attrs->[$vertex_num];
}

sub get_edges_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return [];
    }

    my $edges = $self->{edges};

    return $edges->[$vertex_num];
}

#Return the bitstring segment of the vertex. Each bitstring segment is a (segment no.,
#segment bitstring) pair.
sub get_vertex_bitstring_segment_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return [];
    }

    my $bitstring_segments = $self->{vertex_bitstring_segments};

    return $bitstring_segments->[$vertex_num];
}

#Returns bitstring segments that filter all the subsequent adjacent vertices
sub get_non_adj_vertex_mask_at {
    my ($self, $vertex_num) = @_;

    if ($vertex_num >= $self->{vertex_count}) {
	return 0;
    }

    my $non_adj_vertex_mask_bitstrings = $self->{non_adj_vertex_masks};

    return $non_adj_vertex_mask_bitstrings->[$vertex_num];
}

1;
