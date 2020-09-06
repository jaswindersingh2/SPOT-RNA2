package ScoringFunctions;

use strict;

#my $free_energy_params;
#my $canonical_base_pairs = {'AU' => 0, 'CG' => 0, 'GC' => 0, 'UA' => 0, 'GU' => 0, 'UG' => 0};

#Return a scoring function according to the choice selected
sub get_scoring_function {
    my (undef, $option) = @_;

    if ($option eq 'bp') {
	return \&_base_pair_score, 'max', 0;
    }
    elsif ($option eq 'stem') {
	return \&_stem_score, 'max', 0;
    }
#    elsif ($option eq 'sstab') {
#	$free_energy_params = _init_free_energy_parameters();
#	return \&_stem_bp_stability, 'min', 1;
#    }
    elsif ($option eq 'hb') {
	return \&_hydrogen_bond, 'max', 0;
    }
    elsif ($option eq 'fe') {
	return \&_overall_stability, 'min', 1;
    }
    else {
	return undef, undef, undef;
    }
}

#Number of base pairs in a stem as the stem score
sub _base_pair_score {
    my $chord_attrs = shift;

    my $stem_pair_count = $chord_attrs->{pair_count};
    if (defined($stem_pair_count)) {
	return $stem_pair_count;
    }

    return 0;
}

#Each stem scores equally as 1
sub _stem_score {
    return 1;
}

#GC and CG bonds = 3, other canonical of GU pairs = 2
sub _hydrogen_bond {
    my ($chord_attrs, $base_seq) = @_;

    my $stem_base_pairs = $chord_attrs->{base_pairs};
    my $total_score = 0;

    foreach (@{$stem_base_pairs}) {
	my $base_pair_type = uc($base_seq->[$_->[0] - 1] . $base_seq->[$_->[1] - 1]);
	if ($base_pair_type eq 'GC' || $base_pair_type eq 'CG') {
	    $total_score += 3;
	}
	elsif ($base_pair_type eq 'AU' || $base_pair_type eq 'UA' ||
	       $base_pair_type eq 'GU' || $base_pair_type eq 'UG') {
	    $total_score += 2;
	}
    }

    return $total_score;
}

#This allows all MISs to be reported as MWISs and they will be converted to all possible
#de-knotted structures to determine the minimum free energy (MFE)
sub _overall_stability {
    return 0;
}

1;
