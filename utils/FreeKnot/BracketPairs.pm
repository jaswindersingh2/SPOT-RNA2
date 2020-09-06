#Bracket handler for DPParser

package BracketPairs;
use strict;

my $open_bracket_map = {")" => "(", "]" => "[", "}" => "{", ">" => "<"};

#Check whether a symbol (in dot-parentheses format) is an open bracket
sub is_open_bracket {
    my (undef, $symbol) = @_;

    if ($symbol =~ /^[\(\[{<A-Z]$/) {
	return 1;
    }

    return 0;
}

#Return a corresponding close bracket for an open bracket input  
sub get_open_bracket {
    my (undef, $close_bracket) = @_;

    if ($close_bracket =~ /^[\)\]}>]$/) {
	return $open_bracket_map->{$close_bracket};
    }
    elsif ($close_bracket =~ /^[a-z]$/) {
	return uc $close_bracket;
    }
    else {
	die "Unknown closing bracket\n";
    }
}

1;
