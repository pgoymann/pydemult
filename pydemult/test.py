import Levenshtein

def find_best_match_shift(TAG_seq, tags, maximum_distance):
    """
    Find the best match from the list of tags with sliding window.
    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R1 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.
    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_match = 'unmapped'
    best_score = maximum_distance
    shifts = range(0,len(TAG_seq) - len(max(tags,key=len)))

    for shift in shifts:
        for tag, name in tags.items():
            score = Levenshtein.hamming(tag, TAG_seq[shift:len(tag)+shift])
            if score == 0:
                #Best possible match
                return(name)
            elif score <= best_score:
                best_score = score
                best_match = name
                return(best_match)
    return(best_match)
print(find_best_match_shift('AAAAAAAAAAAAA', {'AAAAAAAAAAAA':1,'CCCCCCCC':2}, 2))

