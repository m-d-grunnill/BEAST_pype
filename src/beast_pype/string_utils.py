

def rstrip_after_nth(string, chars, n):
    """
    Strip everything after the nth occurance of a substring (includes nth occurance).
    """
    return chars.join(string.split(chars)[:n])


    